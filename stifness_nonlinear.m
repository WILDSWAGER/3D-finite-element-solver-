function FEM_Solver_Nonlinear_v10()
% =========================================================================
%  FEM_Solver_Nonlinear_v10.m
%  Nonlinear Static FEM Solver — Tetrahedral Elements (T4)
%
%  Run by typing:  FEM_Solver_Nonlinear_v10()  at the MATLAB prompt.
%
%  MODELS:
%    1 -> Neo-Hookean
%    2 -> Mooney-Rivlin
%    3 -> Yeoh
%    4 -> EESM Arruda-Boyce
%    5 -> Phenomenological  Etan=E0*(1+eps+eps/alpha)
%
%  BUG HISTORY:
%  v6 : lam=1+eps models 1-3,5; F_int element-by-element; model 5 added
%  v7 : EESM beta denominator (1-lamr^3); compression lam convention
%  v8 : EESM secant modulus Et=sigma/eps  (tangent Et<0 -> diverge)
%  v9 : Wrap in function() to fix "misplaced function" error;
%       EESM overflow guard (beta>700); phenomenological fallback
%       when eps>eps_max_safe; NR residual checked after update
%  v10: min NR iterations enforced (>=5 recommended);
%       WT near-cancellation warning printed in EESM diagnostic;
%       physical displacement sanity check after each NR step;
%       max_iter<5 triggers a clear user warning
% =========================================================================

fprintf('=============================================================\n');
fprintf('   NONLINEAR FEM SOLVER  (T4, v10)                           \n');
fprintf('=============================================================\n\n');

% =========================================================================
%  SECTION 1 — USER INPUTS
% =========================================================================

fprintf('--- UNIT SYSTEM ---\n  1->mm  2->cm  3->m\n');
unit_choice = input('  STL unit (1/2/3): ');
while ~ismember(unit_choice,[1 2 3]); unit_choice=input('  Enter 1,2 or 3: '); end
switch unit_choice
    case 1; q=3; unit_label='mm';
    case 2; q=2; unit_label='cm';
    case 3; q=0; unit_label='m';
end
fprintf('  Unit: %s\n\n',unit_label);

fprintf('--- STL GEOMETRY ---\n');
stl_filename = input('  STL filename (e.g. cube3cm.stl): ','s');
stl_filename = strtrim(stl_filename);
if length(stl_filename)<4 || ~strcmpi(stl_filename(end-3:end),'.stl')
    stl_filename=[stl_filename,'.stl'];
end
while exist(stl_filename,'file')~=2
    warning('File "%s" not found.',stl_filename);
    stl_filename=input('  Re-enter: ','s');
    stl_filename=strtrim(stl_filename);
    if length(stl_filename)<4||~strcmpi(stl_filename(end-3:end),'.stl')
        stl_filename=[stl_filename,'.stl'];
    end
end

fprintf('\n--- FORCE DIRECTION ---\n  1->X  2->Y  3->Z\n');
wu=input('  Axis (1/2/3): ');
while ~ismember(wu,[1 2 3]); wu=input('  Enter 1,2 or 3: '); end
dir_labels={'X','Y','Z'};

fprintf('\n--- FIXED FACE ---\n  0->same as force  1->X-min  2->Y-min  3->Z-min\n');
fix_wu_in=input('  Choice (0/1/2/3): ');
while ~ismember(fix_wu_in,[0 1 2 3]); fix_wu_in=input('  Enter 0,1,2 or 3: '); end
fix_wu=fix_wu_in;
if fix_wu==0; fix_wu=wu; end

fprintf('\n--- FORCE ---\n  Negative=compression  Positive=tension\n');
Force=input('  Total force (N) [e.g. -800]: ');

fprintf('\n--- CROSS-SECTIONAL AREA ---\n');
fprintf('  Area perpendicular to loading direction (m^2).\n');
fprintf('  Example: 30mm x 30mm cube -> 0.03*0.03 = 9e-4 m^2\n');
AA=input('  A (m^2): ');
while AA<=0; AA=input('  Must be > 0 (m^2): '); end

fprintf('\n--- LINEAR MATERIAL  (|eps_est| < 10%%) ---\n');
fprintf('  Elastic 50A defaults: E=2.9e6 Pa, nu=0.40\n');
use_def=input('  Use defaults? [y/n]: ','s');
if strcmpi(strtrim(use_def),'n')
    EE=input('  E (Pa): ');
    nu=input('  nu (0<nu<0.5): ');
    while nu<=0||nu>=0.5; nu=input('  nu in (0,0.5): '); end
else
    EE=2.9e6; nu=0.40;
end

eps_est_full = abs(Force)/(AA*EE);
fprintf('\n  eps_est = |F|/(A*E0) = %.6f  (%.2f%%)\n',eps_est_full,eps_est_full*100);
if eps_est_full<0.10
    fprintf('  -> LINEAR regime expected\n\n');
else
    fprintf('  -> NONLINEAR regime expected\n\n');
end

fprintf('--- NONLINEAR MODEL  (|eps_est| >= 10%%) ---\n');
fprintf('  1 -> Neo-Hookean\n');
fprintf('  2 -> Mooney-Rivlin\n');
fprintf('  3 -> Yeoh\n');
fprintf('  4 -> EESM Arruda-Boyce\n');
fprintf('  5 -> Phenomenological\n');
hyper_choice=input('  Model (1/2/3/4/5): ');
while ~ismember(hyper_choice,[1 2 3 4 5]); hyper_choice=input('  Enter 1-5: '); end

% Initialise all params
mu_NH=0; C10=0; C01=0; C1=0; C2=0; C3=0;
mu_AB=0; N_AB=1; A1_AB=0; A2_AB=0; c_AB=0; f_AB=0;
alpha_phen=3;

switch hyper_choice
    case 1
        fprintf('\n  Neo-Hookean: W=mu/2*(I1-3). Default mu=0.57e6 Pa\n');
        q2=input('  Use default? [y/n]: ','s');
        if strcmpi(strtrim(q2),'n'); mu_NH=input('  mu (Pa): ');
        else; mu_NH=0.57e6; end
        C10=mu_NH/2; C1=mu_NH/2;
        hyper_label=sprintf('Neo-Hookean  mu=%.4g Pa',mu_NH);

    case 2
        fprintf('\n  Mooney-Rivlin. Defaults: C10=0.40e6, C01=0.085e6 Pa\n');
        q2=input('  Use defaults? [y/n]: ','s');
        if strcmpi(strtrim(q2),'n')
            C10=input('  C10 (Pa): '); C01=input('  C01 (Pa): ');
        else; C10=0.40e6; C01=0.085e6; end
        mu_NH=2*(C10+C01); C1=C10;
        hyper_label=sprintf('Mooney-Rivlin  C10=%.4g  C01=%.4g Pa',C10,C01);

    case 3
        fprintf('\n  Yeoh: W=C1*(I1-3)+C2*(I1-3)^2+C3*(I1-3)^3\n');
        fprintf('  Defaults: C1=1.6e5, C2=-1.2e4, C3=9e3 Pa  (E0=6*C1=0.96 MPa)\n');
        q2=input('  Use defaults? [y/n]: ','s');
        if strcmpi(strtrim(q2),'n')
            C1=input('  C1 (Pa): '); C2=input('  C2 (Pa): '); C3=input('  C3 (Pa): ');
        else; C1=1.6e5; C2=-1.2e4; C3=9e3; end
        mu_NH=2*C1; C10=C1;
        fprintf('  E0=6*C1=%.4g Pa  (input E0=%.4g Pa)\n',6*C1,EE);
        hyper_label=sprintf('Yeoh  C1=%.4g  C2=%.4g  C3=%.4g Pa',C1,C2,C3);

    case 4
        fprintf('\n  EESM Arruda-Boyce: WT=(1-f)*Wiso+f*Waniso\n');
        fprintf('  Paper Table II defaults (entire cube):\n');
        fprintf('    mu=9e5, N=9.5, A1=-1.89e6, A2=-8.2e5, c=-3.7e5, f=0.79\n');
        q2=input('  Use defaults? [y/n]: ','s');
        if strcmpi(strtrim(q2),'n')
            mu_AB=input('  mu (Pa): '); N_AB=input('  N: ');
            A1_AB=input('  A1 (Pa): '); A2_AB=input('  A2 (Pa): ');
            c_AB =input('  c (Pa): ');
            f_AB =input('  f (0<=f<=1): ');
            while f_AB<0||f_AB>1; f_AB=input('  f in [0,1]: '); end
        else
            mu_AB=9e5; N_AB=9.5; A1_AB=-1.89e6;
            A2_AB=-8.2e5; c_AB=-3.7e5; f_AB=0.79;
        end
        mu_NH=mu_AB;
        hyper_label=sprintf('EESM  mu=%.4g N=%.4g A1=%.4g A2=%.4g c=%.4g f=%.4g',...
                             mu_AB,N_AB,A1_AB,A2_AB,c_AB,f_AB);

    case 5
        fprintf('\n  Phenomenological: Etan=E0*(1+eps+eps/alpha)\n');
        fprintf('  Default alpha=3 (cube) -> Etan=E0*(1+4/3*eps)\n');
        q2=input('  Use default alpha=3? [y/n]: ','s');
        if strcmpi(strtrim(q2),'n')
            alpha_phen=input('  alpha (>0): ');
            while alpha_phen<=0; alpha_phen=input('  alpha>0: '); end
        else; alpha_phen=3; end
        mu_NH=EE/3;
        hyper_label=sprintf('Phenomenological  alpha=%.4g',alpha_phen);
end

fprintf('\n--- NEWTON-RAPHSON SETTINGS ---\n');
fprintf('  NOTE: min 5 iterations recommended for EESM / soft materials.\n');
n_steps=input('  Load increments [default 10]: ');
if isempty(n_steps)||n_steps<1; n_steps=10; end
max_iter=input('  Max NR iterations/step [default 20]: ');
if isempty(max_iter)||max_iter<1; max_iter=20; end
if max_iter<5
    fprintf('  *** WARNING: max_iter=%d is very low.\n',max_iter);
    fprintf('  *** With 1 iteration the NR loop cannot self-correct\n');
    fprintf('  *** bad initial strains. Recommended >= 5.\n');
    fprintf('  *** Results may be unreliable for nonlinear/EESM models.\n');
end
tol_res=input('  Rel. residual tolerance [default 1e-4]: ');
if isempty(tol_res)||tol_res<=0; tol_res=1e-4; end

order_fit=10;

% =========================================================================
%  SECTION 2 — SUMMARY
% =========================================================================

fprintf('\n=============================================================\n');
fprintf('  SIMULATION SUMMARY\n');
fprintf('=============================================================\n');
fprintf('  STL file         : %s  (%s)\n',stl_filename,unit_label);
fprintf('  Force            : %.4g N  (%s axis)\n',Force,dir_labels{wu});
fprintf('  Fixed face       : %s-min\n',dir_labels{fix_wu});
fprintf('  Cross-section A  : %.4g m^2\n',AA);
fprintf('  eps_est (full)   : %.4f  (%.2f%%)\n',eps_est_full,eps_est_full*100);
fprintf('  Linear E0        : %.4g Pa\n',EE);
fprintf('  nu               : %.4g\n',nu);
fprintf('  Nonlinear model  : %s\n',hyper_label);
fprintf('  Load steps       : %d\n',n_steps);
fprintf('  NR max iter/step : %d\n',max_iter);
fprintf('  NR tolerance     : %.2g\n',tol_res);
fprintf('=============================================================\n\n');
prc=input('  Proceed? [y/n]: ','s');
if ~strcmpi(strtrim(prc),'y'); fprintf('  Cancelled.\n'); return; end

% =========================================================================
%  SECTION 3 — GEOMETRY & MESH
% =========================================================================

fprintf('\n  Importing geometry ...\n');
smodel=createpde('structural','static-solid');
importGeometry(smodel,stl_filename);
figure('Name','Geometry','NumberTitle','off');
pdegplot(smodel); view(30,30); title(['Geometry: ',stl_filename]);

fprintf('  Generating mesh ...\n');
model=generateMesh(smodel,'GeometricOrder','linear');
nodeCoordinates=model.Nodes;
connectivity   =model.Elements;
numNodes=size(nodeCoordinates,2);
numElem =size(connectivity,2);
ndof    =3*numNodes;
fprintf('     Nodes=%d   Elements=%d\n\n',numNodes,numElem);

% Bounding-box size in metres (used for sanity check later)
bbox    = max(nodeCoordinates,[],2)-min(nodeCoordinates,[],2);
obj_size= max(bbox)*10^(-q);   % largest dimension in metres

% =========================================================================
%  SECTION 4 — COORDINATE EXTRACTION  (STL units -> metres)
% =========================================================================

xc=zeros(4,numElem); yc=zeros(4,numElem); zc=zeros(4,numElem);
for i=1:numElem
    for j=1:4
        n=connectivity(j,i);
        xc(j,i)=nodeCoordinates(1,n)*10^(-q);
        yc(j,i)=nodeCoordinates(2,n)*10^(-q);
        zc(j,i)=nodeCoordinates(3,n)*10^(-q);
    end
end

% =========================================================================
%  SECTION 5 — BOUNDARY CONDITIONS
% =========================================================================

load_coord   =nodeCoordinates(wu,:);
load_max     =round(max(load_coord),order_fit);
load_Indices =find(load_coord >= load_max - 1e-6*abs(load_max));

fixed_coord  =nodeCoordinates(fix_wu,:);
fixed_min    =round(min(fixed_coord),order_fit);
fixed_Indices=find(fixed_coord <= fixed_min + 1e-6*abs(fixed_min) + 1e-12);

fixedDOFs=reshape([3*fixed_Indices(:)'-2; ...
                   3*fixed_Indices(:)'-1; ...
                   3*fixed_Indices(:)'],1,[]);
freeDOFs=setdiff(1:ndof,fixedDOFs);

fprintf('  Load nodes  : %d  (%s-max)\n',numel(load_Indices),dir_labels{wu});
fprintf('  Fixed nodes : %d  (%s-min)\n\n',numel(fixed_Indices),dir_labels{fix_wu});

force_full=zeros(ndof,1);
for n=load_Indices(:)'
    force_full(3*(n-1)+wu)=Force/numel(load_Indices);
end

% =========================================================================
%  SECTION 6 — PRE-COMPUTE B-MATRICES & VOLUMES
% =========================================================================

fprintf('  Pre-computing B-matrices ...\n');
B_all=cell(numElem,1); V_all=zeros(numElem,1);

for jj=1:numElem
    x1=xc(1,jj); x2=xc(2,jj); x3=xc(3,jj); x4=xc(4,jj);
    y1=yc(1,jj); y2=yc(2,jj); y3=yc(3,jj); y4=yc(4,jj);
    z1=zc(1,jj); z2=zc(2,jj); z3=zc(3,jj); z4=zc(4,jj);

    v6=det([1 x1 y1 z1; 1 x2 y2 z2; 1 x3 y3 z3; 1 x4 y4 z4]);
    Ve=abs(v6)/6; V_all(jj)=Ve; sc=1/(6*Ve);

    b1=-((y3-y2)*(z4-z2)-(y4-y2)*(z3-z2));
    g1=-((x4-x2)*(z3-z2)-(x3-x2)*(z4-z2));
    d1=-((x3-x2)*(y4-y2)-(x4-x2)*(y3-y2));
    b2= (y3-y1)*(z4-z1)-(y4-y1)*(z3-z1);
    g2= (x4-x1)*(z3-z1)-(x3-x1)*(z4-z1);
    d2= (x3-x1)*(y4-y1)-(x4-x1)*(y3-y1);
    b3=-((y2-y1)*(z4-z1)-(y4-y1)*(z2-z1));
    g3=-((x4-x1)*(z2-z1)-(x2-x1)*(z4-z1));
    d3=-((x2-x1)*(y4-y1)-(x4-x1)*(y2-y1));
    b4= (y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);
    g4= (x3-x1)*(z2-z1)-(x2-x1)*(z3-z1);
    d4= (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);

    Bi=@(b,g,d) sc*[b 0 0; 0 g 0; 0 0 d; g b 0; 0 d g; d 0 b];
    B_all{jj}=[Bi(b1,g1,d1) Bi(b2,g2,d2) Bi(b3,g3,d3) Bi(b4,g4,d4)];
end
fprintf('     Done.\n\n');

% =========================================================================
%  SECTION 7 — LINEAR D  (|eps| < 10%)
% =========================================================================

c_lin=EE/((1+nu)*(1-2*nu));
D_lin=c_lin*[1-nu nu    nu    0          0          0;
             nu   1-nu  nu    0          0          0;
             nu   nu    1-nu  0          0          0;
             0    0     0     (1-2*nu)/2 0          0;
             0    0     0     0          (1-2*nu)/2 0;
             0    0     0     0          0          (1-2*nu)/2];

% =========================================================================
%  SECTION 8 — EESM DIAGNOSTIC  (model 4 only, printed before solve)
% =========================================================================

if hyper_choice==4
    fprintf('--- EESM DIAGNOSTIC at full-load operating point ---\n');
    ed = eps_est_full;

    % Safe lam_min: solve lam^2+2/lam = 3*N*0.95^2 by Newton iteration
    tI1=3*N_AB*0.95^2; ls=0.3;
    for kk=1:50
        fl=ls^2+2/ls-tI1; dfl=2*ls-2/ls^2;
        ls=max(0.001,min(0.999,ls-fl/dfl));
        if abs(fl)<1e-12; break; end
    end
    lam_min_diag=max(0.001,ls);
    eps_max_diag=1-lam_min_diag;

    ld =max(lam_min_diag,min(0.999,1-ed));
    l2d=1/sqrt(ld);
    I1d=ld^2+l2d^2+l2d^2;
    I3d=(ld*l2d*l2d)^2;
    lrd=min(sqrt(I1d/(3*N_AB)),0.9999);
    bd =(3*lrd)/(1-lrd^3);

    if isfinite(bd) && bd<=700
        wid=mu_AB*(N_AB*(bd*lrd+log(bd/sinh(bd)))-log(bd/lrd))+c_AB;
    else
        wid=NaN;
    end
    I1m3d=I1d-3;
    wad=(A1_AB/3)*I1m3d+(A2_AB/9)*I1m3d^2-(A1_AB/3)*log(max(I3d,1e-12));

    if isfinite(wid)
        WTd=(1-f_AB)*wid+f_AB*wad;
    else
        WTd=NaN;
    end

    sig_d=eesm_sigma(ld,mu_AB,N_AB,A1_AB,A2_AB,c_AB,f_AB);
    Et_d =tangent_Et(ed,4,EE,nu,mu_NH,C10,C01,C1,C2,C3,...
                     mu_AB,N_AB,A1_AB,A2_AB,c_AB,f_AB,alpha_phen);
    if isfinite(WTd)&&ed>0; Feq28=2*AA*WTd/ed; else; Feq28=NaN; end

    fprintf('  eps_est            = %.6f\n',ed);
    fprintf('  eps_max_safe       = %.4f (%.1f%%) [EESM valid range]\n',...
            eps_max_diag,eps_max_diag*100);
    if ed>=eps_max_diag
        fprintf('  *** eps_est > eps_max_safe: PHENOMENOLOGICAL FALLBACK used ***\n');
    end
    fprintf('  lam1 (= lam)       = %.6f\n',ld);
    fprintf('  lam2 = lam3        = %.6f\n',l2d);
    fprintf('  I1                 = %.6f   (check lam^2+2/lam = %.6f)\n',I1d,ld^2+2/ld);
    fprintf('  I3                 = %.8f   (expected 1.0)\n',I3d);
    fprintf('  lamr               = %.6f\n',lrd);
    fprintf('  beta (1-lamr^3)    = %.6f\n',bd);
    fprintf('  Wiso               = %.4f Pa\n',wid);
    fprintf('  Waniso             = %.4f Pa\n',wad);
    fprintf('  WT = (1-f)*Wiso + f*Waniso\n');
    fprintf('     = %.4f*%.4f + %.4f*%.4f = %.4f Pa\n',...
            1-f_AB,wid,f_AB,wad,WTd);

    % ── Cancellation warning ──────────────────────────────────────
    if isfinite(wid) && abs(wid)>0
        cancel_ratio=abs(WTd)/abs(wid);
        if cancel_ratio<0.05
            fprintf('  *** CANCELLATION WARNING: |WT|/|Wiso| = %.3f%%\n',...
                    cancel_ratio*100);
            fprintf('  *** Wiso and Waniso nearly cancel -> Et_sec << E0\n');
            fprintf('  *** (1-f)*Wiso = %.2f Pa,  f*Waniso = %.2f Pa\n',...
                    (1-f_AB)*wid, f_AB*wad);
            fprintf('  *** Likely cause: A1/A2 too small or f too large\n');
            fprintf('  *** for these mu/N values. Check parameter fitting.\n');
        end
    end

    fprintf('  sigma = lam*dW/dlam= %.4f Pa  (%.4f kPa)\n',sig_d,sig_d/1e3);
    fprintf('  Et used in FEM     = %.4f Pa  (%.4f MPa)\n',Et_d,Et_d/1e6);
    fprintf('  F via Eq.(28)      = %.4f N   (applied = %.4g N)\n',Feq28,Force);

    % Et vs E0 ratio warning
    if Et_d < EE*0.01
        fprintf('  *** Et << E0 WARNING: Et = %.4g Pa,  E0 = %.4g Pa\n',Et_d,EE);
        fprintf('  *** Stiffness matrix will be nearly singular.\n');
        fprintf('  *** Results will be unreliable with these parameters.\n');
    end

    % 1D displacement estimate
    if Et_d>1e4
        fprintf('  1D delta estimate  = %.4f mm\n\n',...
                abs(Force)*obj_size/(Et_d*AA)*1e3);
    else
        fprintf('  1D delta estimate  = unreliable (Et too small)\n\n');
    end
end

% =========================================================================
%  SECTION 9 — NEWTON-RAPHSON SOLVE
% =========================================================================

fprintf('  Newton-Raphson: %d steps, max %d iter, tol=%.2g\n\n',...
        n_steps,max_iter,tol_res);

u_total     =zeros(ndof,1);
F_step_inc  =force_full/n_steps;
conv_history=zeros(n_steps,1);
iter_history=zeros(n_steps,1);
eps_elem    =zeros(numElem,1);

for step=1:n_steps

    F_applied   =F_step_inc*step;
    eps_est_step=eps_est_full*(step/n_steps);

    fprintf('  -- Step %2d/%d  |F|=%.4g N  eps_est=%.4f (%.1f%%) --\n',...
            step,n_steps,abs(Force)*step/n_steps,...
            eps_est_step,eps_est_step*100);

    converged=false;

    for iter=1:max_iter

        % (A) Element axial strain
        if iter==1
            eps_elem(:)=eps_est_step;
            seed_str='seeded from eps_est';
        else
            for jj=1:numElem
                ns   =connectivity(:,jj);
                dof_e=reshape([3*ns-2,3*ns-1,3*ns]',1,[]);
                eps_v=B_all{jj}*u_total(dof_e);
                eps_elem(jj)=max(abs(eps_v(wu)),eps_est_step);
            end
            seed_str='from |B*u|';
        end

        % (B) Build K_T and F_int element-by-element
        global_KT  =sparse(ndof,ndof);
        F_int      =zeros(ndof,1);
        elem_regime=zeros(numElem,1);

        for jj=1:numElem
            ns   =connectivity(:,jj);
            dof_e=reshape([3*ns-2,3*ns-1,3*ns]',1,[]);
            u_e  =u_total(dof_e);
            Ve   =V_all(jj);
            BB_e =B_all{jj};

            De=get_De(eps_elem(jj),D_lin,nu,hyper_choice,EE,...
                      mu_NH,C10,C01,C1,C2,C3,...
                      mu_AB,N_AB,A1_AB,A2_AB,c_AB,f_AB,alpha_phen);
            if abs(eps_elem(jj))>=0.10; elem_regime(jj)=1; end

            global_KT(dof_e,dof_e)=global_KT(dof_e,dof_e)+Ve*(BB_e'*De*BB_e);
            sigma_e     =De*BB_e*u_e;
            F_int(dof_e)=F_int(dof_e)+Ve*(BB_e'*sigma_e);
        end

        % (C) Residual and solve
        R           =F_applied-F_int;
        R(fixedDOFs)=0;

        % (D) Solve and update displacement
        du_free           =global_KT(freeDOFs,freeDOFs)\R(freeDOFs);
        u_total(freeDOFs) =u_total(freeDOFs)+du_free;

        % (E) Re-evaluate residual AFTER update
        F_int2=zeros(ndof,1);
        for jj=1:numElem
            ns   =connectivity(:,jj);
            dof_e=reshape([3*ns-2,3*ns-1,3*ns]',1,[]);
            De   =get_De(eps_elem(jj),D_lin,nu,hyper_choice,EE,...
                         mu_NH,C10,C01,C1,C2,C3,...
                         mu_AB,N_AB,A1_AB,A2_AB,c_AB,f_AB,alpha_phen);
            sig2        =De*B_all{jj}*u_total(dof_e);
            F_int2(dof_e)=F_int2(dof_e)+V_all(jj)*(B_all{jj}'*sig2);
        end
        R2           =F_applied-F_int2;
        R2(fixedDOFs)=0;
        rel_res=norm(R2(freeDOFs))/max(norm(F_applied(freeDOFs)),1e-12);

        n_nl=sum(elem_regime);
        fprintf('    iter %2d [%s]  rel_res=%.4e | NL elems: %d/%d (%.0f%%)\n',...
                iter,seed_str,rel_res,n_nl,numElem,100*n_nl/numElem);

        % (F) Physical displacement sanity check
        max_u_now=max(abs(u_total));
        if max_u_now > 5*obj_size
            fprintf('    *** SANITY WARNING: max disp = %.4g m > 5x object size (%.4g m)\n',...
                    max_u_now,obj_size);
            if hyper_choice==4
                fprintf('    *** EESM near-cancellation likely. Check A1/A2/f values.\n');
                fprintf('    *** Et_sec was %.4g Pa (E0=%.4g Pa). Need more iterations.\n',...
                        tangent_Et(eps_est_step,4,EE,nu,mu_NH,C10,C01,C1,C2,C3,...
                                   mu_AB,N_AB,A1_AB,A2_AB,c_AB,f_AB,alpha_phen), EE);
            end
        end

        if rel_res<tol_res
            fprintf('    -> Converged (%d iter).\n\n',iter);
            conv_history(step)=rel_res;
            iter_history(step)=iter;
            converged=true;
            break;
        end
    end

    if ~converged
        warning('Step %d: not converged after %d iter (rel_res=%.3e)',...
                step,max_iter,rel_res);
        conv_history(step)=rel_res;
        iter_history(step)=max_iter;
    end
end

fprintf('  Solve complete.\n\n');

% =========================================================================
%  SECTION 10 — REACTION FORCES
% =========================================================================

F_int_final=zeros(ndof,1);
for jj=1:numElem
    ns   =connectivity(:,jj);
    dof_e=reshape([3*ns-2,3*ns-1,3*ns]',1,[]);
    De   =get_De(eps_elem(jj),D_lin,nu,hyper_choice,EE,...
                 mu_NH,C10,C01,C1,C2,C3,...
                 mu_AB,N_AB,A1_AB,A2_AB,c_AB,f_AB,alpha_phen);
    sigma_e=De*B_all{jj}*u_total(dof_e);
    F_int_final(dof_e)=F_int_final(dof_e)+V_all(jj)*(B_all{jj}'*sigma_e);
end
rxn_dofs      =3*(fixed_Indices(:)')+fix_wu-3;
reaction_vec  =F_int_final-force_full;
total_reaction=sum(reaction_vec(rxn_dofs));

% =========================================================================
%  SECTION 11 — POST-PROCESSING
% =========================================================================

disp_in_dir=u_total(wu:3:end);
z_disp     =u_total(3:3:end);
max_disp   =max(abs(disp_in_dir));
maxz       =max(abs(z_disp));

for jj=1:numElem
    ns   =connectivity(:,jj);
    dof_e=reshape([3*ns-2,3*ns-1,3*ns]',1,[]);
    eps_v=B_all{jj}*u_total(dof_e);
    eps_elem(jj)=max(abs(eps_v(wu)),eps_est_full);
end

vm_stress      =zeros(numElem,1);
elem_regime_fin=zeros(numElem,1);
for jj=1:numElem
    ns   =connectivity(:,jj);
    dof_e=reshape([3*ns-2,3*ns-1,3*ns]',1,[]);
    De   =get_De(eps_elem(jj),D_lin,nu,hyper_choice,EE,...
                 mu_NH,C10,C01,C1,C2,C3,...
                 mu_AB,N_AB,A1_AB,A2_AB,c_AB,f_AB,alpha_phen);
    if abs(eps_elem(jj))>=0.10; elem_regime_fin(jj)=1; end
    sig =De*B_all{jj}*u_total(dof_e);
    sxx=sig(1); syy=sig(2); szz=sig(3);
    txy=sig(4); tyz=sig(5); txz=sig(6);
    vm_stress(jj)=sqrt(0.5*((sxx-syy)^2+(syy-szz)^2+(szz-sxx)^2 ...
                             +6*(txy^2+tyz^2+txz^2)));
end

node_vm    =zeros(numNodes,1); nc1=zeros(numNodes,1);
node_regime=zeros(numNodes,1); nc2=zeros(numNodes,1);
for jj=1:numElem
    ns=connectivity(:,jj);
    node_vm(ns)    =node_vm(ns)    +vm_stress(jj);       nc1(ns)=nc1(ns)+1;
    node_regime(ns)=node_regime(ns)+elem_regime_fin(jj); nc2(ns)=nc2(ns)+1;
end
node_vm    =node_vm    ./max(nc1,1);
node_regime=node_regime./max(nc2,1);

% =========================================================================
%  SECTION 12 — VISUALIZATIONS
% =========================================================================

figure('Name',sprintf('Displacement %s',dir_labels{wu}),'NumberTitle','off');
pdeplot3D(model,'ColorMapData',disp_in_dir*1e3); colormap jet;
cb=colorbar; cb.Label.String=sprintf('u_%s (mm)',lower(dir_labels{wu}));
title(sprintf('Displacement — %s',dir_labels{wu})); view(30,30);

figure('Name','Z-Displacement','NumberTitle','off');
pdeplot3D(model,'ColorMapData',z_disp*1e3); colormap jet;
cb2=colorbar; cb2.Label.String='u_z (mm)'; title('Z-displacement'); view(30,30);

figure('Name','Von Mises Stress','NumberTitle','off');
pdeplot3D(model,'ColorMapData',node_vm*1e-3); colormap jet;
cb3=colorbar; cb3.Label.String='Von Mises (kPa)'; title('Von Mises Stress'); view(30,30);

figure('Name','Regime Map','NumberTitle','off');
pdeplot3D(model,'ColorMapData',node_regime);
colormap([0.18 0.36 0.80; 0.85 0.20 0.20]);
cb4=colorbar; cb4.Ticks=[0 1];
cb4.TickLabels={'Linear (|eps|<10%)','Nonlinear (|eps|>=10%)'};
title('Regime: Blue=Linear  Red=Hyperelastic'); view(30,30);

figure('Name','NR Convergence','NumberTitle','off');
semilogy(1:n_steps,conv_history,'o-b','LineWidth',1.5,'MarkerFaceColor','b'); hold on;
yline(tol_res,'--r','LineWidth',1.2,'Label',sprintf('tol=%.1e',tol_res));
xlabel('Load step'); ylabel('Rel. residual ||R||/||F||');
title('Newton-Raphson convergence per load step'); grid on;

figure('Name','Z-disp cursor','NumberTitle','off');
pdeplot3D(model,'ColorMapData',z_disp*1e3); colormap jet;
colorbar; title('Z-disp (data cursor on)'); view(30,30);
dcm=datacursormode(gcf); set(dcm,'Enable','on');
set(dcm,'UpdateFcn',...
    @(src,evt) cursor_cb(src,evt,nodeCoordinates,z_disp,unit_label));

% =========================================================================
%  SECTION 13 — SUMMARY TABLE
% =========================================================================

fprintf('=============================================================\n');
fprintf('  NONLINEAR SOLVE SUMMARY (v10)\n');
fprintf('=============================================================\n');
fprintf('  eps_est (full load) : %.4f  (%.2f%%)\n',eps_est_full,eps_est_full*100);
if hyper_choice==4
    fprintf('  EESM lam (1-eps)    : %.4f\n',max(0.001,1-eps_est_full));
    fprintf('  EESM modulus        : SECANT Et=sigma/eps\n');
    fprintf('                        (phenomenological fallback if eps>eps_max)\n');
elseif hyper_choice==3
    fprintf('  lam = 1+eps_est     : %.4f\n',1+eps_est_full);
    fprintf('  Yeoh E0 check       : 6*C1=%.4g Pa (input E0=%.4g Pa)\n',6*C1,EE);
else
    fprintf('  lam = 1+eps_est     : %.4f\n',1+eps_est_full);
end
fprintf('  Nonlinear model     : %s\n',hyper_label);
fprintf('  %-6s  %-6s  %-14s  %s\n','Step','Iters','Rel Residual','Status');
fprintf('  %s\n',repmat('-',1,42));
for s=1:n_steps
    tag='OK'; if conv_history(s)>tol_res; tag='WARN'; end
    fprintf('  %-6d  %-6d  %-14.4e  %s\n',...
            s,iter_history(s),conv_history(s),tag);
end
fprintf('=============================================================\n');
fprintf('  Applied force        : %.4g N\n',Force);
fprintf('  Reaction force       : %.4g N\n',total_reaction);
fprintf('  Max disp (%s)       : %.4g m  (%.4g mm)\n',...
        dir_labels{wu},max_disp,max_disp*1e3);
fprintf('  Max |z_disp|         : %.4g m  (%.4g mm)\n',maxz,maxz*1e3);
fprintf('  Max Von Mises        : %.4g Pa  (%.4g kPa)\n',...
        max(node_vm),max(node_vm)*1e-3);
fprintf('  Nonlinear elements   : %d / %d  (%.1f%%)\n',...
        sum(elem_regime_fin),numElem,100*sum(elem_regime_fin)/numElem);
fprintf('=============================================================\n');
fprintf('\n  maxz = %.6g m  (%.4g mm)\n\n',maxz,maxz*1e3);

% Final sanity check on result
if max_disp > 5*obj_size
    fprintf('  *** RESULT WARNING: displacement (%.4g m) > 5x object size.\n',...
            max_disp);
    fprintf('  *** This is likely unphysical. Possible causes:\n');
    fprintf('  *** 1. EESM parameter cancellation (check diagnostic above)\n');
    fprintf('  *** 2. Too few NR iterations (used %d, try >= 10)\n',max_iter);
    fprintf('  *** 3. Load too large for material parameters\n');
end

end  % closes FEM_Solver_Nonlinear_v10()

% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

% ------------------------------------------------------------------
%  W_eesm — EESM strain energy density  (paper Eqs.21-27)
%  beta = 3*lamr/(1-lamr^3)  per paper Eq.(24)
%  Returns NaN if beta>700 to prevent sinh() overflow
% ------------------------------------------------------------------
function w=W_eesm(lam,mu_AB,N_AB,A1_AB,A2_AB,c_AB,f_AB)
    lam1=max(lam,0.001);
    lam2=1/sqrt(lam1);
    I1  =lam1^2+lam2^2+lam2^2;
    I3  =(lam1*lam2*lam2)^2;
    lamr=min(sqrt(I1/(3*N_AB)),0.9999);
    beta=(3*lamr)/(1-lamr^3);
    if ~isfinite(beta)||beta>700; w=NaN; return; end
    sb=sinh(beta);
    if sb==0||~isfinite(sb); w=NaN; return; end
    w_iso=mu_AB*(N_AB*(beta*lamr+log(beta/sb))-log(beta/lamr))+c_AB;
    if ~isfinite(w_iso); w=NaN; return; end
    w_aniso=(A1_AB/3)*(I1-3)+(A2_AB/9)*(I1-3)^2-(2.*A1_AB/3)*log((I3));
    w=(1-f_AB)*w_iso+f_AB*w_aniso;



end

% ------------------------------------------------------------------
%  eesm_sigma — uniaxial Cauchy stress = lam * dW/dlam
%  Positive for compression (lam<1): material resists loading.
% ------------------------------------------------------------------
function s=eesm_sigma(lam,mu_AB,N_AB,A1_AB,A2_AB,c_AB,f_AB)
    dh=1e-7;
    wp=W_eesm(lam+dh,mu_AB,N_AB,A1_AB,A2_AB,c_AB,f_AB);
    wm=W_eesm(lam-dh,mu_AB,N_AB,A1_AB,A2_AB,c_AB,f_AB);
    if ~isfinite(wp)||~isfinite(wm); s=NaN; return; end
    s=lam*(wp-wm)/(2*dh);
end

% ------------------------------------------------------------------
%  yeoh_sig — uniaxial Cauchy stress from Yeoh  (paper Eq.48)
% ------------------------------------------------------------------
function s=yeoh_sig(lam,C1,C2,C3)
    I1m3=lam^2+2/lam-3;
    dW  =C1+2*C2*I1m3+3*C3*I1m3^2;
    s   =2*dW*(lam-1/lam^2);
end

% ------------------------------------------------------------------
%  tangent_Et — effective Young's modulus for element D matrix
%
%  Models 1,2,3,5 : tangent dσ/dλ at lam=1+|eps|  (tension side)
%  Model 4 EESM   : secant  σ/eps at lam=1-|eps|  (compression)
%    Fallback to phenomenological Etan=E0*(1+eps+eps/3) when:
%      - eps > eps_max_safe (~92% for N=10), OR
%      - W_eesm returns NaN/Inf, OR
%      - sigma <= 0 (wrong sign)
%    Floor: Et >= 1 kPa always.
% ------------------------------------------------------------------
function Et=tangent_Et(eps_e,hc,EE,nu,...
                       mu_NH,C10,C01,C1,C2,C3,...
                       mu_AB,N_AB,A1_AB,A2_AB,c_AB,f_AB,...
                       alpha_phen)
    dh=1e-7;
    switch hc
        case 1
            lam=1+abs(eps_e);
            Et =mu_NH*(1+2/lam^3);

        case 2
            lam=1+abs(eps_e);
            sig=@(l) 2*(l-1/l^2)*(C10+C01/l);
            Et =(sig(lam+dh)-sig(lam-dh))/(2*dh);

        case 3
            lam=1+abs(eps_e);
            sig=@(l) yeoh_sig(l,C1,C2,C3);
            Et =(sig(lam+dh)-sig(lam-dh))/(2*dh);

        case 4
            eps_abs=abs(eps_e);
            if eps_abs<1e-6
                Et=EE;
            else
                % Safe lam_min: solve lam^2+2/lam = 3*N*0.95^2
                tI1=3*N_AB*0.95^2; ls=0.3;
                for kk=1:50
                    fl=ls^2+2/ls-tI1; dfl=2*ls-2/ls^2;
                    ls=max(0.001,min(0.999,ls-fl/dfl));
                    if abs(fl)<1e-12; break; end
                end
                lam_min_safe=max(0.001,ls);
                eps_max_safe=1-lam_min_safe;

                if eps_abs>=eps_max_safe
                    % Beyond physical range: phenomenological fallback
                    Et=EE*(1+eps_abs+eps_abs/3);
                else
                    lam=1-eps_abs;
                    sig_v=eesm_sigma(lam,mu_AB,N_AB,A1_AB,A2_AB,c_AB,f_AB);
                    if ~isfinite(sig_v)||sig_v<=0
                        Et=EE*(1+eps_abs+eps_abs/3);
                    else
                        Et=sig_v/eps_abs;   % secant modulus, always > 0
                    end
                end
            end

        case 5
            Et=EE*(1+abs(eps_e)+abs(eps_e)/alpha_phen);
    end
    Et=max(Et,1e3);   % floor: 1 kPa minimum
end

% ------------------------------------------------------------------
%  build_Dt — isotropic tangent constitutive matrix  (paper Eq.32)
% ------------------------------------------------------------------
function Dt=build_Dt(Et,nu_t)
    ct=Et/((1+nu_t)*(1-2*nu_t));
    Dt=ct*[1-nu_t nu_t   nu_t   0            0            0;
           nu_t   1-nu_t nu_t   0            0            0;
           nu_t   nu_t   1-nu_t 0            0            0;
           0      0      0      (1-2*nu_t)/2 0            0;
           0      0      0      0            (1-2*nu_t)/2 0;
           0      0      0      0            0            (1-2*nu_t)/2];
end

% ------------------------------------------------------------------
%  get_De — pick linear or nonlinear D for one element
% ------------------------------------------------------------------
function De=get_De(eps_ax,D_lin_in,nu_in,hc,EE_in,...
                   mu_NH_in,C10_in,C01_in,C1_in,C2_in,C3_in,...
                   mu_AB_in,N_AB_in,A1_AB_in,A2_AB_in,c_AB_in,f_AB_in,...
                   alpha_in)
    if abs(eps_ax)<0.10
        De=D_lin_in;
    else
        Et=tangent_Et(eps_ax,hc,EE_in,nu_in,...
                      mu_NH_in,C10_in,C01_in,C1_in,C2_in,C3_in,...
                      mu_AB_in,N_AB_in,A1_AB_in,A2_AB_in,c_AB_in,f_AB_in,...
                      alpha_in);
        De=build_Dt(Et,nu_in);
    end
end

% ------------------------------------------------------------------
%  cursor_cb — data cursor callback for Z-disp figure
% ------------------------------------------------------------------
function txt=cursor_cb(~,event_obj,nodeCoords,z_displacement,unit_label)
    pos=event_obj.Position;
    [~,id]=min(sqrt(sum((nodeCoords'-pos).^2,2)));
    txt={sprintf('Node : %d',id),...
         sprintf('X    : %.4f %s',nodeCoords(1,id),unit_label),...
         sprintf('Y    : %.4f %s',nodeCoords(2,id),unit_label),...
         sprintf('Z    : %.4f %s',nodeCoords(3,id),unit_label),...
         sprintf('u_z  : %.4e m',z_displacement(id))};
end

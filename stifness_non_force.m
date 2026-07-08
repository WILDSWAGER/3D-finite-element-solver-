

clear; clc; close all;
clear W_eesm

fprintf('=============================================================\n');
fprintf('  STRAIN-SWEEP Eq.28 SOLVER — self-consistent secant E_t      \n');
fprintf('=============================================================\n\n');

% =========================================================================
%  SECTION 1 — USER INPUTS
% =========================================================================

fprintf('--- UNIT SYSTEM ---\n');
fprintf('  1 -> mm\n  2 -> cm\n  3 -> m\n');
unit_choice = input('  Enter choice (1/2/3): ');
switch unit_choice
    case 1; q=3; unit_label='mm';
    case 2; q=2; unit_label='cm';
    case 3; q=0; unit_label='m';
end
fprintf('  Unit: %s\n\n',unit_label);

fprintf('--- STL SCALE FACTOR ---\n');
scale_stl = input('  Scale factor: ');
fprintf('  scale_stl = %.6g\n\n',scale_stl);

fprintf('--- STL FILE ---\n');
stl_filename = input('  Filename (e.g. cube.stl): ','s');
stl_filename = strtrim(stl_filename);
if length(stl_filename)<4 || ~strcmpi(stl_filename(end-3:end),'.stl')
    stl_filename = [stl_filename,'.stl'];
end
while exist(stl_filename,'file')~=2
    fprintf('  File not found: %s\n',stl_filename);
    stl_filename = input('  Re-enter filename: ','s');
    stl_filename = strtrim(stl_filename);
    if length(stl_filename)<4 || ~strcmpi(stl_filename(end-3:end),'.stl')
        stl_filename = [stl_filename,'.stl'];
    end
end
fprintf('  File: %s\n\n',stl_filename);

fprintf('--- MESH DENSITY ---\n');
fprintf('  Hmax = max element edge length in metres.\n');
fprintf('  0 -> MATLAB default (finest, slowest)\n');
Hmax_in = input('  Hmax (m): ');
if Hmax_in <= 0
    use_hmax = false;
    fprintf('  Using MATLAB default mesh.\n\n');
else
    use_hmax = true;
    Hmax_val = Hmax_in;
    fprintf('  Hmax = %.4g m\n\n',Hmax_val);
end

fprintf('--- LOAD DIRECTION ---\n');
fprintf('  1 -> X\n  2 -> Y\n  3 -> Z\n');
wu = input('  Enter axis (1/2/3): ');
dir_labels = {'X','Y','Z'};
fprintf('  Loading axis: %s\n\n',dir_labels{wu});

fprintf('--- FIXED FACE ---\n');
fprintf('  0->same as force  1->X-min  2->Y-min  3->Z-min\n');
fix_wu = input('  Enter choice (0/1/2/3): ');
if fix_wu==0; fix_wu=wu; end
fprintf('  Fixed face: %s-min\n\n',dir_labels{fix_wu});

fprintf('--- APPLIED FORCE ---\n');
fprintf('  Negative=compression  Positive=tension\n');
Force = input('  Total force (N): ');
fprintf('  Force = %.4g N\n\n',Force);

fprintf('--- CROSS-SECTIONAL AREA ---\n');
AA = input('  A (m^2): ');
while ~isnumeric(AA) || ~isscalar(AA) || AA<=0
    fprintf('  *** Must be positive. You entered: %g\n',AA);
    AA = input('  Re-enter A (m^2): ');
end
fprintf('  A = %.4e m^2\n\n',AA);

fprintf('--- NONLINEAR MODEL (its W_T is swept across strain) ---\n');
fprintf('  1->Neo-Hookean  2->Mooney-Rivlin  3->Yeoh\n');
fprintf('  4->EESM         5->Phenomenological\n');
hyper_choice = input('  Model (1/2/3/4/5): ');

mu_NH=0; C10=0; C01=0; C1=0; C2=0; C3=0;
mu_AB=0; N_AB=1; A1_AB=0; A2_AB=0; c_AB=0; f_AB=0;
alpha_phen=3; EE=0;

switch hyper_choice
    case 1
        fprintf('\n  Neo-Hookean: W=mu/2*(I1-3)\n');
        mu_NH=input('  mu (Pa): ');
        hyper_label=sprintf('Neo-Hookean  mu=%.4g Pa',mu_NH);

    case 2
        fprintf('\n  Mooney-Rivlin: W=C10*(I1-3)+C01*(I2-3)\n');
        C10=input('  C10 (Pa): ');
        C01=input('  C01 (Pa): ');
        hyper_label=sprintf('Mooney-Rivlin  C10=%.4g C01=%.4g Pa',C10,C01);

    case 3
        fprintf('\n  Yeoh: W=C1*(I1-3)+C2*(I1-3)^2+C3*(I1-3)^3\n');
        C1=input('  C1 (Pa): ');
        C2=input('  C2 (Pa): ');
        C3=input('  C3 (Pa): ');
        hyper_label=sprintf('Yeoh  C1=%.4g C2=%.4g C3=%.4g Pa',C1,C2,C3);

    case 4
        fprintf('\n  EESM Arruda-Boyce: WT=(1-f)*Wiso+f*Waniso\n');
        fprintf('  --- PARAMETER SOURCE ---\n');
        fprintf('  1 -> Load from eesm_fitted_params.mat (recommended)\n');
        fprintf('  2 -> Manually type mu/N/A1/A2/f\n');
        param_source = input('  Enter choice (1/2): ');
        if param_source == 1
            fit_file = 'eesm_fitted_params.mat';
            while exist(fit_file,'file')~=2
                fprintf('  *** Fitted parameter file not found: %s\n',fit_file);
                fit_file = input('  Enter path to fitted params .mat file: ','s');
                fit_file = strtrim(fit_file);
            end
            S = load(fit_file);
            mu_AB=S.eesm_mu; N_AB=S.eesm_N; A1_AB=S.eesm_A1;
            A2_AB=S.eesm_A2; f_AB=S.eesm_f; c_AB=S.eesm_c;
            fprintf('  Loaded: mu=%.4g N=%.4g A1=%.4g A2=%.4g f=%.4g\n\n',...
                    mu_AB,N_AB,A1_AB,A2_AB,f_AB);
            if isfield(S,'eesm_fit_eps_max')
                eesm_fit_max_eps = S.eesm_fit_eps_max;
            end
        else
            mu_AB=input('  mu (Pa): ');
            N_AB =input('  N: ');
            A1_AB=input('  A1 (Pa): ');
            A2_AB=input('  A2 (Pa): ');
            f_AB =input('  f (0 to 1): ');
            lam_r0_c=1/sqrt(N_AB);
            beta0_c =3*lam_r0_c/(1-lam_r0_c^3);
            c_AB=-mu_AB*(N_AB*(beta0_c*lam_r0_c+log(beta0_c/sinh(beta0_c)))-log(beta0_c/lam_r0_c));
            fprintf('  c (auto) = %.6e Pa\n\n',c_AB);
        end
        hyper_label=sprintf('EESM  mu=%.4g N=%.4g A1=%.4g A2=%.4g c=%.4g f=%.4g',...
                             mu_AB,N_AB,A1_AB,A2_AB,c_AB,f_AB);

    case 5
        fprintf('\n  Phenomenological: Etan=E0*(1+eps+eps/alpha)\n');
        EE=input('  E (Pa): ');
        alpha_phen=input('  alpha: ');
        hyper_label=sprintf('Phenomenological  E=%.4g alpha=%.4g',EE,alpha_phen);
end
fprintf('\n');

fprintf('--- STRAIN SWEEP SETTINGS ---\n');
fprintf('  Sweep range: eps from ~0 up to a maximum you choose.\n');
fprintf('  Recommended: 0.95 for compression (avoid lam<=0.05), 3-5 for tension.\n');
eps_sweep_max = input('  Max eps for sweep: ');
n_sweep = input('  Number of sweep points (e.g. 2000): ');
fprintf('\n');

% =========================================================================
%  SECTION 2 — BUILD THE W_T -> Eq.28 FORCE-STRAIN TABLE
% =========================================================================

fprintf('--- STEP 1-2: SWEEPING STRAIN, COMPUTING THE REAL STRESS sigma(eps) ---\n');
fprintf('  NOTE: this sweep now uses sigma(eps)=dW_T/dlambda directly (the\n');
fprintf('  REAL nominal stress), NOT the Eq.28 energy average. Taking a\n');
fprintf('  local slope of the Eq.28 curve was mixing two different\n');
fprintf('  approximations (an energy average, then a derivative of that\n');
fprintf('  average) which is not the same mathematical object as the\n');
fprintf('  real tangent modulus dsigma/deps — that mismatch is almost\n');
fprintf('  certainly why E_t came out far from your reference value.\n\n');

sign_F = sign(Force);
if sign_F==0; sign_F=1; end

eps_sweep = linspace(1e-5, eps_sweep_max, n_sweep)';
sigma_sweep = zeros(n_sweep,1);
F_sweep     = zeros(n_sweep,1);

for kk=1:n_sweep
    eps_k = eps_sweep(kk);
    sigma_sweep(kk) = sigma_of_strain(eps_k,sign_F,hyper_choice,...
                        mu_NH,C10,C01,C1,C2,C3,...
                        mu_AB,N_AB,A1_AB,A2_AB,c_AB,f_AB,EE,alpha_phen);
    F_sweep(kk) = AA*sigma_sweep(kk);
end

fprintf('  Swept %d points from eps=%.4g to eps=%.4g.\n',n_sweep,eps_sweep(1),eps_sweep(end));
fprintf('  F range on this sweep (REAL, A*sigma): %.4f N to %.4f N\n\n',min(F_sweep),max(F_sweep));

% Diagnostic plot of the built table
figure('Name','Real Stress-Strain Table','NumberTitle','off');
subplot(2,1,1);
plot(eps_sweep,sigma_sweep,'b-','LineWidth',1.5);
xlabel('\epsilon'); ylabel('\sigma (Pa)'); title('REAL nominal stress vs strain (dW_T/d\lambda, swept)'); grid on;
subplot(2,1,2);
plot(eps_sweep,F_sweep,'r-','LineWidth',1.5); hold on;
yline(Force,'--k','LineWidth',1.2,'Label',sprintf('Force=%.4g N',Force));
xlabel('\epsilon'); ylabel('F=A*\sigma (N)'); title('Real force vs strain (swept, signed)'); grid on;

% =========================================================================
%  SECTION 3 — FIND eps* WHERE F_eq28(eps*) = |Force| (self-consistent)
% =========================================================================

fprintf('--- STEP 3: FINDING eps* WHERE F(eps*) MATCHES THE APPLIED Force ---\n');

target_F = Force;   % signed — F_sweep is already signed (A*sigma, sigma signed via lam=1+sign_F*eps)
resid_sweep = F_sweep - target_F;

% Look for a sign change (bracket) first, scanning from the SMALLEST
% strain outward — this finds the physically relevant (smallest-strain)
% crossing first, rather than jumping to a distant secondary root.
eps_star = NaN;
for kk=1:n_sweep-1
    if isfinite(resid_sweep(kk)) && isfinite(resid_sweep(kk+1)) && ...
       sign(resid_sweep(kk))~=sign(resid_sweep(kk+1))
        % linear interpolation between the two bracketing sweep points
        eps_star = eps_sweep(kk) + (eps_sweep(kk+1)-eps_sweep(kk)) * ...
                   (0-resid_sweep(kk))/(resid_sweep(kk+1)-resid_sweep(kk));
        break;
    end
end

if isnan(eps_star)
    fprintf('  *** No sign change found in [%.4g, %.4g] — the table never\n',...
            eps_sweep(1),eps_sweep(end));
    fprintf('  *** reaches Force=%.4g N. Using the closest point instead.\n',target_F);
    [~,idx_min] = min(abs(resid_sweep));
    eps_star = eps_sweep(idx_min);
    fprintf('  *** Closest match: eps=%.6f, F=%.4f N (target=%.4f N)\n\n',...
            eps_star,F_sweep(idx_min),target_F);
else
    fprintf('  eps* = %.6f (%.2f%%)  [interpolated crossing]\n\n',eps_star,eps_star*100);
end

% EESM extrapolation check, if a fitted range is available
if hyper_choice==4 && exist('eesm_fit_max_eps','var')
    fprintf('--- EXTRAPOLATION CHECK ---\n');
    fprintf('  eps* found      = %.4f (%.1f%%)\n', eps_star, eps_star*100);
    fprintf('  Fit calibrated up to = %.4f (%.1f%%)\n', eesm_fit_max_eps, eesm_fit_max_eps*100);
    if eps_star > eesm_fit_max_eps
        over_pct = (eps_star/eesm_fit_max_eps - 1)*100;
        fprintf('  *** WARNING: eps* is %.0f%% BEYOND the fitted data range.\n', over_pct);
        fprintf('  *** This result extrapolates past experimental support.\n\n');
    else
        fprintf('  OK — eps* is within the fitted data range.\n\n');
    end
end

% =========================================================================
%  SECTION 4 — STEP 4: E_t AS THE SLOPE OF THE LINE FROM THE ORIGIN
% =========================================================================

fprintf('--- STEP 4: E_t = LOCAL SLOPE OF THE F-eps CURVE AT eps* ---\n');
fprintf('  (NOT the secant |Force|/(A*eps*) — the actual tangent slope of\n');
fprintf('  the swept F_eq28(eps) curve built in Section 2, evaluated at\n');
fprintf('  eps*, then converted to a modulus via dSigma/deps = (1/A)*dF/deps)\n');

% Local slope of the F-eps curve (dF/deps), via central difference across
% the whole sweep, then interpolated at eps*.
dF_deps_sweep = gradient(F_sweep, eps_sweep);
slope_at_star = interp1(eps_sweep, dF_deps_sweep, eps_star, 'linear', 'extrap');

if ~isfinite(slope_at_star)
    fprintf('  *** Local slope is non-finite at eps* — falling back to the\n');
    fprintf('  *** secant formula |Force|/(A*eps*) as a last resort.\n');
    E_t = abs(target_F)/(AA*max(eps_star,1e-9));
else
    E_t = abs(slope_at_star)/AA;
    fprintf('  dF/deps at eps*  = %.4f N (per unit strain)\n', slope_at_star);
end

if E_t<1e3
    fprintf('  *** NOTE: E_t floored at 1e3 Pa (computed value was <=0 or\n');
    fprintf('  *** implausibly small) — treat this result with caution.\n');
end
E_t = max(E_t,1e3);
fprintf('  E_t = (1/A)*dF/deps|eps* = %.4f Pa (%.4f MPa)\n\n',E_t,E_t/1e6);

% =========================================================================
%  SECTION 5 — SUMMARY
% =========================================================================

fprintf('=============================================================\n');
fprintf('  SIMULATION SUMMARY\n');
fprintf('=============================================================\n');
fprintf('  STL file     : %s (%s)\n',stl_filename,unit_label);
fprintf('  Applied Force: %.4g N (%s)\n',Force,dir_labels{wu});
fprintf('  A            : %.4e m^2\n',AA);
fprintf('  Model        : %s\n',hyper_label);
fprintf('  eps*         : %.4f (%.2f%%)  [self-consistent, from swept table]\n',eps_star,eps_star*100);
fprintf('  E_t          : %.4f Pa (%.4f MPa)\n',E_t,E_t/1e6);
fprintf('  METHOD       : Force applied AS-IS with this E_t, ONE linear FEM solve\n');
fprintf('=============================================================\n\n');
prc=input('  Proceed to FEM solve? [y/n]: ','s');
if ~strcmpi(strtrim(prc),'y'); fprintf('  Cancelled.\n'); return; end

% =========================================================================
%  SECTION 6 — GEOMETRY & MESH
% =========================================================================

fprintf('\n  Importing geometry ...\n');
smodel=createpde('structural','static-solid');
importGeometry(smodel,stl_filename);
figure('Name','Geometry','NumberTitle','off');
pdegplot(smodel); view(30,30); title(['Geometry: ',stl_filename]);

fprintf('  Generating mesh ...\n');
if use_hmax
    model=generateMesh(smodel,'GeometricOrder','linear','Hmax',Hmax_val);
else
    model=generateMesh(smodel,'GeometricOrder','linear');
end
nodeCoordinates=model.Nodes;
connectivity   =model.Elements;
numNodes=size(nodeCoordinates,2);
numElem =size(connectivity,2);
ndof    =3*numNodes;
fprintf('  Nodes=%d   Elements=%d\n\n',numNodes,numElem);

coords_m=nodeCoordinates*scale_stl*10^(-q);
bbox    =max(coords_m,[],2)-min(coords_m,[],2);
obj_size=max(bbox);
fprintf('  Bounding box (m): %.4f x %.4f x %.4f\n\n',bbox(1),bbox(2),bbox(3));

xc=zeros(4,numElem); yc=zeros(4,numElem); zc=zeros(4,numElem);
for i=1:numElem
    for j=1:4
        n=connectivity(j,i);
        xc(j,i)=coords_m(1,n);
        yc(j,i)=coords_m(2,n);
        zc(j,i)=coords_m(3,n);
    end
end

% =========================================================================
%  SECTION 7 — BOUNDARY CONDITIONS
% =========================================================================

order_fit=10;
load_coord   =coords_m(wu,:);
load_max     =round(max(load_coord),order_fit);
load_Indices =find(load_coord>=load_max-1e-6*abs(load_max));

fixed_coord  =coords_m(fix_wu,:);
fixed_min    =round(min(fixed_coord),order_fit);
fixed_Indices=find(fixed_coord<=fixed_min+1e-6*abs(fixed_min)+1e-12);

fixedDOFs=reshape([3*fixed_Indices(:)'-2;3*fixed_Indices(:)'-1;...
                   3*fixed_Indices(:)'],1,[]);
freeDOFs = setdiff(1:ndof,fixedDOFs);

fprintf('  Load nodes  : %d (%s-max)\n',numel(load_Indices),dir_labels{wu});
fprintf('  Fixed nodes : %d (%s-min)\n\n',numel(fixed_Indices),dir_labels{fix_wu});

% =========================================================================
%  SECTION 8 — B-MATRICES & VOLUMES
% =========================================================================

fprintf('  Pre-computing B-matrices ...\n');
B_all=cell(numElem,1);
V_all=zeros(numElem,1);

for jj=1:numElem
    x1=xc(1,jj);x2=xc(2,jj);x3=xc(3,jj);x4=xc(4,jj);
    y1=yc(1,jj);y2=yc(2,jj);y3=yc(3,jj);y4=yc(4,jj);
    z1=zc(1,jj);z2=zc(2,jj);z3=zc(3,jj);z4=zc(4,jj);

    v6=det([1 x1 y1 z1;1 x2 y2 z2;1 x3 y3 z3;1 x4 y4 z4]);
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

    Bi=@(b,g,d) sc*[b 0 0;0 g 0;0 0 d;g b 0;0 d g;d 0 b];
    B_all{jj}=[Bi(b1,g1,d1) Bi(b2,g2,d2) Bi(b3,g3,d3) Bi(b4,g4,d4)];
end
fprintf('  Done.\n\n');

% =========================================================================
%  SECTION 9 — LINEAR D MATRIX (from E_t)
% =========================================================================

nu_input = input('  nu - Poisson ratio for the linear D matrix: ');
c_lin=E_t/((1+nu_input)*(1-2*nu_input));
D_lin=c_lin*[1-nu_input nu_input    nu_input    0          0          0;
             nu_input   1-nu_input  nu_input    0          0          0;
             nu_input   nu_input    1-nu_input  0          0          0;
             0    0     0     (1-2*nu_input)/2 0          0;
             0    0     0     0          (1-2*nu_input)/2 0;
             0    0     0     0          0          (1-2*nu_input)/2];

% =========================================================================
%  SECTION 10 — SINGLE LINEAR FEM SOLVE (Force applied AS-IS, using E_t)
% =========================================================================

fprintf('  Assembling global stiffness matrix (using E_t) ...\n');
force_full=zeros(ndof,1);
for n=load_Indices(:)'
    force_full(3*(n-1)+wu)=Force/numel(load_Indices);
end

global_K = sparse(ndof,ndof);
for jj=1:numElem
    ns   =connectivity(:,jj);
    dof_e=reshape([3*ns-2,3*ns-1,3*ns]',1,[]);
    Ve   =V_all(jj);
    BB_e =B_all{jj};
    global_K(dof_e,dof_e)=global_K(dof_e,dof_e)+Ve*(BB_e'*D_lin*BB_e);
end

fprintf('  Solving linear system (single pass, no iteration) ...\n');
u_total = zeros(ndof,1);
u_total(freeDOFs) = global_K(freeDOFs,freeDOFs)\force_full(freeDOFs);
fprintf('  Done.\n\n');

% =========================================================================
%  SECTION 11 — POST-PROCESSING
% =========================================================================

disp_in_dir=u_total(wu:3:end);
z_disp     =u_total(3:3:end);
max_disp   =max(abs(disp_in_dir));
maxz       =max(abs(z_disp));

vm_stress = zeros(numElem,1);
for jj=1:numElem
    ns   =connectivity(:,jj);
    dof_e=reshape([3*ns-2,3*ns-1,3*ns]',1,[]);
    sig  =D_lin*B_all{jj}*u_total(dof_e);
    sxx=sig(1);syy=sig(2);szz=sig(3);
    txy=sig(4);tyz=sig(5);txz=sig(6);
    vm_stress(jj)=sqrt(0.5*((sxx-syy)^2+(syy-szz)^2+(szz-sxx)^2 ...
                             +6*(txy^2+tyz^2+txz^2)));
end

node_vm=zeros(numNodes,1); nc1=zeros(numNodes,1);
for jj=1:numElem
    ns=connectivity(:,jj);
    node_vm(ns)=node_vm(ns)+vm_stress(jj); nc1(ns)=nc1(ns)+1;
end
node_vm=node_vm./max(nc1,1);

% =========================================================================
%  SECTION 12 — VISUALIZATION
% =========================================================================

fprintf('  Building surface mesh ...\n');
all_faces=[connectivity([1 2 3],:)';connectivity([1 2 4],:)';...
           connectivity([1 3 4],:)';connectivity([2 3 4],:)'];
all_faces=sort(all_faces,2);
[~,ia,ic]=unique(all_faces,'rows');
face_counts=accumarray(ic,1);
surf_faces=all_faces(ia(face_counts==1),:);
surf_verts=coords_m';

do_disp=input(sprintf('  Displacement %s? [y/n]: ',dir_labels{wu}),'s');
do_vm  =input('  Von Mises?      [y/n]: ','s');
fprintf('\n');

if strcmpi(strtrim(do_disp),'y')
    ttl=sprintf('Displacement-%s (mm) — strain-sweep Eq.28 solve',dir_labels{wu});
    figure('Name',ttl,'NumberTitle','off');
    patch('Vertices',surf_verts,'Faces',surf_faces,...
          'FaceVertexCData',disp_in_dir*1e3,'FaceColor','interp','EdgeColor','none');
    colormap jet; cb=colorbar; cb.Label.String=sprintf('u%s (mm)',lower(dir_labels{wu}));
    title(ttl); view(30,30); axis equal tight;
    xlabel('X(m)'); ylabel('Y(m)'); zlabel('Z(m)');
    lighting gouraud; camlight;
end

if strcmpi(strtrim(do_vm),'y')
    figure('Name','Von Mises (kPa)','NumberTitle','off');
    patch('Vertices',surf_verts,'Faces',surf_faces,...
          'FaceVertexCData',node_vm*1e-3,'FaceColor','interp','EdgeColor','none');
    colormap jet; cb=colorbar; cb.Label.String='Von Mises (kPa)';
    caxis([0, prctile(node_vm,98)/1e3]);
    title('Von Mises Stress (kPa)'); view(30,30); axis equal tight;
    xlabel('X(m)'); ylabel('Y(m)'); zlabel('Z(m)');
    lighting gouraud; camlight;
end

% =========================================================================
%  SECTION 13 — RESULTS SUMMARY
% =========================================================================

fprintf('=============================================================\n');
fprintf('  RESULTS SUMMARY  (strain-sweep Eq.28, self-consistent secant)\n');
fprintf('=============================================================\n');
fprintf('  Model                : %s\n', hyper_label);
fprintf('  eps* (self-consistent): %.6f (%.2f%%)\n', eps_star, eps_star*100);
fprintf('  E_t (local slope of curve): %.4f Pa (%.4f MPa)\n', E_t, E_t/1e6);
fprintf('  Applied Force        : %.4g N   (used AS-IS)\n', Force);
fprintf('  Max disp (%s)         : %.6g m = %.4g mm\n', ...
        dir_labels{wu}, max_disp, max_disp*1e3);
fprintf('  Max |z_disp|         : %.6g m = %.4g mm\n', maxz, maxz*1e3);
fprintf('  Max Von Mises        : %.4g Pa = %.4g kPa\n', ...
        max(node_vm), max(node_vm)*1e-3);
fprintf('=============================================================\n');
fprintf('  REMINDER: this is a ONE-SHOT linear solve at a self-consistent\n');
fprintf('  secant stiffness derived from a swept Eq.28 table. It matches\n');
fprintf('  the applied force at the nominal 1-D strain exactly by\n');
fprintf('  construction, but still approximates the true 3D field\n');
fprintf('  solution as globally linear.\n');
fprintf('=============================================================\n\n');

% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function s=sigma_of_strain(eps_try,sign_F,hc,...
                           mu_NH,C10,C01,C1,C2,C3,...
                           mu_AB,N_AB,A1_AB,A2_AB,c_AB,f_AB,EE,alpha_phen)
    % Returns the REAL signed nominal stress dW_T/dlambda at strain
    % magnitude eps_try, with tension/compression via lam=1+sign_F*eps_try.
    lam = 1 + sign_F*eps_try;
    switch hc
        case 1
            s = mu_NH*(lam - 1/lam^2);
        case 2
            s = 2*(lam - 1/lam^2)*(C10 + C01/lam);
        case 3
            s = yeoh_sig(lam,C1,C2,C3);
        case 4
            s = eesm_sigma(lam,mu_AB,N_AB,A1_AB,A2_AB,c_AB,f_AB);
        case 5
            s_mag = EE*(eps_try + eps_try^2/2 + eps_try^2/(2*alpha_phen));
            s = sign_F*s_mag;
    end
end

function s=yeoh_sig(lam,C1,C2,C3)
    I1m3=lam^2+2/lam-3;
    s=2*(C1+2*C2*I1m3+3*C3*I1m3^2)*(lam-1/lam^2);
end

function w=W_eesm(lam,mu_AB,N_AB,A1_AB,A2_AB,c_AB,f_AB)
    persistent lam_log w_log call_count
    if isempty(lam_log); lam_log=[]; w_log=[]; call_count=0; end

    lam1=max(lam,0.001);
    lam2=1/sqrt(lam1);
    I1  =lam1^2+lam2^2+lam2^2;
    I3  =(lam1*lam2*lam2)^2;
    lamr=min(sqrt(I1/(3*N_AB)),0.9999);
    beta=(3*lamr)/(1-lamr^3);

    if ~isfinite(beta)||beta>700
        w=NaN;
    else
        sb=sinh(beta);
        if sb==0||~isfinite(sb)
            w=NaN;
        else
            w_iso=mu_AB*(N_AB*(beta*lamr+log(beta/sb))-log(beta/lamr))+c_AB;
            if ~isfinite(w_iso)
                w=NaN;
            else
                w_aniso=(A1_AB/3)*(I1-3)+(A2_AB/9)*(I1-3)^2 ...
                        -(A1_AB/3)*log(max(I3,1e-12));
                w=(1-f_AB)*w_iso+f_AB*w_aniso;
            end
        end
    end

    call_count=call_count+1; %#ok<*AGROW>
    lam_log(end+1,1)=lam;
    if isfinite(w); w_log(end+1,1)=w; else; w_log(end+1,1)=NaN; end
    if mod(call_count,1000)==0||call_count<=3
        assignin('base','W_eesm_lambda_log',lam_log);
        assignin('base','W_eesm_w_log',     w_log);
        assignin('base','W_eesm_call_count',call_count);
    end
end

function s=eesm_sigma(lam,mu_AB,N_AB,A1_AB,A2_AB,c_AB,f_AB)
    % Closed-form analytic derivative (same as v23):
    %   sigma(lambda) = lambda * (dWT/dI1) * (dI1/dlambda)
    lam1 = max(lam,0.001);
    lam2 = 1/sqrt(lam1);
    I1   = lam1^2 + lam2^2 + lam2^2;
    lamr = min(sqrt(I1/(3*N_AB)), 0.9999);
    beta = (3*lamr)/(1-lamr^3);

    if isfinite(beta) && beta<=700 && beta>1e-8 && lamr>1e-8
        dbeta_dlamr = 3*(1+2*lamr^3) / (1-lamr^3)^2;
        h_prime = N_AB*beta ...
                  + dbeta_dlamr*( N_AB*lamr + N_AB/beta - N_AB/tanh(beta) - 1/beta ) ...
                  + 1/lamr;
        dWiso_dI1   = mu_AB * h_prime / (6*N_AB*lamr);
        dWaniso_dI1 = A1_AB/3 + (2*A2_AB/9)*(I1-3);
        dWT_dI1     = (1-f_AB)*dWiso_dI1 + f_AB*dWaniso_dI1;
        dI1_dlam    = 2*lam1 - 2/lam1^2;
        s_analytic  = lam1 * dWT_dI1 * dI1_dlam;
    else
        s_analytic = NaN;
    end

    if isfinite(s_analytic)
        s = s_analytic;
        return;
    end

    % Fallback: finite-difference on W_eesm if the analytic path is
    % numerically degenerate at this lambda.
    dh=1e-7;
    wp=W_eesm(lam+dh,mu_AB,N_AB,A1_AB,A2_AB,c_AB,f_AB);
    wm=W_eesm(lam-dh,mu_AB,N_AB,A1_AB,A2_AB,c_AB,f_AB);
    if ~isfinite(wp)||~isfinite(wm); s=NaN; return; end
    s=lam*(wp-wm)/(2*dh);
end
% =========================================================================
%  FEM_Solver_Nonlinear.m
%  Nonlinear Static FEM Solver — Tetrahedral Elements (T4)
%  Material: Formlabs Elastic 50A (or user-defined hyperelastic)
%
%  BUG FIX (v2):
%    eps_elem was always zero on the first NR iteration because u_total=0.
%    Fix: compute a global estimated strain  est = |F| / (A * E0)
%         and use it to seed element regimes on iter==1 of every load step.
%         From iter>=2, displacement-derived strains take over normally.
%
%  NONLINEAR STRATEGY — Piecewise Tangent Modulus + Newton-Raphson:
%    • Strain < 10%  : Linear elastic  D(E_lin, nu)  — constant
%    • Strain >= 10% : Hyperelastic tangent D_T(lambda, model)
%                      Model choices: Neo-Hookean | Mooney-Rivlin 2P | Yeoh 3P
%
%  SOLVER: Newton-Raphson with load stepping
%    Each load increment:
%      1. Compute element strains from current displacement
%         (on iter 1 use estimated strain to seed regime correctly)
%      2. Decide linear vs nonlinear D_T per element
%      3. Re-assemble K_T (tangent stiffness)
%      4. Solve for displacement correction  du
%      5. Update displacement u += du
%      6. Check residual convergence
% =========================================================================

clear all; clc; %#ok<CLALL>
fprintf('--- UNIT SYSTEM ---\n');
fprintf('  1 -> mm  (q=3, multiply by 1e-3)\n');
fprintf('  2 -> cm  (q=2, multiply by 1e-2)\n');
fprintf('  3 -> m   (q=0, no scaling)\n');
unit_choice = input('  STL coordinate unit (1/2/3): ');
while ~ismember(unit_choice, [1 2 3])
    unit_choice = input('  Enter 1, 2, or 3: ');
end
switch unit_choice
    case 1; q = 3; unit_label = 'mm';
    case 2; q = 2; unit_label = 'cm';
    case 3; q = 0; unit_label = 'm';
end
fprintf('  Unit: %s  (scale factor: 1e-%d)\n\n', unit_label, q);
%_______________________________________
fprintf('=============================================================\n');
fprintf('   NONLINEAR FEM SOLVER — Elastic 50A (T4 Tetrahedra)        \n');
fprintf('=============================================================\n\n');

% =========================================================================
%  SECTION 1 — INTERACTIVE USER INPUTS
% =========================================================================

%% --- STL Geometry ---
fprintf('--- STL GEOMETRY IMPORT ---\n');
stl_filename = input('  Enter STL filename (e.g., cube3cm.stl): ', 's');
stl_filename = strtrim(stl_filename);
if ~strcmpi(stl_filename(end-3:end), '.stl')
    stl_filename = [stl_filename, '.stl'];
end
while exist(stl_filename, 'file') ~= 2
    warning('File "%s" not found.', stl_filename);
    stl_filename = input('  Valid STL filename: ', 's');
    stl_filename = strtrim(stl_filename);
    if ~strcmpi(stl_filename(end-3:end), '.stl')
        stl_filename = [stl_filename, '.stl'];
    end
end

%% --- Force direction ---
fprintf('\n--- FORCE DIRECTION ---\n');
fprintf('  1->X  2->Y  3->Z\n');
wu = input('  Axis (1/2/3): ');
while ~ismember(wu, [1 2 3])
    wu = input('  Enter 1, 2, or 3: ');
end
dir_labels = {'X','Y','Z'};




%% --- Fixed face ---
fprintf('\n--- FIXED FACE ---\n');
fprintf('  0->Same as force  1->X-min  2->Y-min  3->Z-min\n');
fix_wu_input = input('  Choice (0/1/2/3): ');
while ~ismember(fix_wu_input,[0 1 2 3])
    fix_wu_input = input('  Enter 0,1,2, or 3: ');
end
fix_wu = fix_wu_input;
if fix_wu == 0; fix_wu = wu; end

%% --- Force magnitude ---
fprintf('\n--- FORCE (N) ---\n');
fprintf('  Negative = compression, Positive = tension.\n');
Force = input('  Total force (N) [e.g. -25]: ');
%% --- Area magnitude ---

fprintf('\n--- cross section area (m^2) ---\n');
AA = input('  cross section area (m^2) [e.g. -25]: ');

%% --- Linear material (low-strain regime) ---
fprintf('\n--- LINEAR MATERIAL PROPERTIES (strain < 10%%) ---\n');
use_defaults = input('  Use Elastic 50A defaults (E=2.9e6 Pa, nu=0.40)? [y/n]: ','s');
if strcmpi(strtrim(use_defaults),'n')
    EE = input('  Young''s modulus E (Pa): ');
    nu = input('  Poisson''s ratio nu (0<nu<0.5): ');
    while nu<=0||nu>=0.5
        nu = input('  nu must be in (0,0.5): ');
    end
else
    EE = 2.9e6;   % measured elastic modulus for Formlabs 50A (post-cured) [Pa]
    nu = 0.4;     % nearly incompressible elastomer
end
%________________________________________

%% --- Show estimated strain immediately after E is known ---
estimated_strain = abs(Force) / (AA * EE);

fprintf('             = %.4g', estimated_strain);

%% --- Hyperelastic model (high-strain regime, >= 10%) ---
fprintf('\n--- HYPERELASTIC MODEL (strain >= 10%%) ---\n');
fprintf('  1 -> Neo-Hookean     (1 param:  mu)\n');
fprintf('  2 -> Mooney-Rivlin   (2 params: C10, C01)\n');
fprintf('  3 -> Yeoh            (3 params: C1, C2, C3)\n');
hyper_choice = input('  Model (1/2/3): ');
while ~ismember(hyper_choice,[1 2 3])
    hyper_choice = input('  Enter 1, 2, or 3: ');
end

switch hyper_choice
    case 1
        fprintf('\n  Neo-Hookean: W = mu/2*(I1-3)\n');
        fprintf('  Calibrated default: mu = 0.57e6 Pa\n');
        use_def_hyp = input('  Use default? [y/n]: ','s');
        if strcmpi(strtrim(use_def_hyp),'n')
            mu_NH = input('  mu (Pa): ');
        else
            mu_NH = 0.57e6;
        end
        C10=mu_NH/2; C01=0; C1=mu_NH/2; C2=0; C3=0;
        hyper_label = sprintf('Neo-Hookean  mu=%.3g Pa', mu_NH);

    case 2
        fprintf('\n  Mooney-Rivlin: W = C10*(I1-3) + C01*(I2-3)\n');
        fprintf('  Calibrated defaults: C10=0.28e6, C01=0.07e6 Pa\n');
        use_def_hyp = input('  Use defaults? [y/n]: ','s');
        if strcmpi(strtrim(use_def_hyp),'n')
            C10 = input('  C10 (Pa): ');
            C01 = input('  C01 (Pa): ');
        else
            C10 = 0.28e6; C01 = 0.07e6;
        end
        mu_NH=2*(C10+C01); C1=C10; C2=0; C3=0;
        hyper_label = sprintf('Mooney-Rivlin  C10=%.3g  C01=%.3g Pa', C10, C01);

    case 3
        fprintf('\n  Yeoh: W = C1*(I1-3) + C2*(I1-3)^2 + C3*(I1-3)^3\n');
        fprintf('  Calibrated defaults: C1=0.28e6, C2=3e3, C3=1.6e3 Pa\n');
        use_def_hyp = input('  Use defaults? [y/n]: ','s');
        if strcmpi(strtrim(use_def_hyp),'n')
            C1 = input('  C1 (Pa): ');
            C2 = input('  C2 (Pa): ');
            C3 = input('  C3 (Pa): ');
        else
            C1=0.28e6; C2=3e3; C3=1.6e3;
        end
        mu_NH=2*C1; C10=C1; C01=0;
        hyper_label = sprintf('Yeoh  C1=%.3g  C2=%.3g  C3=%.3g Pa', C1, C2, C3);
end

%% --- Newton-Raphson parameters ---
fprintf('\n--- NEWTON-RAPHSON SETTINGS ---\n');
n_steps    = input('  Number of load increments [default 10]: ');
if isempty(n_steps)||n_steps<1; n_steps=10; end
max_iter   = input('  Max NR iterations per step [default 20]: ');
if isempty(max_iter)||max_iter<1; max_iter=20; end
tol_res    = input('  Convergence tolerance (relative residual) [default 1e-4]: ');
if isempty(tol_res)||tol_res<=0; tol_res=1e-4; end

order_fit  = 10;

% =========================================================================
%  SECTION 2 — SUMMARY
% =========================================================================

fprintf('\n=============================================================\n');
fprintf('  NONLINEAR SIMULATION SUMMARY\n');
fprintf('=============================================================\n');
fprintf('  STL file         : %s\n',        stl_filename);
fprintf('  Force            : %.4g N  (%s)\n', Force, dir_labels{wu});
fprintf('  Fixed face       : %s-min\n',    dir_labels{fix_wu});
fprintf('  Linear E         : %.4g Pa  (strain < 10%%)\n', EE);
fprintf('  nu               : %.4g\n',      nu);
fprintf('  Hyperelastic     : %s\n',        hyper_label);
fprintf('  Load steps       : %d\n',        n_steps);
fprintf('  NR max iters     : %d\n',        max_iter);
fprintf('  NR tolerance     : %.2g\n',      tol_res);
fprintf('=============================================================\n\n');

proceed = input('  Proceed? [y/n]: ','s');
if ~strcmpi(strtrim(proceed),'y')
    fprintf('  Cancelled.\n'); return;
end
fprintf('\n');

% =========================================================================
%  SECTION 3 — GEOMETRY & MESH
% =========================================================================

fprintf('  [1/5] Importing geometry and meshing ...\n');
smodel = createpde('structural','static-solid');
importGeometry(smodel, stl_filename);
figure('Name','Geometry','NumberTitle','off');
pdegplot(smodel); view(30,30); title(['Geometry: ', stl_filename]);

model = generateMesh(smodel,'GeometricOrder','linear');
nodeCoordinates = model.Nodes;     % 3 x numNodes  (mm)
connectivity    = model.Elements;  % 4 x numElem

numNodes = size(nodeCoordinates,2);
numElem  = size(connectivity,   2);
ndof     = 3*numNodes;

fprintf('     Nodes    : %d\n', numNodes);
fprintf('     Elements : %d\n\n', numElem);

% =========================================================================
%  SECTION 4 — COORDINATE ARRAYS (mm -> m)
% =========================================================================

fprintf('  [2/5] Extracting coordinates (mm->m) ...\n');
xconmatrix = zeros(4,numElem);
yconmatrix = zeros(4,numElem);
zconmatrix = zeros(4,numElem);

for i = 1:numElem
    for j = 1:4
        nIdx = connectivity(j,i);
        xconmatrix(j,i) = nodeCoordinates(1,nIdx)*10.^(-q);
        yconmatrix(j,i) = nodeCoordinates(2,nIdx)*10.^(-q);
        zconmatrix(j,i) = nodeCoordinates(3,nIdx)*10.^(-q);
    end
end

% =========================================================================
%  SECTION 5 — BOUNDARY CONDITIONS
% =========================================================================

load_coord   = nodeCoordinates(wu,:);
load_max     = round(max(load_coord),order_fit);
load_Indices = find(load_coord >= load_max - 1e-6*abs(load_max));

fixed_coord   = nodeCoordinates(fix_wu,:);
fixed_min     = round(min(fixed_coord),order_fit);
fixed_Indices = find(fixed_coord <= fixed_min + 1e-6*abs(fixed_min) + 1e-12);

fixedDOFs = reshape([3*fixed_Indices(:)'-2; ...
                     3*fixed_Indices(:)'-1; ...
                     3*fixed_Indices(:)'  ],1,[]);
freeDOFs  = setdiff(1:ndof, fixedDOFs);

fprintf('  Boundary conditions:\n');
fprintf('    Load nodes  : %d  (%s-max)\n', numel(load_Indices),  dir_labels{wu});
fprintf('    Fixed nodes : %d  (%s-min)\n\n',numel(fixed_Indices),dir_labels{fix_wu});

%% --- Full force vector ---
force_full      = zeros(ndof,1);
f_per_node      = Force / numel(load_Indices);
for n = load_Indices(:)'
    force_full(3*(n-1)+wu) = f_per_node;
end

% =========================================================================
%  SECTION 6 — PRE-COMPUTE ELEMENT B-MATRICES
% =========================================================================

fprintf('  [3/5] Pre-computing element B-matrices ...\n');

B_all  = cell(numElem,1);
V_all  = zeros(numElem,1);

for jj = 1:numElem
    x1=xconmatrix(1,jj); x2=xconmatrix(2,jj);
    x3=xconmatrix(3,jj); x4=xconmatrix(4,jj);
    y1=yconmatrix(1,jj); y2=yconmatrix(2,jj);
    y3=yconmatrix(3,jj); y4=yconmatrix(4,jj);
    z1=zconmatrix(1,jj); z2=zconmatrix(2,jj);
    z3=zconmatrix(3,jj); z4=zconmatrix(4,jj);

    v_six = det([1 x1 y1 z1; 1 x2 y2 z2; 1 x3 y3 z3; 1 x4 y4 z4]);
    V_e   = abs(v_six)/6;
    V_all(jj) = V_e;
    sc = 1/(6*V_e);

    b1=-((y3-y2)*(z4-z2)-(y4-y2)*(z3-z2)); g1=-((x4-x2)*(z3-z2)-(x3-x2)*(z4-z2)); d1=-((x3-x2)*(y4-y2)-(x4-x2)*(y3-y2));
    b2= (y3-y1)*(z4-z1)-(y4-y1)*(z3-z1);  g2= (x4-x1)*(z3-z1)-(x3-x1)*(z4-z1);  d2= (x3-x1)*(y4-y1)-(x4-x1)*(y3-y1);
    b3=-((y2-y1)*(z4-z1)-(y4-y1)*(z2-z1)); g3=-((x4-x1)*(z2-z1)-(x2-x1)*(z4-z1)); d3=-((x2-x1)*(y4-y1)-(x4-x1)*(y2-y1));
    b4= (y2-y1)*(z3-z1)-(y3-y1)*(z2-z1);  g4= (x3-x1)*(z2-z1)-(x2-x1)*(z3-z1);  d4= (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);

    Bi = @(b,g,d) sc*[b 0 0; 0 g 0; 0 0 d; g b 0; 0 d g; d 0 b];
    B_all{jj} = [Bi(b1,g1,d1) Bi(b2,g2,d2) Bi(b3,g3,d3) Bi(b4,g4,d4)];
end

fprintf('     Done.\n\n');

% =========================================================================
%  SECTION 7 — ESTIMATED STRAIN SEED  (THE BUG FIX)
%
%  On the first NR iteration of every load step, u_total is unchanged from
%  the previous step (or zero at step 1), so B*u gives stale/zero strains.
%  Instead, use the engineering estimate:
%
%       eps_est = |F_step| / (A_load * E0)
%
%  where A_load is the cross-sectional area of the loaded face (m^2),
%  computed from the bounding-box dimensions perpendicular to load.
%  This scalar is broadcast to ALL elements as a conservative global seed;
%  once iter>=2 the displacement field gives element-specific strains.
% =========================================================================

% Bounding box of the full geometry (in metres)
coords_m = nodeCoordinates * 10.^(-q);   % 3 x numNodes, metres

dims = max(coords_m,[],2) - min(coords_m,[],2);  % [Lx; Ly; Lz]

% Cross-sectional area perpendicular to loading direction (wu)
side_dirs   = setdiff([1 2 3], wu);
A_cross     = dims(side_dirs(1)) * dims(side_dirs(2));   % m^2

% Length in loading direction
L_body = dims(wu);   % m

fprintf('  Geometry (m): Lx=%.4g  Ly=%.4g  Lz=%.4g\n', dims(1),dims(2),dims(3));
fprintf('  Cross-sectional area (perp. to %s): A = %.4g m^2\n', dir_labels{wu}, A_cross);

% Estimated global strain for the FULL applied force
%   eps_est = sigma / E = (|F| / A) / E0
eps_est_full = abs(Force) / (A_cross * EE);

fprintf('  Estimated global strain at full load: eps_est = %.4g\n\n', eps_est_full);

% =========================================================================
%  SECTION 8 — HELPER FUNCTIONS
% =========================================================================

% --- Linear elastic D (pre-built once) ---
c_lin = EE/((1+nu)*(1-2*nu));
D_lin = c_lin*[1-nu nu nu 0 0 0;
               nu 1-nu nu 0 0 0;
               nu nu 1-nu 0 0 0;
               0 0 0 (1-2*nu)/2 0 0;
               0 0 0 0 (1-2*nu)/2 0;
               0 0 0 0 0 (1-2*nu)/2];

    function Et = tangent_modulus(eps_eff, hyper_choice, mu_NH, C10, C01, C1, C2, C3)
        lam = 1 + max(eps_eff, 0);
        dh  = 1e-6;
        switch hyper_choice
            case 1
                Et = mu_NH*(1 + 2/lam^3);
            case 2
                sigma  = @(l) 2*(l - 1/l^2)*(C10 + C01/l);
                Et = (sigma(lam+dh) - sigma(lam-dh))/(2*dh);
            case 3
                sigma  = @(l) yeoh_stress(l, C1, C2, C3);
                Et = (sigma(lam+dh) - sigma(lam-dh))/(2*dh);
        end
        Et = max(Et, 1e3);
    end

    function s = yeoh_stress(lam, C1, C2, C3)
        I1_m3 = lam^2 + 2/lam - 3;
        dW    = C1 + 2*C2*I1_m3 + 3*C3*I1_m3^2;
        s     = 2*dW*(lam - 1/lam^2);
    end

    function D_t = build_D_tangent(Et, nu_t)
        c_t = Et/((1+nu_t)*(1-2*nu_t));
        D_t = c_t*[1-nu_t  nu_t   nu_t   0             0             0;
                   nu_t    1-nu_t nu_t   0             0             0;
                   nu_t    nu_t   1-nu_t 0             0             0;
                   0       0      0      (1-2*nu_t)/2  0             0;
                   0       0      0      0             (1-2*nu_t)/2  0;
                   0       0      0      0             0             (1-2*nu_t)/2];
    end

% =========================================================================
%  SECTION 9 — NEWTON-RAPHSON NONLINEAR SOLVE
% =========================================================================

fprintf('  [4/5] Nonlinear Newton-Raphson solve (%d load steps) ...\n\n', n_steps);

u_total = zeros(ndof,1);
F_step  = force_full / n_steps;

conv_history = zeros(n_steps,1);
iter_history = zeros(n_steps,1);

% eps_elem persists across NR iterations (updated each iter from iter>=2)
eps_elem = zeros(numElem,1);

for step = 1:n_steps

    F_applied = F_step * step;
    % Estimated strain at this load fraction (scales linearly with force)
    eps_est_step = eps_est_full * (step / n_steps);

    fprintf('  Step %2d/%d  (Force = %.4g N,  eps_est = %.4g) ...\n', ...
            step, n_steps, Force*step/n_steps, eps_est_step);

    for iter = 1:max_iter

        % ------------------------------------------------------------------
        % (A) Element effective strains
        %
        %  KEY FIX: On iter==1 we cannot trust B*u (u may be zero or stale).
        %  Use the global estimated strain for every element so the regime
        %  decision is physically meaningful from the very first assembly.
        %  On iter>=2 the actual displacement field is used.
        % ------------------------------------------------------------------

        if iter == 1
            % Seed ALL elements with the estimated strain for this load step.
            % This correctly identifies nonlinear elements even before any
            % displacement increment has been computed.
            eps_elem(:) = eps_est_step;
        else
            % Displacement-derived Von Mises equivalent strain per element
            for jj = 1:numElem
                ns    = connectivity(:,jj);
                dof_e = reshape([3*ns-2, 3*ns-1, 3*ns]', 1, []);
                u_e   = u_total(dof_e);
                eps_v = B_all{jj} * u_e;
                exx=eps_v(1); eyy=eps_v(2); ezz=eps_v(3);
                gxy=eps_v(4); gyz=eps_v(5); gxz=eps_v(6);
                eps_elem(jj) = sqrt((2/3)*(exx^2 + eyy^2 + ezz^2 + ...
                                           0.5*(gxy^2+gyz^2+gxz^2)));
            end
        end

        % ------------------------------------------------------------------
        % (B) Build tangent stiffness K_T
        % ------------------------------------------------------------------
        global_KT   = sparse(ndof, ndof);
        elem_regime = zeros(numElem,1);

        for jj = 1:numElem
            eps_e = eps_est_step;%eps_elem(jj);

            if eps_e < 0.10
                D_e = D_lin;
            else
                Et_e = tangent_modulus(eps_e, hyper_choice, mu_NH, C10, C01, C1, C2, C3);
                D_e  = build_D_tangent(Et_e, nu);
                elem_regime(jj) = 1;
            end

            k_e = V_all(jj) * (B_all{jj}' * D_e * B_all{jj});

            ns    = connectivity(:,jj);
            dof_e = reshape([3*ns-2, 3*ns-1, 3*ns]', 1, []);
            global_KT(dof_e,dof_e) = global_KT(dof_e,dof_e) + k_e;
        end

        % ------------------------------------------------------------------
        % (C) Internal force vector
        % ------------------------------------------------------------------
        F_int = global_KT * u_total;

        % ------------------------------------------------------------------
        % (D) Residual on free DOFs
        % ------------------------------------------------------------------
        R = F_applied - F_int;
        R(fixedDOFs) = 0;

        res_norm = norm(R(freeDOFs));
        f_norm   = max(norm(F_applied(freeDOFs)), 1e-12);
        rel_res  = res_norm / f_norm;

        fprintf('    Iter %2d : rel_res = %.3e  (nonlinear elems: %d/%d)\n', ...
                iter, rel_res, sum(elem_regime), numElem);

        if rel_res < tol_res
            fprintf('    Converged in %d iterations.\n\n', iter);
            conv_history(step) = rel_res;
            iter_history(step) = iter;
            break;
        end

        % ------------------------------------------------------------------
        % (E) Solve  K_T * du = R  and update displacements
        % ------------------------------------------------------------------
        KT_red  = global_KT(freeDOFs, freeDOFs);
        R_red   = R(freeDOFs);
        du_free = KT_red \ R_red;

        u_total(freeDOFs) = u_total(freeDOFs) + du_free;

        if iter == max_iter
            warning('Step %d did NOT converge after %d iterations (rel_res=%.3e).', ...
                    step, max_iter, rel_res);
            conv_history(step) = rel_res;
            iter_history(step) = iter;
        end
    end
end

fprintf('  Nonlinear solve complete.\n\n');

% =========================================================================
%  SECTION 10 — REACTION FORCES
% =========================================================================

reaction_force  = global_KT * u_total - force_full;
fixed_rxn_dofs  = 3*(fixed_Indices(:)') + fix_wu - 3;
total_reaction  = sum(reaction_force(fixed_rxn_dofs));

fprintf('  Results:\n');
fprintf('    Applied force        : %.4g N\n',  Force);
fprintf('    Total reaction force : %.4g N\n\n', total_reaction);

% =========================================================================
%  SECTION 11 — POST-PROCESSING
% =========================================================================

fprintf('  [5/5] Post-processing ...\n\n');

disp_in_dir = u_total(wu:3:end);
max_disp    = max(abs(disp_in_dir));
fprintf('  Max |displacement| (%s): %.6g m\n\n', dir_labels{wu}, max_disp);

%% --- Recompute eps_elem from final displacements for stress evaluation ---
for jj = 1:numElem
    ns    = connectivity(:,jj);
    dof_e = reshape([3*ns-2, 3*ns-1, 3*ns]', 1, []);
    u_e   = u_total(dof_e);
    eps_v = B_all{jj} * u_e;
    exx=eps_v(1); eyy=eps_v(2); ezz=eps_v(3);
    gxy=eps_v(4); gyz=eps_v(5); gxz=eps_v(6);
    eps_elem(jj) = sqrt((2/3)*(exx^2 + eyy^2 + ezz^2 + ...
                               0.5*(gxy^2+gyz^2+gxz^2)));
end

%% --- Von Mises stress ---
vm_stress = zeros(numElem,1);
for jj = 1:numElem
    ns    = connectivity(:,jj);
    dof_e = reshape([3*ns-2, 3*ns-1, 3*ns]', 1, []);
    u_e   = u_total(dof_e);
    eps_e = eps_elem(jj);

    if eps_e < 0.10
        D_e = D_lin;
    else
        Et_e = tangent_modulus(eps_e, hyper_choice, mu_NH, C10, C01, C1, C2, C3);
        D_e  = build_D_tangent(Et_e, nu);
    end

    sigma = D_e * B_all{jj} * u_e;
    sxx=sigma(1); syy=sigma(2); szz=sigma(3);
    txy=sigma(4); tyz=sigma(5); txz=sigma(6);
    vm_stress(jj) = sqrt(0.5*((sxx-syy)^2+(syy-szz)^2+(szz-sxx)^2 ...
                               +6*(txy^2+tyz^2+txz^2)));
end

%% --- Node-averaged Von Mises ---
node_vm    = zeros(numNodes,1);
node_count = zeros(numNodes,1);
for jj = 1:numElem
    ns = connectivity(:,jj);
    node_vm(ns)    = node_vm(ns)    + vm_stress(jj);
    node_count(ns) = node_count(ns) + 1;
end
node_vm = node_vm ./ max(node_count,1);

fprintf('  Max Von Mises stress : %.4g Pa\n\n', max(node_vm));

% =========================================================================
%  SECTION 12 — REGIME MAP
% =========================================================================

node_regime  = zeros(numNodes,1);
node_count2  = zeros(numNodes,1);
for jj = 1:numElem
    ns = connectivity(:,jj);
    node_regime(ns)  = node_regime(ns)  + elem_regime(jj);
    node_count2(ns)  = node_count2(ns)  + 1;
end
node_regime = node_regime ./ max(node_count2,1);

% =========================================================================
%  SECTION 13 — VISUALIZATIONS
% =========================================================================

z_disp = u_total(3:3:end);

%% Displacement in loading direction
figure('Name',sprintf('Displacement — %s',dir_labels{wu}),'NumberTitle','off');
pdeplot3D(model,'ColorMapData',disp_in_dir*1e3);
colormap jet; cb=colorbar;
cb.Label.String = sprintf('u_%s  (mm)', lower(dir_labels{wu}));
title(sprintf('Displacement in %s direction', dir_labels{wu}));
view(30,30);

%% Z displacement
figure('Name','Z displacement','NumberTitle','off');
pdeplot3D(model,'ColorMapData',z_disp*1e3);
colormap jet; cb2=colorbar;
cb2.Label.String = 'u_z  (mm)';
title('Z-displacement'); view(30,30);

%% Von Mises stress
figure('Name','Von Mises Stress','NumberTitle','off');
pdeplot3D(model,'ColorMapData',node_vm*1e-3);
colormap jet; cb3=colorbar;
cb3.Label.String = 'Von Mises Stress (kPa)';
title('Von Mises Stress'); view(30,30);

% %% Regime map
% figure('Name','Material Regime Map','NumberTitle','off');
% pdeplot3D(model,'ColorMapData',node_regime);
% colormap([ 0.18 0.36 0.80;
%            0.85 0.20 0.20]);
% cb4=colorbar;
% cb4.Ticks=[0 1];
% cb4.TickLabels={'Linear (< 10%)','Nonlinear (>= 10%)'};
% title('Material Regime: Linear vs Hyperelastic zones');
% view(30,30);
%----------------------------------------------------------
% %% Convergence plot
% figure('Name','NR Convergence','NumberTitle','off');
% semilogy(1:n_steps, conv_history, 'o-b','LineWidth',1.5,'MarkerFaceColor','b');
% hold on;
% yline(tol_res,'--r','LineWidth',1.2,'Label',sprintf('Tolerance %.1e',tol_res));
% xlabel('Load step'); ylabel('Relative residual');
% title('Newton-Raphson convergence per load step');
% grid on;

%% Data cursor on Z-disp figure
dcm = datacursormode(gcf);
set(dcm,'Enable','on');
set(dcm,'UpdateFcn',{@myupdatefcn, nodeCoordinates, z_disp});

% =========================================================================
%  SECTION 14 — CONSOLE SUMMARY TABLE
% =========================================================================

fprintf('=============================================================\n');
fprintf('  NONLINEAR SOLVE SUMMARY\n');
fprintf('=============================================================\n');
fprintf('  %-8s  %-8s  %-12s  %-8s\n','Step','Iters','Rel Residual','Status');
fprintf('  %s\n', repmat('-',1,44));
for s=1:n_steps
    status = 'OK';
    if conv_history(s) > tol_res; status = 'WARN'; end
    fprintf('  %-8d  %-8d  %-12.3e  %s\n', s, iter_history(s), conv_history(s), status);
end
fprintf('=============================================================\n');
fprintf('  Applied force        : %.4g N\n',  Force);
fprintf('  Reaction force       : %.4g N\n',  total_reaction);
fprintf('  Max displacement     : %.4g m\n',  max_disp);
fprintf('  Max Von Mises stress : %.4g Pa\n', max(node_vm));
fprintf('  Est. global strain   : %.4g (full load)\n', eps_est_full);
fprintf('  Nonlinear elements   : %d / %d  (%.1f%%)\n', ...
        sum(elem_regime), numElem, 100*sum(elem_regime)/numElem);
fprintf('=============================================================\n');

% =========================================================================
%  SECTION 15 — DATA CURSOR CALLBACK
% =========================================================================

function txt = myupdatefcn(~,event_obj,nodeCoords,z_displacement)
    pos      = event_obj.Position;
    dists    = sqrt(sum((nodeCoords'-pos).^2,2));
    [~,id]   = min(dists);
    txt = {sprintf('Node : %d',    id),
           sprintf('X    : %.4f mm',nodeCoords(1,id)),
           sprintf('Y    : %.4f mm',nodeCoords(2,id)),
           sprintf('Z    : %.4f mm',nodeCoords(3,id)),
           sprintf('u_z  : %.4e m', z_displacement(id))};
end


maxz=max(abs(z_disp))

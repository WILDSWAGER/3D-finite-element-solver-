% =========================================================================
%  FEM_Solver_Fixed.m
%  Linear Static FEM Solver — Tetrahedral Elements (T4/CST)
%  Description:
%    Solves a linear-elastic static problem on a solid tetrahedral mesh
%    imported from an STL file.

% =========================================================================

clear all; clc; %#ok<CLALL>

fprintf('=============================================================\n');
fprintf('       LINEAR STATIC FEM SOLVER — Tetrahedral Elements       \n');
fprintf('=============================================================\n\n');

% =========================================================================
%  SECTION 1 — INTERACTIVE USER INPUTS
% =========================================================================

%% --- Force direction ---
fprintf('--- FORCE DIRECTION ---\n');
fprintf('  Select the axis along which the force is applied:\n');
fprintf('    1  ->  X axis\n');
fprintf('    2  ->  Y axis\n');
fprintf('    3  ->  Z axis\n');
wu = input('  Enter choice (1 / 2 / 3): ');
while ~ismember(wu, [1 2 3])
    warning('Invalid input. Please enter 1, 2, or 3.');
    wu = input('  Enter choice (1 / 2 / 3): ');
end
dir_labels = {'X', 'Y', 'Z'};
fprintf('  Force direction : %s axis (wu = %d)\n\n', dir_labels{wu}, wu);

%% --- Fixed face direction ---
fprintf('--- FIXED (CLAMPED) FACE ---\n');
fprintf('  Select the axis whose MINIMUM face will be fully fixed (all DOFs = 0):\n');
fprintf('    1  ->  X-min face\n');
fprintf('    2  ->  Y-min face\n');
fprintf('    3  ->  Z-min face\n');
fprintf('    0  ->  Same axis as force (standard uniaxial setup)\n');
fix_wu_input = input('  Enter choice (0 / 1 / 2 / 3): ');
while ~ismember(fix_wu_input, [0 1 2 3])
    warning('Invalid input. Please enter 0, 1, 2, or 3.');
    fix_wu_input = input('  Enter choice (0 / 1 / 2 / 3): ');
end
if fix_wu_input == 0
    fix_wu = wu;
else
    fix_wu = fix_wu_input;
end
fprintf('  Fixed face      : %s-min face (fix_wu = %d)\n\n', dir_labels{fix_wu}, fix_wu);

%% --- Force magnitude ---
fprintf('--- FORCE MAGNITUDE ---\n');
fprintf('  Enter total applied force in Newtons.\n');
fprintf('  Convention: Negative = compression,  Positive = tension.\n');
Force = input('  Force (N) [e.g. -25]: ');
fprintf('  Force           : %.4g N\n\n', Force);

%% --- Material properties ---
fprintf('--- MATERIAL PROPERTIES ---\n');
EE_default  = 3e6;
nu_default  = 0.46;
use_defaults = input('  Use default material (E = 3e6 Pa, nu = 0.46)? [y/n]: ', 's');

if strcmpi(strtrim(use_defaults), 'n')
    EE = input('  Young''s modulus E (Pa): ');
    nu = input('  Poisson''s ratio  nu (0 < nu < 0.5): ');
    while nu <= 0 || nu >= 0.5
        warning('Poisson''s ratio must be in the range (0, 0.5).');
        nu = input('  Poisson''s ratio  nu (0 < nu < 0.5): ');
    end
else
    EE = EE_default;
    nu = nu_default;
end
fprintf('  E               : %.4g Pa\n', EE);
fprintf('  nu              : %.4g\n\n', nu);

%% --- Precision for BC coordinate rounding ---
order_fit = 10;   % decimal places for min/max coordinate comparison

% =========================================================================
%  SECTION 2 — SUMMARY
% =========================================================================

fprintf('=============================================================\n');
fprintf('  SIMULATION SUMMARY\n');
fprintf('=============================================================\n');
fprintf('  STL file         : cube3cm.stl\n');
fprintf('  Force            : %.4g N  along  %s axis\n', Force, dir_labels{wu});
fprintf('  Load face        : %s-max face\n', dir_labels{wu});
fprintf('  Fixed face       : %s-min face  (all DOFs clamped)\n', dir_labels{fix_wu});
fprintf('  Young''s modulus  : %.4g Pa\n', EE);
fprintf('  Poisson''s ratio  : %.4g\n', nu);
fprintf('=============================================================\n\n');

proceed = input('  Proceed with analysis? [y/n]: ', 's');
if ~strcmpi(strtrim(proceed), 'y')
    fprintf('  Analysis cancelled by user.\n');
    return;
end
fprintf('\n');

% =========================================================================
%  SECTION 3 — GEOMETRY IMPORT & MESH GENERATION
% =========================================================================

fprintf('  [1/6] Importing geometry ...\n');
smodel = createpde('structural', 'static-solid');
importGeometry(smodel, 'cube3cm.stl');
figure('Name','Geometry','NumberTitle','off');
pdegplot(smodel);
view(30, 30);
title('Imported Geometry');

fprintf('  [2/6] Generating mesh ...\n');
model = generateMesh(smodel, 'GeometricOrder', 'linear');

nodeCoordinates = model.Nodes;      % 3 x numNodes  (mm, from STL)
connectivity    = model.Elements;   % 4 x numElem

numNodes = size(nodeCoordinates, 2);
numElem  = size(connectivity,    2);

fprintf('         Nodes    : %d\n', numNodes);
fprintf('         Elements : %d\n\n', numElem);

% =========================================================================
%  SECTION 4 — COORDINATE EXTRACTION  (mm -> m)
% =========================================================================

fprintf('  [3/6] Extracting element coordinates (mm -> m) ...\n');

xconmatrix = zeros(4, numElem);
yconmatrix = zeros(4, numElem);
zconmatrix = zeros(4, numElem);

for i = 1:numElem
    for j = 1:4
        nodeIdx = connectivity(j, i);
        xconmatrix(j, i) = nodeCoordinates(1, nodeIdx) * 1e-3;   % mm -> m
        yconmatrix(j, i) = nodeCoordinates(2, nodeIdx) * 1e-3;
        zconmatrix(j, i) = nodeCoordinates(3, nodeIdx) * 1e-3;
    end
end

% =========================================================================
%  SECTION 5 — MATERIAL CONSTITUTIVE MATRIX  [D]
% =========================================================================

fprintf('  [4/6] Building material matrix [D] ...\n');

c = EE / ((1 + nu) * (1 - 2*nu));

D_matrix = c * ...
    [1-nu   nu    nu    0              0              0;
     nu     1-nu  nu    0              0              0;
     nu     nu    1-nu  0              0              0;
     0      0     0     (1-2*nu)/2     0              0;
     0      0     0     0              (1-2*nu)/2     0;
     0      0     0     0              0              (1-2*nu)/2];

% =========================================================================
%  SECTION 6 — PARALLEL ELEMENT STIFFNESS MATRICES  [k_e]
% =========================================================================

fprintf('  [5/6] Assembling element stiffness matrices (parallel) ...\n');

if isempty(gcp('nocreate'))
    parpool('local', feature('numcores'));
end

k_matrix = cell(numElem, 1);

parfor jj = 1:numElem

    x1 = xconmatrix(1,jj);  x2 = xconmatrix(2,jj);
    x3 = xconmatrix(3,jj);  x4 = xconmatrix(4,jj);

    y1 = yconmatrix(1,jj);  y2 = yconmatrix(2,jj);
    y3 = yconmatrix(3,jj);  y4 = yconmatrix(4,jj);

    z1 = zconmatrix(1,jj);  z2 = zconmatrix(2,jj);
    z3 = zconmatrix(3,jj);  z4 = zconmatrix(4,jj);

    % ------------------------------------------------------------------
    %  6V = det of coordinate matrix (signed volume × 6)
    %  NOTE: sign depends on node ordering convention.
    %        We use abs() for the volume scalar, and divide B by v_six
    %        (with its sign) so that B'DB remains positive definite.
    % ------------------------------------------------------------------
    v_six = det([1 x1 y1 z1;
                 1 x2 y2 z2;
                 1 x3 y3 z3;
                 1 x4 y4 z4]);

    % Volume of element (must be positive)
    V_elem = abs(v_six) / 6;

    % ------------------------------------------------------------------
    %  Shape-function gradient coefficients
    %  Standard T4 formulation: N_i = (alpha_i + beta_i*x + gamma_i*y + delta_i*z) / (6V)
    %  The coefficients below are the SIGNED cofactors of the coord matrix.
    % ------------------------------------------------------------------

    % --- Node 1 ---
    beta_1  = -( (y3-y2)*(z4-z2) - (y4-y2)*(z3-z2) );
    gama_1  = -( (x4-x2)*(z3-z2) - (x3-x2)*(z4-z2) );
    del1    = -( (x3-x2)*(y4-y2) - (x4-x2)*(y3-y2) );

    % --- Node 2 ---
    beta_2  =  (y3-y1)*(z4-z1) - (y4-y1)*(z3-z1);
    gama_2  =  (x4-x1)*(z3-z1) - (x3-x1)*(z4-z1);
    del2    =  (x3-x1)*(y4-y1) - (x4-x1)*(y3-y1);

    % --- Node 3 ---
    beta_3  = -( (y2-y1)*(z4-z1) - (y4-y1)*(z2-z1) );
    gama_3  = -( (x4-x1)*(z2-z1) - (x2-x1)*(z4-z1) );
    del3    = -( (x2-x1)*(y4-y1) - (x4-x1)*(y2-y1) );

    % --- Node 4 ---
    beta_4  =  (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1);
    gama_4  =  (x3-x1)*(z2-z1) - (x2-x1)*(z3-z1);
    del4    =  (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1);

    % ------------------------------------------------------------------
    %  Strain-displacement matrix B (6 x 12)
    %  Strain vector: [eps_xx, eps_yy, eps_zz, gam_xy, gam_yz, gam_xz]'
    %
    %  For node i:
    %    B_i = 1/(6V) * [ beta_i   0        0      ]
    %                   [ 0        gama_i   0      ]
    %                   [ 0        0        del_i  ]
    %                   [ gama_i   beta_i   0      ]
    %                   [ 0        del_i    gama_i ]
    %                   [ del_i    0        beta_i ]
    %
    %  NOTE: Rows 5 & 6 were swapped/incorrect in the original code.
    %        Row 5 = gamma_yz  (dv/dz + dw/dy) -> uses del_i  and gama_i on v,w cols
    %        Row 6 = gamma_xz  (du/dz + dw/dx) -> uses del_i  and beta_i on u,w cols
    % ------------------------------------------------------------------

    scale = 1 / (6 * V_elem);   % always positive

    B1 = scale * [beta_1   0        0     ;
                  0        gama_1   0     ;
                  0        0        del1  ;
                  gama_1   beta_1   0     ;
                  0        del1     gama_1;
                  del1     0        beta_1];

    B2 = scale * [beta_2   0        0     ;
                  0        gama_2   0     ;
                  0        0        del2  ;
                  gama_2   beta_2   0     ;
                  0        del2     gama_2;
                  del2     0        beta_2];

    B3 = scale * [beta_3   0        0     ;
                  0        gama_3   0     ;
                  0        0        del3  ;
                  gama_3   beta_3   0     ;
                  0        del3     gama_3;
                  del3     0        beta_3];

    B4 = scale * [beta_4   0        0     ;
                  0        gama_4   0     ;
                  0        0        del4  ;
                  gama_4   beta_4   0     ;
                  0        del4     gama_4;
                  del4     0        beta_4];

    BB = [B1  B2  B3  B4];   % 6 x 12

    % Element stiffness:  k_e = V * B' * D * B
    k_matrix{jj} = V_elem * (BB' * D_matrix * BB);

end

% =========================================================================
%  SECTION 7 — GLOBAL STIFFNESS ASSEMBLY  [K]
% =========================================================================

fprintf('  [6/6] Assembling global stiffness matrix [K] ...\n\n');

ndof     = 3 * numNodes;
global_k = sparse(ndof, ndof);

for w = 1:numElem
    nodes = connectivity(:, w);
    % DOF mapping: node n -> DOFs [3n-2, 3n-1, 3n]  (u, v, w)
    dof = reshape([3*nodes-2,  3*nodes-1,  3*nodes]', 1, []);
    global_k(dof, dof) = global_k(dof, dof) + k_matrix{w};
end

% =========================================================================
%  SECTION 8 — BOUNDARY CONDITIONS
% =========================================================================

%% --- Load nodes: MAX face along force direction ---
load_coord   = nodeCoordinates(wu, :);
load_max     = round(max(load_coord), order_fit);
load_Indices = find(load_coord >= load_max - 1e-6 * abs(load_max));

%% --- Fixed nodes: MIN face along fix_wu direction ---
fixed_coord   = nodeCoordinates(fix_wu, :);
fixed_min     = round(min(fixed_coord), order_fit);
fixed_Indices = find(fixed_coord <= fixed_min + 1e-6 * abs(fixed_min) + 1e-12);

fprintf('  Boundary conditions:\n');
fprintf('    Load nodes   : %d  (on %s-max face)\n', numel(load_Indices),  dir_labels{wu});
fprintf('    Fixed nodes  : %d  (on %s-min face)\n\n', numel(fixed_Indices), dir_labels{fix_wu});

% =========================================================================
%  SECTION 9 — FORCE VECTOR  {F}
% =========================================================================

force_v    = zeros(ndof, 1);
f_per_node = Force / numel(load_Indices);

for n = load_Indices(:)'
    % FIX: correct DOF index for direction wu at node n
    % DOF layout: node n -> [3(n-1)+1, 3(n-1)+2, 3(n-1)+3] = [3n-2, 3n-1, 3n]
    % Direction wu=1 -> 3(n-1)+1 = 3n-2
    %           wu=2 -> 3(n-1)+2 = 3n-1
    %           wu=3 -> 3(n-1)+3 = 3n
    dof_idx = 3*(n-1) + wu;   % FIXED (was: 3*n + wu - 3, same but clearer)
    force_v(dof_idx) = f_per_node;
end

% =========================================================================
%  SECTION 10 — APPLY FIXED BOUNDARY CONDITIONS & SOLVE
% =========================================================================

%% Collect ALL 3 DOFs for every fixed node
fixedDOFs = reshape([3*fixed_Indices(:)'-2; ...
                     3*fixed_Indices(:)'-1; ...
                     3*fixed_Indices(:)'  ], 1, []);

%% Save full system for reaction-force computation
total_global_k = global_k;
total_force_v  = force_v;

%% Reduce system (eliminate fixed DOFs)
freeDOFs = setdiff(1:ndof, fixedDOFs);

global_k_red = global_k(freeDOFs, freeDOFs);
force_v_red  = force_v(freeDOFs);

%% Solve  [K_red]{u_free} = {F_free}
fprintf('  Solving system (%d free DOFs) ...\n\n', numel(freeDOFs));

free_disp = global_k_red \ force_v_red;

% =========================================================================
%  SECTION 11 — RECONSTRUCT FULL DISPLACEMENT VECTOR
% =========================================================================

total_disp = zeros(ndof, 1);
total_disp(freeDOFs) = free_disp;

% =========================================================================
%  SECTION 12 — REACTION FORCES
% =========================================================================

reaction_force = total_global_k * total_disp - total_force_v;

% DOFs at fixed nodes in the fix_wu direction
fixed_reaction_dofs = 3*(fixed_Indices(:)') + fix_wu - 3;   % = 3(n-1) + fix_wu
total_reaction      = sum(reaction_force(fixed_reaction_dofs));

fprintf('  Results:\n');
fprintf('    Applied force       : %.4g N\n',  Force);
fprintf('    Total reaction force: %.4g N\n\n', total_reaction);

% =========================================================================
%  SECTION 13 — POST-PROCESSING & VISUALIZATION
% =========================================================================

%% Displacement in the force direction (wu)
%  DOF wu for node n is at index 3*(n-1)+wu.
%  Extracting every 3rd starting at wu gives nodes 1,2,...,numNodes in order.
disp_in_dir = total_disp(wu : 3 : end);   % numNodes x 1

max_disp = max(abs(disp_in_dir));
fprintf('  Max displacement (|%s|): %.6g m\n\n', dir_labels{wu}, max_disp);

%% --- Plot: Displacement in force direction ---
figure('Name', sprintf('Displacement — %s direction', dir_labels{wu}), ...
       'NumberTitle', 'off');
pdeplot3D(model, 'ColorMapData', disp_in_dir * 1e3);   % convert to mm
colormap jet;
cb = colorbar;
cb.Label.String = sprintf('Displacement u_%s  (mm)', lower(dir_labels{wu}));
title(sprintf('Displacement in %s direction', dir_labels{wu}));
view(30, 30);

%% --- Plot: Z displacement (always shown for reference) ---
z_disp = total_disp(3:3:end);   % every 3rd DOF starting at index 3 = uz

figure('Name', 'Displacement — Z direction (reference)', 'NumberTitle','off');
pdeplot3D(model, 'ColorMapData', z_disp * 1e3);
colormap jet;
cb2 = colorbar;
cb2.Label.String = 'Displacement u_z  (mm)';
title('Z-displacement (reference)');
view(30, 30);

%% --- Von Mises stress (optional post-processing) ---
fprintf('  Computing Von Mises stress ...\n');
vm_stress = zeros(numElem, 1);

for jj = 1:numElem
    nodes_e = connectivity(:, jj);
    dof_e   = reshape([3*nodes_e-2, 3*nodes_e-1, 3*nodes_e]', 1, []);
    u_e     = total_disp(dof_e);   % 12 x 1

    x1 = xconmatrix(1,jj); x2 = xconmatrix(2,jj);
    x3 = xconmatrix(3,jj); x4 = xconmatrix(4,jj);
    y1 = yconmatrix(1,jj); y2 = yconmatrix(2,jj);
    y3 = yconmatrix(3,jj); y4 = yconmatrix(4,jj);
    z1 = zconmatrix(1,jj); z2 = zconmatrix(2,jj);
    z3 = zconmatrix(3,jj); z4 = zconmatrix(4,jj);

    v_six_e = det([1 x1 y1 z1; 1 x2 y2 z2; 1 x3 y3 z3; 1 x4 y4 z4]);
    V_e     = abs(v_six_e) / 6;
    sc      = 1 / (6 * V_e);

    beta_1e  = -( (y3-y2)*(z4-z2) - (y4-y2)*(z3-z2) );
    gama_1e  = -( (x4-x2)*(z3-z2) - (x3-x2)*(z4-z2) );
    del1e    = -( (x3-x2)*(y4-y2) - (x4-x2)*(y3-y2) );
    beta_2e  =  (y3-y1)*(z4-z1) - (y4-y1)*(z3-z1);
    gama_2e  =  (x4-x1)*(z3-z1) - (x3-x1)*(z4-z1);
    del2e    =  (x3-x1)*(y4-y1) - (x4-x1)*(y3-y1);
    beta_3e  = -( (y2-y1)*(z4-z1) - (y4-y1)*(z2-z1) );
    gama_3e  = -( (x4-x1)*(z2-z1) - (x2-x1)*(z4-z1) );
    del3e    = -( (x2-x1)*(y4-y1) - (x4-x1)*(y2-y1) );
    beta_4e  =  (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1);
    gama_4e  =  (x3-x1)*(z2-z1) - (x2-x1)*(z3-z1);
    del4e    =  (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1);

    B1e = sc*[beta_1e 0 0; 0 gama_1e 0; 0 0 del1e; gama_1e beta_1e 0; 0 del1e gama_1e; del1e 0 beta_1e];
    B2e = sc*[beta_2e 0 0; 0 gama_2e 0; 0 0 del2e; gama_2e beta_2e 0; 0 del2e gama_2e; del2e 0 beta_2e];
    B3e = sc*[beta_3e 0 0; 0 gama_3e 0; 0 0 del3e; gama_3e beta_3e 0; 0 del3e gama_3e; del3e 0 beta_3e];
    B4e = sc*[beta_4e 0 0; 0 gama_4e 0; 0 0 del4e; gama_4e beta_4e 0; 0 del4e gama_4e; del4e 0 beta_4e];

    BBe   = [B1e B2e B3e B4e];
    sigma = D_matrix * BBe * u_e;   % [sxx syy szz txy tyz txz]'

    sxx = sigma(1); syy = sigma(2); szz = sigma(3);
    txy = sigma(4); tyz = sigma(5); txz = sigma(6);

    vm_stress(jj) = sqrt(0.5*((sxx-syy)^2 + (syy-szz)^2 + (szz-sxx)^2 ...
                               + 6*(txy^2 + tyz^2 + txz^2)));
end

%% Map element stress to nodes (simple averaging)
node_vm = zeros(numNodes, 1);
node_count = zeros(numNodes, 1);
for jj = 1:numElem
    ns = connectivity(:, jj);
    node_vm(ns)    = node_vm(ns)    + vm_stress(jj);
    node_count(ns) = node_count(ns) + 1;
end
node_vm = node_vm ./ max(node_count, 1);

figure('Name','Von Mises Stress','NumberTitle','off');
pdeplot3D(model, 'ColorMapData', node_vm * 1e-3);   % display in kPa
colormap jet;
cb3 = colorbar;
cb3.Label.String = 'Von Mises Stress (kPa)';
title('Von Mises Stress');
view(30, 30);

fprintf('  Max Von Mises stress: %.4g Pa\n\n', max(node_vm));

%% --- Enable data cursor for node inspection (Z-displacement plot) ---
dcm = datacursormode(gcf);
set(dcm, 'Enable', 'on');
set(dcm, 'UpdateFcn', {@myupdatefcn, nodeCoordinates, z_disp});

% =========================================================================
%  SECTION 14 — DATA CURSOR CALLBACK
% =========================================================================

function txt = myupdatefcn(~, event_obj, nodeCoords, z_displacement)
% MYUPDATEFCN  Custom data-cursor tooltip for FEM displacement plots.

    pos = event_obj.Position;

    % Find nearest node
    dists   = sqrt(sum((nodeCoords' - pos) .^ 2, 2));
    [~, id] = min(dists);

    txt = {
        sprintf('Node      : %d',     id),
        sprintf('X         : %.4f mm', nodeCoords(1, id)),
        sprintf('Y         : %.4f mm', nodeCoords(2, id)),
        sprintf('Z         : %.4f mm', nodeCoords(3, id)),
        sprintf('u_z       : %.4e m',  z_displacement(id))
    };
end

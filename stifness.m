

clear all; clc;

%% ================= INPUTS =================
Force = -25;                 % Total force (compression < 0)
EE = 2.8e6;                   % Young's modulus
nu = 0.38;                    % Poisson's ratio
wu = 3;                       % Load direction (1=x,2=y,3=z)
order_fit = 20;                % Decimal rounding

%% ================= IMPORT STL =================
smodel = createpde('structural','static-solid');
importGeometry(smodel,'cube3cm.stl');
pdegplot(smodel); view(30,30);

%% ================= MESH =================
model = generateMesh(smodel,'GeometricOrder','linear');
nodeCoordinates = model.Nodes;        % 3 x N
connectivity    = model.Elements;     % 4 x Nel

numNodes = size(nodeCoordinates,2);
numElem  = size(connectivity,2);

%% ================= PARALLEL POOL =================
if isempty(gcp('nocreate'))
    parpool('local',feature('numcores'));
end

%% ================= COORDINATE MATRICES =================
xconmatrix = zeros(4,numElem);
yconmatrix = zeros(4,numElem);
zconmatrix = zeros(4,numElem);

for i = 1:numElem
    for j = 1:4
        hh = connectivity(j,i);
        xconmatrix(j,i) = nodeCoordinates(1,hh)*1e-3;%conver unit of mm(acording my stl file) to m
        yconmatrix(j,i) = nodeCoordinates(2,hh)*1e-3;
        zconmatrix(j,i) = nodeCoordinates(3,hh)*1e-3;
    end
end

%% ================= MATERIAL MATRIX =================
D_matrix = (EE/((1+nu)*(1-2*nu))) * ...
 [1-nu nu nu 0 0 0;
  nu 1-nu nu 0 0 0;
  nu nu 1-nu 0 0 0;
  0 0 0 (1-2*nu)/2 0 0;
  0 0 0 0 (1-2*nu)/2 0;
  0 0 0 0 0 (1-2*nu)/2];

%% ================= ELEMENT STIFFNESS =================
k_matrix = cell(numElem,1);

parfor jj = 1:numElem

    x = xconmatrix(:,jj); 
    y = yconmatrix(:,jj); 
    z = zconmatrix(:,jj);

    x1=x(1); x2=x(2); x3=x(3); x4=x(4);
    y1=y(1); y2=y(2); y3=y(3); y4=y(4);
    z1=z(1); z2=z(2); z3=z(3); z4=z(4);

    alpha_1 = x2*(y3*z4-y4*z3)+y2*(x4*z3-x3*z4)+z2*(x3*y4-x4*y3);
    beta_1  = -((y3*z4-z3*y4)+y2*(z3-z4)+z2*(y4-y3));
    gama_1  = (x3*z4-z3*x4)+(z3-z4)*x2+(x4-x3)*z2;
    del1    = -((x3*y4-y3*x4)+x2*(y3-y4)+y2*(x4-x3));

    alpha_2 = -(x1*(y3*z4-y4*z3)+y1*(x4*z3-x3*z4)+z1*(x3*y4-x4*y3));
    beta_2  = ((y3*z4-z3*y4)+y1*(z3-z4)+z1*(y4-y3));
    gama_2  = -((x3*z4-z3*x4)+(z3-z4)*x1+(x4-x3)*z1);
    del2    = (x3*y4-y3*x4)+x1*(y3-y4)+y1*(x4-x3);

    alpha_3 = x1*(y2*z4-y4*z2)+y1*(x4*z2-x2*z4)+z1*(x2*y4-x4*y2);
    beta_3  = -((y2*z4-z2*y4)+y1*(z2-z4)+z1*(y4-y2));
    gama_3  = (x2*z4-z2*x4)+(z2-z4)*x1+(x4-x2)*z1;
    del3    = -((x2*y4-y2*x4)+x1*(y2-y4)+y1*(x4-x2));

    alpha_4 = -(x1*(y2*z3-y3*z2)+y1*(x3*z2-x2*z3)+z1*(x2*y3-y2*x3));
    beta_4  = ((y2*z3-z2*y3)+y1*(z2-z3)+z1*(y3-y2));
    gama_4  = -((x2*z3-z2*x3)+(z2-z3)*x1+(x3-x2)*z1);
    del4    = ((x2*y3-y2*x3)+x1*(y2-y3)+y1*(x3-x2));

    v_six = det([1 x1 y1 z1;
                 1 x2 y2 z2;
                 1 x3 y3 z3;
                 1 x4 y4 z4]);

    B1=(1/v_six)*[beta_1 0 0;0 gama_1 0;0 0 del1;gama_1 beta_1 0;0 del1 gama_1;del1 0 beta_1];
    B2=(1/v_six)*[beta_2 0 0;0 gama_2 0;0 0 del2;gama_2 beta_2 0;0 del2 gama_2;del2 0 beta_2];
    B3=(1/v_six)*[beta_3 0 0;0 gama_3 0;0 0 del3;gama_3 beta_3 0;0 del3 gama_3;del3 0 beta_3];
    B4=(1/v_six)*[beta_4 0 0;0 gama_4 0;0 0 del4;gama_4 beta_4 0;0 del4 gama_4;del4 0 beta_4];

    BB = [B1 B2 B3 B4];
    k_matrix{jj} = (v_six/6)*(BB'*D_matrix*BB);
end

%% ================= GLOBAL STIFFNESS =================
ndof = 3*numNodes;
global_k = sparse(ndof,ndof);

for w = 1:numElem
    nodes = connectivity(:,w);
    dof = reshape([3*nodes-2 3*nodes-1 3*nodes]',1,[]);
    global_k(dof,dof) = global_k(dof,dof) + k_matrix{w};
end

%% ================= BOUNDARY CONDITIONS =================
coord = nodeCoordinates(wu,:);
fixed_Indices = find(coord <= round(min(coord),order_fit));
load_Indices  = find(coord >= round(max(coord),order_fit));

%% ================= FORCE VECTOR =================
force_v = zeros(ndof,1);
fnode = Force/numel(load_Indices);
for n = load_Indices'
    force_v(3*n+wu-3) = fnode;
end

%% ================= FIXED DOF =================
fixedpoits = reshape([3*fixed_Indices-2;3*fixed_Indices-1;3*fixed_Indices],1,[]);

%% ================= SOLUTION =================
total_global_k = global_k;
total_force_v  = force_v;

global_k(fixedpoits,:) = [];
global_k(:,fixedpoits) = [];
force_v(fixedpoits)    = [];

free_force_v = force_v;

%% obtain invers of global matrix
identity_matrix = eye(size(global_k));
invers_k = global_k \ identity_matrix;

%% obtain total displacement vector 
free_dis = invers_k * free_force_v;

%% ================= RECONSTRUCT =================
total_dis = zeros(ndof,1);
freeDOF = setdiff(1:ndof,fixedpoits);
total_dis(freeDOF) = free_dis;

%% ================= REACTION FORCE =================
reaction_force = total_global_k*total_dis - total_force_v;

%% ================= POST =================
z_dis = total_dis(3:3:end);
pdeplot3D(model,'ColorMapData',z_dis.*100);
colormap jet; colorbar;


%%%%%%%%%%%%%%FIND DISPLACMENT OF ESPATIAL NODE%%%%%%%%%%%%%%%%
z_dis = total_dis(3:3:end);

pdeplot3D(model,'ColorMapData',z_dis.*100);

colormap jet;
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%




pdeplot3D(model,'ColorMapData',z_dis.*100);

colormap jet;
colorbar;

%% Enable data cursor
dcm = datacursormode(gcf);

set(dcm,'Enable','on');

set(dcm,'UpdateFcn',{@myupdatefcn,nodeCoordinates,z_dis});

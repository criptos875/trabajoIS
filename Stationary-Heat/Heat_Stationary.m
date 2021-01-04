
close all;
clear;

inp_file_name='Sphere';

%% Parameters
T_exterior = 100;
T_interior = 0;
%% load mesh
mesh=ofem.mesh;
mesh.load_from_inp(inp_file_name);


%% Information abut the radious of the sphere
co = double(permute(mesh.co,[3,1,2]));

boundary1_index = unique(mesh.el(mesh.bd{2,1}{2,1}));
R_exterior = max(dot(co(boundary1_index,:)',co(boundary1_index,:)'));
R_exterior = sqrt(R_exterior);

boundary2_index = unique(mesh.el(mesh.bd{2,2}{2,1}));
R_interior = max(dot(co(boundary2_index,:)',co(boundary2_index,:)'));
R_interior = sqrt(R_interior);

%% Plot boundaries

 figure
 trimesh(mesh.el(mesh.bd{2,1}{2,1},:),co(:,1),co(:,2),co(:,3));
 hold on
 trimesh(mesh.el(mesh.bd{2,2}{2,1},:),co(:,1),co(:,2),co(:,3));
 title('Boundaries')
 xlabel('X')
 ylabel('Y')
 zlabel('Z')
 hold off
 

%% Plot parts
figure
part2_index = unique(mesh.el(mesh.parts{3,2}));
scatter3(co(part2_index,1),co(part2_index,2),co(part2_index,3));
hold on
part1_index = unique(mesh.el(mesh.parts{3,1}));
scatter3(co(part1_index,1),co(part1_index,2),co(part1_index,3));
 title('Parts')
 xlabel('X')
 ylabel('Y')
 zlabel('Z')
 hold off
 

%% define function space discretization
fe=ofem.finiteelement.P1;

%% define equation type in oFEM
eq=ofem.elliptic(mesh,fe,ofem.gaussianquadrature(mesh,fe));


%% stiffness matrix and boundary conditions
 opt.S = 1;
 opt.f = 0;
 opt.dirichlet{1} = struct('idx',1,'f',T_exterior);
 opt.dirichlet{2} = struct('idx',2,'f',T_interior);


[asm,info,~]=eq.assemble(opt);


%% total assembly and scalar material incorporation
b    = asm.b;
S = asm.S;
DOFs = asm.DOFs;

%% Solution
u = full(asm.dirichlet);
u(DOFs) = S(DOFs,DOFs) \b(DOFs);

%% Plot of the solution
figure
scatter3(co(:,1),co(:,2),co(:,3),[],u);
colorbar

   title('Temperature in the sphere')
   xlabel('X')
   ylabel('Y')
   zlabel('Z')
%% Analytic solution
alpha = (T_exterior-T_interior)/(1/R_exterior-1/R_interior);
beta = T_interior-alpha/R_interior;
% Auxiliary variables to represent radius
R_e = linspace(R_interior,R_exterior,1000);
R_i = linspace(0,R_interior);
% Analytic solution
T_e = alpha./R_e+beta;
T_i = T_interior.*ones(size(R_i));
%% Plot solution in terms of the distance to the center
figure
R = sqrt(co(:,1).*co(:,1)+co(:,2).*co(:,2)+co(:,3).*co(:,3));
scatter(R,u,30,u);
hold on
% Exterior analytical solution
plot(R_e,T_e);
% Interior analytical solution
plot(R_i,T_i);
grid;
colorbar;
xlabel("Radius");
ylabel("Temperature");
ylabel("Temperature in the sphere")





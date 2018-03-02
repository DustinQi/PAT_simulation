clear all
clc
% run $TOASTDIR/mtoast2_install.m if required
% mesh_file = 'validation000.msh';
% mesh_file = 'testDistance_2mm.msh';

[filename, pathname] = uigetfile( ...
    {'*.msh'},...
    '请选择mesh文件（不可多选）', ...
    'MultiSelect', 'off');
mesh_file = strcat(pathname,filename);
mesh=toastMesh(mesh_file,'gmsh');
%mesh.Display
%mesh.Write('ball_inner_toast.msh');%write out the mesh in Toast format  

% The internal structure is not obvious from the surface display of the mesh, but each 
% element contains information about the sub-volume it belongs to via its region index.
num_el = mesh.ElementCount;
num_nd = mesh.NodeCount;


%----------------------------------------------------------------------------
% 挑选两个有限元代替两个小球作为target,（这里以“ball.msh”为例）
%[vtx,idx,eltp] = mesh.Data;% vtx：每个格点的坐标，idx：每个有限元对应的四个格点
% 这里选取了第1509个格点和第3232个格点，它们的坐标分别为p1:（0,0）和p2:（0.3404,1.0276,0.3859）
% [p1_x,p1_y]=find(idx==1509);
% [p2_x,p2_y]=find(idx==3232);
% 选取包含有这两个格点的有限元306和3217，
%----------------------------------------------------------------------------


reg = mesh.Region;
region_kind = unique(reg);

for i=1:length(region_kind)
    eval(['region',num2str(i), '=', 'find(reg == region_kind(i))',';']);
 end

target1 = find(reg == region_kind(2));
target2 = find(reg == region_kind(3));
% target3 = find(reg == region_kind(4));

% target = find(reg == region_kind(2));

for i=1:length(region_kind)
    temp = find(reg == region_kind(i));
    eval(['num_region', num2str(i), '=', 'length(temp)', ';']);
end
clear temp i

% some parameters
refind = 1.4;   % refractive index
c0 = 0.3;       % speed of light in vacuum [mm/ps]
cm = c0/refind; % speed of light in the medium [mm/ps]
mua_bkg = 0.01; % background absorption [1/mm]
% mua_bkg = 0; % for 3D validation, background absorption = 0  [1/mm]
mus_bkg =1 ;    % background scattering [1/mm];
itrmax = 100;   % CG iteration limit
tolCG = 1e-6;   % convergence criterion

ref = ones(num_el,1)*refind;
mua = ones(num_el,1)*mua_bkg;
mua_0 = zeros(num_el,1);
mus = ones(num_el,1)*mus_bkg;

mua(target1) = mua_bkg*80; % target's mua[1/mm]
mua(target2) = mua_bkg*75; 
% mua(target3) = mua_bkg*30; 
% mus(target1) = mus_bkg*0.54; % target's mus[1/mm]
% mus(target2) = mus_bkg*0.54; % target's mus[1/mm]
% mus(target3) = mus_bkg*0.54; % target's mus[1/mm]

% mua(target1) = 0.49; % target's mua[1/mm]
% mua(target2) = 0.49; 
% mus(target1) = mus_bkg*0.54; % target's mus[1/mm]
% mus(target2) = mus_bkg*0.54; % target's mus[1/mm]

% mua(target) = 0.49; 
% mus(target) = mus_bkg*0.24; 

% % visualise the regions in a cross section
% grd = [256 256 256];
% basis = toastBasis(mesh,grd);
% elref = basis.GridElref;
% regim = zeros(size(elref));
% for i=1:length(elref)
%   el = elref(i);
%   if el>0
%     regim(i) = reg(el);
%   end
% end
% regim = reshape(regim,grd);
% imagesc(squeeze(regim(:,128,:)),[min(reg)-1,max(reg)]);
% axis equal tight; colorbar

%define the position of sources and detectors
num_detector = 128;
num_source =64;

% %48个平面光源（上）
% Q1 = makeCartCircle(7, 24, [0,0]);
% Q2 = makeCartCircle(5, 16, [0,0]);
% Q3 = makeCartCircle(3, 8, [0,0]);
% Q=[Q1';Q2';Q3'];
% Qz=ones(48,1).*(6);
% Q=[Q,Qz];
% clear Q1 Q2 Q3 Qz

% %48个平面光源（侧）
% Qx = ones(48,1).*(8);
% Qy = [-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8];
% Qy = [Qy,Qy,Qy];
% Qy=Qy';
% Qz = [ones(16,1).*(-5);zeros(16,1);zeros(16,1).*5];
% Q=[Qx,Qy,Qz];
% clear Qx Qy Qz

Q=1000*makeSphere_dustin(0.02,num_source,[0,0,0]);
Q=Q';
M=1000*makeSphere_dustin(0.02,num_detector,[0,0,0]);
M=M';
mesh.SetQM(Q,M);
qvec = mesh.Qvec ('Neumann', 'Gaussian', 2);
mvec = mesh.Mvec ('Gaussian', 2);

%run forward solver
K = dotSysmat (mesh,mua,mus,ref,'EL');
Phi = K\qvec;
K_0 = dotSysmat(mesh,mua_0,mus,ref,'EL');
Phi_0 = K_0\qvec;
%forward solver running time

% mesh-->grid(may depends on the size of k-wave)
% M = 64;
M1 = 30;
M2 = 30;
M3 = 30;
grd = [M1 M2 M3];
basis = toastBasis(mesh, grd);

% get absorb
PPP_mua = zeros(M1*M2*M3,num_detector);
for i=1:num_source
    PPP_mua(:,i) = basis.Map('M->B',Phi(:,i));
end

PPP_mua_0 = zeros(M1*M2*M3,num_detector);
for i=1:num_source
    PPP_mua_0(:,i) = basis.Map('M->B',Phi_0(:,i));
end

absorb = zeros(M1,M2,M3,num_source);
for i=1:num_source
    absorb(:,:,:,i) = reshape((PPP_mua_0(:,i)-PPP_mua(:,i)),[M1,M2,M3]);
end
absorb = sum(absorb,4);

% Visualisation of 3D matrix
xs = 1:2:M1;
ys = 1:2:M2;
zs = 1:2:M3;
absorb_copy = absorb;
% absorb_copy(absorb_copy<0.4)=NaN;
absorb_copy(absorb_copy==0)=NaN;
%absorb_copy(1:M/2,1:M,2/M:M)=NaN;
absorb_copy(1:M1,1:M2/2,1:M3)=NaN;
h = slice(absorb_copy,xs,ys,zs);
set(h,'FaceColor','interp','EdgeColor','none')
camproj perspective
box on
view(-90,0)
axis([1 M1 1 M2 1 M3])
colormap jet
colorbar

clearvars -except absorb

% %return position of maximun of the absorb
% [Amax, indmax] = max(absorb(:));
% [max_x, max_y, max_z] = ind2sub(size(absorb), indmax);
% max_absorb = [max_x,max_y,max_z];

% clearvars -except absorb PPP num_el num_nd region_kind max_absorb grd mesh_file M f_time



% % Visualisation of 3D matrix
% xs = 1:2:50;
% ys = xs;
% zs = xs;
% absorb_copy = absorb;
% % absorb_copy(absorb_copy<0.4)=NaN;
% absorb_copy(absorb_copy==0)=NaN;
% %absorb_copy(1:M/2,1:M,2/M:M)=NaN;
% absorb_copy(1:M/2,1:M,1:M)=NaN;
% h = slice(absorb_copy,xs,ys,zs);
% set(h,'FaceColor','interp','EdgeColor','none')
% camproj perspective
% box on
% view(-35,35)
% axis([1 M 1 M 1 M])
% colormap jet
% colorbar


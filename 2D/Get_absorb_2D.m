% K:parameters matrix(absorption, scattering, refractive index)
% Q:matrix of column vectors,each column consists of one source distribution
% M:matrix of column vectors for the detector response distributions
% Phi: matrix of photon density distributions for each source
% Y:measurements for each source and detector combination.

%===========================%
%                           %
%      INHOMOGENEOUS        %
%                           %
%===========================%
close all;
clear all;
clc;
rad = 25;   % mesh radius [mm]
nsect = 8;  % number of sectors 
nring = 32; % number of rings
nbnd = 4;   % number of boundary rings
[vtx,idx,eltp] = mkcircle (rad, nsect, nring, nbnd);
mesh = toastMesh(vtx,idx,eltp);

% [vtx,idx,eltp] = mesh.Data;
% nnode = mesh.NodeCount;
% nel = mesh.ElementCount;
% **************************************
% Step 1: Creating parameter maps
% **************************************
% bmua = imread('demo_matlab_fwd2_mua.png');
% bmus = imread('demo_matlab_fwd2_mus.png');
% bmua = imread('dustin_mua.png');
% bmus = imread('dustin_mus.png');

% (76,56) (59,73) 
bmua = imread('mua_2Ps.png');
bmus = imread('mua_2Ps.png');

% %(62,63)
% bmua = imread('mua_1P.png');
% bmus = imread('mus_1P.png');

% bmua = zeros(128,128);
% bmua(64,32) = 741/15;
% bmua(64,96) = 741/15; 
% bmus = ones(128,128);
% bmus(64,32) = 900/1670;
% bmus(64,96) = 900/1670; 

% bmua = double(bmua)./255.*0.02 + 0.01;
% bmus = double(bmus)./255.*1.0 + 1.0;

bmua = double(bmua)./255.*0.12;
bmus = double(bmus)./255.*1.0+1;
bmua_0 = zeros(128,128);
bmua_0 = double(bmua_0)./255.*0.12;

% ***************************************
% Step 2: Mapping the images in the mesh basis
% ***************************************
grd = size(bmua);
basis = toastBasis(mesh,grd);
mua = basis.Map('B->M',bmua);
mus = basis.Map('B->M',bmus);
%===============================================================
% Map a scalar function defined over the problem domain from
%   one basis representation to another.
%  
%   Syntax: tgt = basis.Map(mapstr,src)
%  
%   Arguments:
%           mapstr [string]:
%               specifies source and target basis representations
%               (see notes)
%           src [real or complex array]:
%               function coefficients in source basis
%  
%   Return values:
%           tgt [real or complex array]:
%               function coefficients in target basis
%  
%   Notes:  The mapstr argument has the format 'X->Y', where X
%           is a character representing the source basis, and Y a
%           character representing the target basis. The
%           following are supported for X and Y:
%  
%           X/Y : basis
%  
%           M   : mesh basis
%           B   : regular grid basis
%           S   : regular grid basis, excluding external points
%           G   : intermediate grid basis
%  
%   Example:
%           simg = basis.Map('M->S', mimg)
%           maps function 'mimg', represented by nodal values of
%           the unstructured mesh, to regular grid basis 'simg'
%           excluding grid points external to the mesh.
%===============================================================

% figure('Name','mua');
% mesh.Display;
% hold on;
% mesh.Display(mua);
% figure('Name','mus');
% mesh.Display;
% hold on;
% mesh.Display(mus);

nnd = mesh.NodeCount;
ref_bkg = 1.4;
ref = ones(nnd,1) * ref_bkg;

% *****************************************
% Step 3: Invoking the forward solver
% *****************************************
% Create the source and detector positions
nq = 16;
for i=1:nq
  phi_q = 2*pi*(i-1)/nq;
  Q(i,:) = rad * [cos(phi_q) sin(phi_q)];
  phi_m = 2*pi*(i-0.5)/nq;
  M(i,:) = rad * [cos(phi_m) sin(phi_m)];
end
mesh.SetQM(Q,M);
% figure('Name','source and detector positions');
% mesh.Display;
% hold on
% plot(Q(:,1),Q(:,2),'ro','MarkerFaceColor','r');
% plot(M(:,1),M(:,2),'bs','MarkerFaceColor','b');

clear rad nsect nring nbnd
% Create the source and boundary projection vectors
qvec = mesh.Qvec ('Neumann', 'Gaussian', 2);
mvec = mesh.Mvec ('Gaussian', 2);

% Solve the FEM linear system
K = dotSysmat (mesh,mua,mus,ref,0);
Phi = K\qvec;  % number of node * number of source
%Y = mvec.' * Phi;

% % Display sinogram
% figure('Name','display the measurements as a sinogram');
% imagesc(log(Y));
% xlabel('source index q');
% ylabel('detector index m');
% axis equal tight;
% colorbar

% % Display boundary profile
% figure('Name','display the measurements as a boundary profile');
% hold on
% angle = [360/32:360/16:360];
% for i=1:size(Y,2)
%     ywrap = [Y(i:end,i); Y(1:i-1,i)];
%     plot(angle,log(ywrap),'o-');
% end
% axis([0 360 -14 -2]);
% xlabel('angular source-detector separation');
% ylabel('log intensity');
% 
% % Write solver results to file
% data = reshape(log(Y'),[],1);
% toastWriteVector('demo_matlab_fwd2.dat', data);

%Show Phi
% Phi_1 = basis.Map('M->B',Phi_1);
% Phi_1=reshape(Phi_1,128,128);
PPP_mua=zeros(128,128,16);
for i=1:16
    PPP_mua(:,:,i) = reshape(basis.Map('M->B',Phi(:,i)),128,128);
end
% figure('Name','Phi_1');
% imagesc(PPP_mua(:,:,1));

%Generate optical density when Phi=0(Phi_0)
mua_0 = basis.Map('B->M',bmua_0);
% Solve the FEM linear system
K_0 = dotSysmat (mesh,mua_0,mus,ref,0);
Phi_0 = K_0\qvec;
%Y_0 = mvec.' * Phi_0;

PPP_mua_0=zeros(128,128,16);
for i=1:16
    PPP_mua_0(:,:,i) = reshape(basis.Map('M->B',Phi_0(:,i)),128,128);
end
% figure('Name','Phi_0_1');
% imagesc(PPP_mua_0(:,:,1));
absorb = zeros(128,128,16);
for i=1:16
    absorb(:,:,i) = PPP_mua_0(:,:,i) - PPP_mua(:,:,i);
end
absorb = sum(absorb,3);
clearvars -except absorb

figure('Name','absorb');
imagesc(absorb);
function Display3D(target,mode)
% DESCRIPTION:
%       % Visualisation of 3D matrix
%
% INPUTS:
%       target          - target
%       Nx,Ny,Nz        - size of target
%       viewA,viewB     - view direction
%       mode            - display mode 
% ABOUT:
%       author          - Dustin Qi
%       date            - 8th Dec 2017

% Copyright (C) 2015-2018 Dustin Qi
    target_d=double(target);
    tmp=size(target_d);
    Nx=tmp(1,1);
    Ny=tmp(1,2);
    Nz=tmp(1,3);
    xs = 1:2:Nx;
    ys = 1:2:Ny;
    zs = 1:2:Nz;
    absorb_copy = target_d;
    
    if mode==1
        absorb_copy(1:Nx,1:Ny/2,1:Nz)=NaN;
%         viewA=-90;
%         viewB=0;
        viewA=-115;
        viewB=23;        
    elseif mode==2
        absorb_copy(1:Nx,1:(Ny/2-1),1:Nz)=NaN;
        viewA=-115;
        viewB=23;  
    end

%     absorb_copy(1:Nx,1:Ny,Nz/part:Nz)=NaN;
    h = slice(absorb_copy,xs,ys,zs);
    set(h,'FaceColor','interp','EdgeColor','none')
    camproj perspective
    box on
    view(viewA,viewB)
    axis([1 Nx 1 Ny 1 Nz])
    colormap jet
    colorbar
function index2cart3D(index,imageSize)
% DESCRIPTION:
%       % display positions of every layer's sensors 
%
% INPUTS:
%       index          - index information, array of interger, size:[M,1] 
%       size           - size of 3Dimage, like[30,30,30]
% ABOUT:
%       author          - Dustin Qi
%       date            - 28th Dec 2017
% Copyright (C) 2015-2018 Dustin Qi

Nx=imageSize(1,1);
Ny=imageSize(1,2);
Nz=imageSize(1,3);
len = size(index,1);
layer = zeros(1,len);
for i=1:len
    tmp = floor(index(i,1)/(Nx*Ny));
    layer(1,i) = tmp;
    clear tmp
end
layer = unique(layer);

for j = 1:size(layer,2)
    tmp = zeros(len,1);
    for i=1:len
        if layer(1,j) == floor(index(i,1)/(Nx*Ny))
            tmp(i,1)=index(i,1)-layer(1,j)*Nx*Ny;
            tmp(find(tmp==0))=[];
            figure(j);
            index2cart2D(tmp,Nx,Ny);
        end
    end
    clear tmp
end


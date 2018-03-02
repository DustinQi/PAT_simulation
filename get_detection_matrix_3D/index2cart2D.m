%根据探测器位置信息画出探测器真实位置
function index2cart(index,size1,size2)
% DESCRIPTION:
%       % index information to cart information and display
%
% INPUTS:
%       index          - index information, array of interger, size:[M,1] 
%       size           - size of image
% ABOUT:
%       author          - Dustin Qi
%       date            - 20th Dec 2017
% Copyright (C) 2015-2018 Dustin Qi

cartinfo=zeros(size1,size2);
index = double(index);
for i=1:size(index,1)
    carty = floor(index(i,1)/size1);
    cartx = index(i,1) - carty*size1;
    cartinfo(cartx,carty+1) = 1; 
end
imagesc(cartinfo);
colormap(hot);
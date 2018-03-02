%% ART based on FLA
clear all
close all
clc;
L = 128;%图像大小
nIter = 0;% 迭代次数
iterOK = 0;
iter0 = ones(4*L)*0.2;%元素为0.2的矩阵
theta0 = [0:14; 15:29; 30:44; 45:59; 60:74; 75:89; 90:104; 105:119; 120:134; 135:149; 150:164; 165:179];
theta = reshape (theta0, 180, 1);%theta = [0,15,30,...165,1,16,31,...166,2,17,32,...167,......]'
nTheta = length(theta);% 180

%% 生成投影，以及初始迭代矩阵
% I = phantom (L);% 原始图像
I = p0;

R1 = zeros(128,128);% 迭代初值为0

% P1 = radon (I, 0:30:179);
% R1 = iradon (P1, 0:30:179, 'linear', 'hamming');% 迭代初始值为R1(效果好)

for aa = 1:L %%
    for bb = 1:L
        iter0 ((4*aa-3):4*aa, (4*bb-3):4*bb) = R1 (aa, bb)/16;
        I1 ((4*aa-3):4*aa, (4*bb-3):4*bb) = I (aa, bb)/16;% 把原图扩大4倍
    end
end

%% 开始迭代  每考虑一个投影角度更新一次图像
while (nIter <= 10)% 迭代次数
    nIter = nIter + 1;
    temp = mod(nIter-1, nTheta)+1;
%一次迭代过程
    for ii = 1:L
        %确定射线束的点
        proj = zeros (4*L);
        proj (:,(4*ii-3):4*ii) = 1;% (4*ii-3)~4*ii 列的值为1
        proj1 = imrotate (proj, theta (temp), 'crop');% 确定射线束
        P_proj = sum(sum(I1.*proj1));
        R_proj = sum(sum(iter0.*proj1));
        pp_proj = sum(sum(proj1 .* proj1));
        iter1 = iter0 - proj1.*(R_proj-P_proj)/pp_proj;% sum(sum(proj))
        iter0 = iter1;    
    end
%一次迭代过程结束
end
% imshow(iter0);

%% 重建后图像
R = zeros(L);
for aa = 1:L;
    for bb = 1:L;
        R (aa, bb) = sum(sum(iter0((4*aa-3):4*aa,(4*bb-3):4*bb)));
    end
end
figure, imshow(R), title('迭代重建图像')

function [W,b,sPos] = randSensorRecon(p0,allPos,sNum,timeStep)
% DESCRIPTION:
%       % get W and b with random sensors
%
% INPUTS:
%       p0              - initial pressure
%       allPos          - all potential positions of sensors
%       sNum            - number of sensors used to recon
%       timeStep        - number of time steps
% OUTPUTS:
%       W               - detection matrix
%       b               - sensor data
%       sPos            - positions of sensors used to recon
% ABOUT:
%       author          - Dustin Qi
%       date            - 26th Dec 2017

% Copyright (C) 2015-2018 Dustin Qi

c = 1500; %sound speed
dt = 3.2e-8; %s
dx = 0.25e-3; %m

% % 得到采样变换矩阵
absize=size(p0);
Nx = absize(1,1);
Ny = absize(1,2);
Nz = absize(1,3);
T=timeStep; % 采样点数

det_pos=double(allPos); %探测器信息
det_num=size(det_pos,1); % 探测器个数
    
    disp('生成随机探测器');
    %选取随机个数探测器
    rdet = zeros(1,sNum);
    for i=1:sNum
        tmp = unidrnd(det_num);
        while ismember(tmp,rdet)
            tmp = unidrnd(det_num);
        end
        rdet(1,i) = tmp;
        clear tmp i
    end    
    
    disp('开始构建权重矩阵');
    disp('@_@ ');
    %单探测器测量矩阵
    for di=1:size(rdet,2)
        det_index = rdet(1,di);
        c_z=floor(det_pos(det_index,1)/(Nx*Ny));
        cxy=det_pos(det_index,1)-(Nx*Ny)*c_z;
        c_y=floor(cxy/Nx);
        c_x=cxy-Nx*c_y;
        center_p=[c_x,c_y+1,c_z+1]; %探测器位置 
        clear c_x c_y c_z cxy   
        mask_T=zeros(Nx,Ny,Nz,T);

        %单个探测点的测量矩阵
        for t=1:T         
            mask=zeros(Nx,Ny,Nz);  
            for i=1:Nx
                for j=1:Ny
                    for k=1:Nz
                        d = sqrt((i-center_p(1,1))^2+(j-center_p(1,2))^2+(k-center_p(1,3))^2);     
                        ab=abs((d*dx)/(c*dt)-t);
                        if ab<1
                            mask(i,j,k)=1-ab;
                        end 
                    end
                end
            end
            mask_T(:,:,:,t)=mask;
            clear mask
         end

        eval([ 'A',num2str(di), '=', 'mask_T', ';']);
        clear mask_T 
        display([num2str(di),'/',num2str(size(rdet,2)),' 已完成...']);
    end    
    clear center_p d di k t ab i j

    %将累加权重归一化为一列
    for i=1:size(rdet,2)
        eval([ 'A' '=', 'A',num2str(i), ';']);
        tmp=zeros(Nx*Ny*Nz,T);
        for j = 1:T
            tmp(:,j)=reshape(A(:,:,:,j),[Nx*Ny*Nz,1]); % size:(Nx*Ny*Nz,T)
        end
        eval([ 'A' ,num2str(i),'=', 'tmp',';']);
        clear tmp
    end
    clear A i j
    disp('构建权重矩阵完成！');
        
    %实际探测器数据g
    p = reshape(p0,[Nx*Ny*Nz,1]);
    for i=1:size(rdet,2)
        g=zeros(Nx*Ny*Nz,T); %size(16384,300)
        eval([ 'A', '=', 'A',num2str(i),';']);
        for j=1:T
            g(:,j)=A(:,j).*p;  
        end

        g=sum(g);
        g=g';
        eval([ 'g',num2str(i) '=', 'g',';']);
        clear g A
    end
    clear i j  

    %合并所有测量矩阵
    M=A1';
    for i=1:(size(rdet,2)-1)
        eval([ 'L', '=', 'A',num2str(i+1),';']);
        L=L';
        M=[M;L]; 
    end
    clear L i 

    %合并所有探测器数据
    g=g1;
    for i=1:(size(rdet,2)-1)
        eval([ 'h', '=', 'g',num2str(i+1),';']);
        g=[g;h]; 
    end
    clear h i

    W = M;
    b = g;
    sPos = det_pos(rdet);

% %尝试lsqr进行重建
% X=lsqr(W,b,1e-3);
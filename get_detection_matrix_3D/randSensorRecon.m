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

% % �õ������任����
absize=size(p0);
Nx = absize(1,1);
Ny = absize(1,2);
Nz = absize(1,3);
T=timeStep; % ��������

det_pos=double(allPos); %̽������Ϣ
det_num=size(det_pos,1); % ̽��������
    
    disp('�������̽����');
    %ѡȡ�������̽����
    rdet = zeros(1,sNum);
    for i=1:sNum
        tmp = unidrnd(det_num);
        while ismember(tmp,rdet)
            tmp = unidrnd(det_num);
        end
        rdet(1,i) = tmp;
        clear tmp i
    end    
    
    disp('��ʼ����Ȩ�ؾ���');
    disp('@_@ ');
    %��̽������������
    for di=1:size(rdet,2)
        det_index = rdet(1,di);
        c_z=floor(det_pos(det_index,1)/(Nx*Ny));
        cxy=det_pos(det_index,1)-(Nx*Ny)*c_z;
        c_y=floor(cxy/Nx);
        c_x=cxy-Nx*c_y;
        center_p=[c_x,c_y+1,c_z+1]; %̽����λ�� 
        clear c_x c_y c_z cxy   
        mask_T=zeros(Nx,Ny,Nz,T);

        %����̽���Ĳ�������
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
        display([num2str(di),'/',num2str(size(rdet,2)),' �����...']);
    end    
    clear center_p d di k t ab i j

    %���ۼ�Ȩ�ع�һ��Ϊһ��
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
    disp('����Ȩ�ؾ�����ɣ�');
        
    %ʵ��̽��������g
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

    %�ϲ����в�������
    M=A1';
    for i=1:(size(rdet,2)-1)
        eval([ 'L', '=', 'A',num2str(i+1),';']);
        L=L';
        M=[M;L]; 
    end
    clear L i 

    %�ϲ�����̽��������
    g=g1;
    for i=1:(size(rdet,2)-1)
        eval([ 'h', '=', 'g',num2str(i+1),';']);
        g=[g;h]; 
    end
    clear h i

    W = M;
    b = g;
    sPos = det_pos(rdet);

% %����lsqr�����ؽ�
% X=lsqr(W,b,1e-3);
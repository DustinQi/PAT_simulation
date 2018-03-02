%% 40^3=64000, PML = 5, (40+5*2)^3=125000
%������άreshape��˳�򣺲ο��������ӣ�
% 
% >> a=rand(2,2,2)
% a(:,:,1) =
%     0.6787    0.7431
%     0.7577    0.3922
% a(:,:,2) =
%     0.6555    0.7060
%     0.1712    0.0318
% >> reshape(a,1,8)
% ans =
%     0.6787    0.7577    0.7431    0.3922    0.6555    0.1712    0.7060    0.0318

%% p0 = resize(absorb, [Nx,Ny,Nz]);
c = 1500; %sound speed
dt = 5e-8; %s
% dx = 0.16e-3; %m
dx = 0.25e-3;

% % �õ������任����
p0 = absorb;
absize=size(p0);
Nx = absize(1,1);
Ny = absize(1,2);
Nz = absize(1,3); 
T=300; % ��������

det_pos=double(sensor_norm8); %̽������Ϣ
det_num=size(det_pos,1); % ̽��������

%% �õ��ۼ�Ȩ��det_num��A      
for di=1:det_num
    di
    
    c_z=floor(det_pos(di,1)/(Nx*Ny));
    cxy=det_pos(di,1)-(Nx*Ny)*c_z;
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
    disp('Done');
end    
clear center_p d di k t ab i j

%���ۼ�Ȩ�ع�һ��Ϊһ��
for i=1:det_num
    eval([ 'A' '=', 'A',num2str(i), ';']);
    tmp=zeros(Nx*Ny*Nz,T);
    for j = 1:T
        tmp(:,j)=reshape(A(:,:,:,j),[Nx*Ny*Nz,1]); % size:(Nx*Ny*Nz,T)
    end
    eval([ 'A' ,num2str(i),'=', 'tmp',';']);
    clear tmp
end
clear A i j

%% ʵ��̽��������g
p = reshape(p0,[Nx*Ny*Nz,1]);
for i=1:det_num
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

%% �ϲ����в�������
M=A1';
for i=1:(det_num-1)
    eval([ 'L', '=', 'A',num2str(i+1),';']);
    L=L';
    M=[M;L]; 
end
clear L i 

%% �ϲ�����̽��������
g=g1;
for i=1:(det_num-1)
    eval([ 'h', '=', 'g',num2str(i+1),';']);
    g=[g;h]; 
end
clear h i


W = M;
b = g;
% %����l1_ls�����ؽ�
% lambda = 0.01;
% rel_tol = 0.01;
% tic;
% [X_ls,status] = l1_ls(W,b,lambda,rel_tol);
% [X_ls,status] = l1_ls_nonneg(W,b,lambda,rel_tol);
% toc;

% %����lsqr�����ؽ�
% disp('����lsqr�����ؽ�...');
% tic;
% X=lsqr(W,b,1e-3);
% toc;
% 
% xx=reshape(X,[30,30,30]);
% Display3D(xx,1);
% disp('ȫ������...');

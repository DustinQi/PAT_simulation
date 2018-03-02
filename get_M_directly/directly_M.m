% % from ����ϡ�����ԵĹ���ͼ�񷽷��о� (��ȷ�汾)

c = 1500; %sound speed
dt = 2e-8; %s
dx = 1e-4; %m

% �õ������任����
N=128;
T=250; % ��������
mask_T=zeros(N,N,T);

record_info=record_8_R58; %̽������Ϣ
det_num=size(record_info.tri,1); % ̽��������
det_pos=zeros(det_num,1); % ��ʼ��̽������λ��

%ȷ��̽����λ��
for i = 1:det_num 
    tmp=find(max(record_info.bc(i,:))==record_info.bc(i,:));
    det_pos(i,:)=record_info.tri(i,tmp);
    clear tmp i
end

%�õ��ۼ�Ȩ��det_num��A      
for di=1:det_num
    di
    
    c_col=floor(det_pos(di,1)/N);
    c_row=det_pos(di,1)-128*c_col;
    center_p=[c_row,c_col+1]; %Բ��λ�� 
    clear c_row c_col
%     mask=zeros(N,N);       
    %����̽���Ĳ�������
    for t=1:T     
        mask=zeros(N,N);  
        for i=1:N
            for j=1:N
                d = sqrt((i-center_p(1,1))^2+(j-center_p(1,2))^2);     
                ab=abs((d*dx)/(c*dt)-t);
                if ab<1
                    mask(i,j)=1-ab;
                end 
            end
            clear i j l Rij
        end
        mask_T(:,:,t)=mask;
    end
    
    eval([ 'A',num2str(di), '=', 'mask_T', ';']);
    clear mask_T mask
    disp('Done');
end    
clear center_p d di k t ab

%���ۼ�Ȩ�ع�һ��Ϊһ��
for i=1:det_num
    eval([ 'A' '=', 'A',num2str(i), ';']);
    tmp=zeros(N*N,T);
    for j = 1:T
        tmp(:,j)=reshape(A(:,:,j),N*N,1); % size:(16384,T)
    end
    eval([ 'A' ,num2str(i),'=', 'tmp',';']);
    clear tmp
end
clear A i j


%ʵ��̽��������g
p = reshape(p1,16384,1);
for i=1:det_num
    g=zeros(N*N,T); %size(16384,300)
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
for i=1:(det_num-1)
    eval([ 'L', '=', 'A',num2str(i+1),';']);
    L=L';
    M=[M;L]; 
end
clear L i

%�ϲ�����̽��������
g=g1;
for i=1:(det_num-1)
    eval([ 'h', '=', 'g',num2str(i+1),';']);
    g=[g;h];
end
clear h i

% ===========ʹ��l1_ls�����ؽ�==========[X,status] = l1_ls(A,y,lambda,rel_tol,quiet)======
W=M;
b=g;
lambda = 0.01;
rel_tol = 0.01;
tic;
[X_ls,status] = l1_ls(W,b,lambda,rel_tol);
toc;

X_ls = reshape(X_ls,128,128);
imagesc(X_ls);
colormap(jet);

% % ===========ʹ��lsqr�����ؽ�==========[X,status] = l1_ls(A,y,lambda,rel_tol,quiet)======
% W=M;
% b=g;
% X_lsqr = lsqr(W,b);
% 
% tic;
% X_lsqr = reshape(X_lsqr,128,128);
% toc;
% imagesc(X_lsqr);
% colormap(jet);
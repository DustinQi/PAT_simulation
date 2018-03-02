%=====================
% b=W*x
% sesor_data=A9*p0;
%=====================

% %拼接sensor信号
% f=zeros(2400,1);
% for i=1:8
%     f(1+300*(i-1):300*i,1)=f_all(i,:);
% end
% clear i

% %拼接测量矩阵
% M=[A1;A2;A3;A4;A5;A6;A7;A8;A9;A10;A11;A12;A13;A14;A15;A16];



% % ====================使用art进行重建========================
% clc;
% % b = g;
% X0 = zeros(16384,1)';
% W = M;
% Niter = 100;
% e=0.01;
% [X_art1]=art1(W,b,X0,Niter,e);
% 
% X_art1=reshape(X_art1,128,128);
% imagesc(X_art1);
% colormap(jet);

% b = g1;
% W = A1;
% Niter = 100;
% [X_art2,rho,eta]=art(W,b,Niter);
% X_art2=reshape(X_art2,128,128);
% imagesc(X_art2);
% colormap(jet);

% % ====================使用yall1进行重建=====================
% W=A1;
% b=g1;
% opts.rho=5e-3;
% opts.tol=5e-3;
% opts.nonneg=1;
% X_yall1=yall1(W,b,opts);
% X_yall1=reshape(X_yall1,128,128);
% imagesc(X_yall1);
% colormap(jet);

% ===========使用l1_ls进行重建==========[X,status] = l1_ls(A,y,lambda,rel_tol,quiet)======
W = M;
b = g;
lambda = 0.01;
rel_tol = 0.01;

tic;
[X_ls,status] = l1_ls(W,b,lambda,rel_tol);
% [X_ls,status] = l1_ls_nonneg(W,b,lambda,rel_tol);
toc;

X_ls = reshape(X_ls,[128,128]);
imagesc(X_ls);
colormap(jet);

%展示重建图像和原始图像
% figure(1);
% imagesc(X);
% figure(2);
% imagesc(p0);
% subplot(1,2,1);imagesc(X);title('ART重建');
% subplot(1,2,2);imagesc(p0);title('原始图像');

% %画X和p0得斜对角线值作比较
% diag_X=zeros(1,128);
% diag_p0=zeros(1,128);
% X=X_art;
% p0=X_art_pure;
% for i=1:128
%     diag_X(1,i)=X(129-i,i);
%     diag_p0(1,i)=p0(129-i,i);
% end
% figure(2);
% plot(diag_X,'r*');
% axis([0 128 0 1.2]);
% hold on;
% plot(diag_p0,'bo');
% axis([0 128 0 1.2]);
% clear diag_X diag_p0
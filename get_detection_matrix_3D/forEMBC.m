%% lsqr�����ؽ�
% disp('==========lsqr�ؽ�==========');
% tic;
% X_lsqr=lsqr(W,b,1e-6);
% toc;
% xx_lsqr=reshape(X_lsqr,[30,30,30]);
% Display3D(xx_lsqr,1);
% disp('lsqr�ؽ�����');

%% yall1�����ؽ� 
% clc;
% disp('==========yall1�ؽ�==========');
% opts.rho=1e-1; 
% opts.tol=5e-3;
% % opts.weights = ws2;
% opts.maxit = 5000;
% tic;
% X_yall1=yall1(W,b,opts);
% toc;
% xx_yall1=reshape(X_yall1,[30,30,30]);
% Display3D(xx_yall1,1);
% disp('yall1�ؽ�����');

%% ISD�����ؽ�
% opts.rho=1e-1;
% opts.tol=5e-3;
% opts.maxit = 1;
% % opts.weights = x_real;
% opts.suppDetct = 0;
% tic;
% [X_isd,Out] = Threshold_ISD_1D(W,b,opts);
% toc;
% xx_isd=reshape(X_isd,[30,30,30]);  
% Display3D(xx_isd,1);

%% yall_group
% beta = [1,1]/mean(abs(b));
% tol=5e-4;
% 
% ppp1=ppp(1:13500,1);
% ppp2=ppp(13501:27000,1);
% ppp1(find(ppp1<8.5e-3))=3;
% ppp1(find(ppp1~=3))=1;
% ppp2(find(ppp2<8.5e-3))=3;
% ppp2(find(ppp2~=3))=2;
% groups=[ppp1;ppp2];
% 
% weights = [1,1,100]';
% [x_yp,Out_yp] = YALL1_group(W,b,groups,...
%                                'StopTolerance', tol, ...
%                                'GrpWeights', weights, ...
%                                'overlap', false, ...
%                                'nonorthA', true, ...
%                                'ExactLinSolve', true, ...
%                                'Solver', 1, ...             % 1 - primal; 2 - dual
%                                'QuadPenaltyPar', beta, ...
%                                'maxIter', 500);
% xx_yp=reshape(x_yp,[30,30,30]);
% Display3D(xx_yp,1);

% slice = 8;
% s_lsqr=xx_lsqr(:,:,slice);
% subplot(1,3,1);
% imagesc(s_lsqr);
% 
% s_yall1=xx_yall1(:,:,slice);
% subplot(1,3,2);
% imagesc(s_yall1);
% 
% s_w_yall1=xx_w_yall1(:,:,slice);
% subplot(1,3,3);
% imagesc(s_w_yall1);






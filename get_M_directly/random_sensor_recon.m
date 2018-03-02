% % from 基于稀疏特性的光声图像方法研究 (正确版本)

c = 1500; %sound speed
dt = 2e-8; %s
dx = 1e-4; %m

% 得到测样变换矩阵
N=128;
T=370; % 采样点数
mask_T=zeros(N,N,T);

record_info=record_32_R58; %探测器信息
det_num=size(record_info.tri,1); % 探测器总个数

randomRconTime = 4; %随机重建的次数
rdet_num = 8; %选取随机n个探测器用于重建

for rrect=1:randomRconTime
    display(['开始第 ',num2str(rrect),' 次随机重建,共',num2str(randomRconTime),'次请蛋定。']);
    
    %选取随机个数探测器
    rdet = zeros(1,rdet_num);
    for i=1:rdet_num
        tmp = unidrnd(32);
        while ismember(tmp,rdet)
            tmp = unidrnd(32);
        end
        rdet(1,i) = tmp;
        clear tmp i
    end
    
    %确定探测器位置
    det_pos=zeros(rdet_num,1); % 初始化探测器的位置
    for i = 1:size(rdet,2)
        det_index = rdet(1,i);
        tmp=find(max(record_info.bc(det_index,:))==record_info.bc(det_index,:));
        det_pos(i,:)=record_info.tri(det_index,tmp);
        clear tmp i det_index
    end
    
    disp('开始构建权重矩阵，请等待......');
    %得到累加权重det_num个A      
    for di=1:size(rdet,2)
        c_col=floor(det_pos(di,1)/N);
        c_row=det_pos(di,1)-128*c_col;
        center_p=[c_row,c_col+1]; %圆心位置 
        clear c_row c_col
    %     mask=zeros(N,N);       
        %单个探测点的测量矩阵
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
    end    
    clear center_p d di k t ab
    disp('构建权重矩阵完成！');
    disp('关于权重矩阵可参考《基于全变分法重建光声图像 张砚》或参考本段代码。');

    %将累加权重归一化为一列
    for i=1:size(rdet,2)
        eval([ 'A' '=', 'A',num2str(i), ';']);
        tmp=zeros(N*N,T);
        for j = 1:T
            tmp(:,j)=reshape(A(:,:,j),N*N,1); % size:(16384,T)
        end
        eval([ 'A' ,num2str(i),'=', 'tmp',';']);
        clear tmp
    end
    clear A i j 

    %实际探测器数据g
    p = reshape(p1,16384,1);
    for i=1:size(rdet,2)
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

%     % ===========使用l1_ls进行重建==========[X,status] = l1_ls(A,y,lambda,rel_tol,quiet)======
%     W=M;
%     b=g;
%     lambda = 0.01;
%     rel_tol = 0.01;
%     [X,status] = l1_ls(W,b,lambda,rel_tol);
    
    % ===========使用lsqr进行重建==========[X,status] = l1_ls(A,y,lambda,rel_tol,quiet)======
    W=M;
    b=g;
    X = lsqr(W,b);    
    
    eval([ 'rdet_pos_',num2str(rrect) '=', 'det_pos',';']);
    eval([ 'X',num2str(rrect) '=', 'X',';']);
    clear X
end
clear rrect status p det_pos A1 A2 A3 A4 A5 A6 A7 A8 g g1 g2 g3 g4 g5 g6 g7 g8 
clear lambda randomReconTime rdet rdet_num record_info W b rel_tol X_ls
disp('随机重建全部完成');

% %随机选取8个探测器
% rdet_num = 8;
% rdet = zeros(1,rdet_num);
% for i=1:rdet_num
%     tmp = unidrnd(32);
%     while ismember(tmp,rdet)
%         tmp = unidrnd(32);
%     end
%     rdet(1,i) = tmp;
%     clear tmp
% end

% %确定探测器位置
% det_pos=zeros(rdet_num,1); % 初始化探测器的位置
% for i = 1:size(rdet,2)
%     det_index = rdet(1,i);
%     tmp=find(max(record_info.bc(det_index,:))==record_info.bc(det_index,:));
%     det_pos(i,:)=record_info.tri(det_index,tmp);
%     clear tmp i det_index
% end

% %得到累加权重det_num个A      
% for di=1:det_num
%     c_col=floor(det_pos(di,1)/N);
%     c_row=det_pos(di,1)-128*c_col;
%     center_p=[c_row,c_col+1]; %圆心位置 
%     clear c_row c_col
% %     mask=zeros(N,N);       
%     %单个探测点的测量矩阵
%     for t=1:T     
%         mask=zeros(N,N);  
%         for i=1:N
%             for j=1:N
%                 d = sqrt((i-center_p(1,1))^2+(j-center_p(1,2))^2);     
%                 ab=abs((d*dx)/(c*dt)-t);
%                 if ab<1
%                     mask(i,j)=1-ab;
%                 end 
%             end
%             clear i j l Rij
%         end
%         mask_T(:,:,t)=mask;
%     end
%     
%     eval([ 'A',num2str(di), '=', 'mask_T', ';']);
%     clear mask_T mask
%     disp('Done');
% end    
% clear center_p d di k t ab

% %将累加权重归一化为一列
% for i=1:det_num
%     eval([ 'A' '=', 'A',num2str(i), ';']);
%     tmp=zeros(N*N,T);
%     for j = 1:T
%         tmp(:,j)=reshape(A(:,:,j),N*N,1); % size:(16384,T)
%     end
%     eval([ 'A' ,num2str(i),'=', 'tmp',';']);
%     clear tmp
% end
% clear A i j
% 
% %实际探测器数据g
% p = reshape(p1,16384,1);
% for i=1:det_num
%     g=zeros(N*N,T); %size(16384,300)
%     eval([ 'A', '=', 'A',num2str(i),';']);
%     for j=1:T
%         g(:,j)=A(:,j).*p;  
%     end
%     
%     g=sum(g);
%     g=g';
%     eval([ 'g',num2str(i) '=', 'g',';']);
%     clear g A
% end
% clear i j 

% %合并所有测量矩阵
% M=A1';
% for i=1:(det_num-1)
%     eval([ 'L', '=', 'A',num2str(i+1),';']);
%     L=L';
%     M=[M;L]; 
% end
% clear L i
% 
% %合并所有探测器数据
% g=g1;
% for i=1:(det_num-1)
%     eval([ 'h', '=', 'g',num2str(i+1),';']);
%     g=[g;h];
% end
% clear h i
% 
% % % ===========使用l1_ls进行重建==========[X,status] = l1_ls(A,y,lambda,rel_tol,quiet)======
% W=M;
% b=g;
% lambda = 0.01;
% rel_tol = 0.01;
% [X_ls,status] = l1_ls(W,b,lambda,rel_tol);
% eval()
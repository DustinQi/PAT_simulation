% �����ֵ�����PSNR�;������MSE

%����ģ��
template_v = zeros(30,30);
template_v(20:25,6:11)=255;
template_v(5:10,19:24)=255;

template_h_up=zeros(30,30);
template_h_up(7:11,14:17) = 255;

template_h_down = zeros(30,30);
template_h_down(20:24,14:17)=255;

%ͼ���һ����0~255
inImg=horizon_wyall1_up;
Max=255;
Min=0;
pmin=min(min(inImg));
pmax=max(max(inImg));
outImg = round((Max-Min)*(inImg-pmin)/(pmax-pmin)+Min);

%�����ֵ����Ⱥ;��������
X=outImg;
Y=template_h_up;
D=X-Y;
MSE=sqrt(sum(D(:).*D(:))/numel(X));
PSNR=10*log10(255^2/MSE);
display(MSE);
display(PSNR);
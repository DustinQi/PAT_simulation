% 得到第i组实验的探测器卦限分布情况
% position(1)~position(8)分别对应8个卦限内探测器数目
load('C:\Users\T804\Desktop\光声仿真简单结果\暑假组会报告\第三次\随机6探测器\对应的6探测器位置.mat');
i=4;
pp = out_sensor_mask(:,:,i);
l=size(pp);
l=l(2);
position=[0,0,0,0,0,0,0,0];
for j=1:l
    temp=pp(:,j);
    temp=temp';
    if temp(1)>0 && temp(2)>0 && temp(3)>0
        position=position+[1,0,0,0,0,0,0,0];
    elseif temp(1)<0 && temp(2)>0 && temp(3)>0
        position=position+[0,1,0,0,0,0,0,0];
    elseif temp(1)<0 && temp(2)<0 && temp(3)>0
        position=position+[0,0,1,0,0,0,0,0];
    elseif temp(1)>0 && temp(2)<0 && temp(3)>0
        position=position+[0,0,0,1,0,0,0,0];
    elseif temp(1)>0 && temp(2)>0 && temp(3)<0
        position=position+[0,0,0,0,1,0,0,0];
    elseif temp(1)<0 && temp(2)>0 && temp(3)<0
        position=position+[0,0,0,0,0,1,0,0];
    elseif temp(1)<0 && temp(2)<0 && temp(3)<0
        position=position+[0,0,0,0,0,0,1,0];
    elseif temp(1)>0 && temp(2)<0 && temp(3)<0
        position=position+[0,0,0,0,0,0,0,1];
    end
end
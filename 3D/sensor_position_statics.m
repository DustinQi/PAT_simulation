% �õ���i��ʵ���̽�������޷ֲ����
% position(1)~position(8)�ֱ��Ӧ8��������̽������Ŀ
load('C:\Users\T804\Desktop\��������򵥽��\�����ᱨ��\������\���6̽����\��Ӧ��6̽����λ��.mat');
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
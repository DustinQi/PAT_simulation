clc;
p = p0;
allPos = sensor64;
sNum = 8;
timeStep = 250;
reconTime=8;

for i=1:reconTime
    
    [W,g,pos]=randSensorRecon(p,allPos,sNum,timeStep);
    X=lsqr(W,g);
    eval([ 'X' ,num2str(i),'=', 'X',';']);
    eval([ 'pos' ,num2str(i),'=', 'pos',';']);
    clear X pos W g    
end
clear i
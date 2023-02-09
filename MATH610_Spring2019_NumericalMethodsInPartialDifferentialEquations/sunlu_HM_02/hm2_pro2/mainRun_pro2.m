%numIterations = 7; 
%baseNumElements = 10; %num elements on first solution
%numElements = baseNumElements*2.^(0:(numIterations-1));

numIterations = 8;
%numElements = [70,80,90,100,110,120,130,140];
numElements = 1:1:8;



epsillion = 1;
InfError = zeros(numIterations,3);
for i = 1:1:3
    fprintf('method %1d\n',i);
    [mesh,Inf_error] = mainErrorAnalysis(i,epsillion,numIterations,numElements);
    InfError(:,i)= Inf_error;
    
end

Y = 0:0.1:1;
X = epsillion * ones(size(Y));

i = 1:1:numIterations;
lmbd=figure;
loglog(mesh(i),InfError(i,1),'-+',mesh(i),InfError(i,2),'-+',mesh(i),InfError(i,3));%,'-+',X,Y,':');
%plot(mesh(i),InfError(i,1),'-+',mesh(i),InfError(i,2),'-+',mesh(i),InfError(i,3),'-+',X,Y,':');
xlabel('mesh size');
ylabel('L^{Inf} error');
legend('unpwind','downwind','centered');%,'|hb|=1');
title({['homework ',num2str(2)],['epsillion= ',num2str(epsillion),',1-8']});
saveas(lmbd, ['\\blender\homes\l\u\lusun8825\nt\AccountSettings\Desktop\610\FEM1D_For_Students\sunlu\result\',...
    'EX2_hm2_2','_epmax.png'],'png');



epsillion = 10^(-2);
InfError = zeros(numIterations,3);
for i = 1:1:3
    fprintf('method %1d\n',i);
    [mesh,Inf_error] = mainErrorAnalysis(i,epsillion,numIterations,numElements);
    InfError(:,i)= Inf_error;
    
end

Y = 0:0.1:1;
X = epsillion * ones(size(Y));

i = 1:1:numIterations;
lmbd=figure;
loglog(mesh(i),InfError(i,1),'-+',mesh(i),InfError(i,2),'-+',mesh(i),InfError(i,3),'-+');%,X,Y,':');

%plot(mesh(i),InfError(i,1),'-+',mesh(i),InfError(i,2),'-+',mesh(i),InfError(i,3),'-+',X,Y,':');
xlabel('mesh size');
ylabel('L^{Inf} error');
legend('unpwind','downwind','centered');%,'|hb|=1');
title({['homework ',num2str(2)],['epsillion= ',num2str(epsillion),',1-8']});
saveas(lmbd, ['\\blender\homes\l\u\lusun8825\nt\AccountSettings\Desktop\610\FEM1D_For_Students\sunlu\result\',...
    'EX2_hm2_2','_epmin.png'],'png');
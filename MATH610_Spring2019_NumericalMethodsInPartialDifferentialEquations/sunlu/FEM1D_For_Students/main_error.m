choose=[1,100;1,1000;1,10000;2,100;2,1000;2,10000;3,100;3,1000;3,10000;4,100;4,1000;4,10000];
for i=1:1:12
   fprintf('homework%d,S=%d\n',choose(i,1),choose(i,2));
   %num_elements=25;
   %for j=1:1:4
   %    mainFEM1D_Dirichlet(num_elements,choose(i,:));
   %    num_elements=2*num_elements;
   %end
   mainErrorAnalysis(choose(i,:));
end

choose=[0,false];
fprintf('homework%d\n',choose(1));
%num_elements=25;
%for j=1:1:4
%    mainFEM1D_Dirichlet(num_elements,choose);
%    num_elements=2*num_elements;
%end
mainErrorAnalysis(choose);





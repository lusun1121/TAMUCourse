choose=[1,3;2,3;2,1;2,2;2,3];

for i=1:1:5
   fprintf('homework 1,Stiffness qua=%d,mass qua=%d\n',choose(i,1),choose(i,2));
   %num_elements=10;
   %for j=1:1:4
   %    mainFEM1D_Dirichlet(num_elements,1,choose(i,:));
   %    num_elements=2*num_elements;
   %end
   mainErrorAnalysis(1,choose(i,:));
end

for i=1:1:5
   fprintf('homework 2,Stiffness qua=%d,mass qua=%d\n',choose(i,1),choose(i,2));
   %num_elements=10;
   %for j=1:1:4
   %    mainFEM1D_Dirichlet(num_elements,2,choose(i,:));
   %    num_elements=2*num_elements;
   %end
   mainErrorAnalysis(2,choose(i,:));
end




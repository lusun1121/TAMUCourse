function E = mainErrorAnalysis(problem_case,choose)%[problem_case,S]

%problem_case=choose(1);
%S=choose(2);
p = 2; % fe polynomial degree
quad_n_points = 4;%max(choose(1),choose(2));%4;

numIterations = 4; 
baseNumElements = 10; %num elements on first solution
numElements = baseNumElements*2.^(0:(numIterations-1));

E = zeros(numIterations,8);
for i = 1:numIterations
    
   [mesh_size, n_dofs, L2error, H1error, Linferror] = runFEM1D(quad_n_points,numElements(i),...
       p,choose, false, problem_case);
   E(i,:) = [mesh_size, n_dofs, Linferror, 0, L2error, 0, H1error, 0];
   
   % calculate the one step convergence rate using the mesh size, h
 %  if (i>1)
 %      E(i,4) = log(E(i,3)/E(i-1,3))/log(E(i,1)/E(i-1,1)); %Linf error
 %      E(i,6) = log(E(i,5)/E(i-1,5))/log(E(i,1)/E(i-1,1)); %L2 error
 %      E(i,8) = log(E(i,7)/E(i-1,7))/log(E(i,1)/E(i-1,1)); %H1 error
 %  end
   
%       % calculate the one step convergence rate using the N_dofs
%     if (i>1)
%        dim = 1; % space dimension
%        E(i,4) = -1/dim*log(E(i,3)/E(i-1,3))/log(E(i,2)/E(i-1,2)); %Linf error
%        E(i,6) = -1/dim*log(E(i,5)/E(i-1,5))/log(E(i,2)/E(i-1,2)); %L2 error
%        E(i,8) = -1/dim*log(E(i,7)/E(i-1,7))/log(E(i,2)/E(i-1,2)); %H1 error
%     end
   
   % output table of errors and rates

   if (i>1)
       E(i,4) = E(i-1,3)/E(i,3); %Linf error
       E(i,6) = E(i-1,5)/E(i,5); %L2 error
       E(i,8) = E(i-1,7)/E(i,7); %H1 error
       %E(i,4) = log(E(i,3)/E(i-1,3))/log(E(i,1)/E(i-1,1)); %Linf error
       %E(i,6) = log(E(i,5)/E(i-1,5))/log(E(i,1)/E(i-1,1)); %L2 error
       %E(i,8) = log(E(i,7)/E(i-1,7))/log(E(i,1)/E(i-1,1)); %H1 error
   end
   fprintf('%4d & %4d & %1.6e & %1.6e & %1.6f & %1.6e & %1.6f & %1.6e & %1.6f \\\\ \n', E(i,2),(E(i,2)-1)/2+1,E(i,1),E(i,3),E(i,4),E(i,5),E(i,6),E(i,7),E(i,8));
        %fprintf('n_dofs = %4d, h = %1.6e, Einf = %1.6e, Rinf = %1.6f, EL2 = %1.6e, RL2 = %1.6f, EH1 = %1.6e, RH1 = %1.6f \n', E(i,2),E(i,1),E(i,3),E(i,4),E(i,5),E(i,6),E(i,7),E(i,8));

       %fprintf('%4d & %1.2e & %1.2e & %1.2f & %1.2e & %1.2f & %1.2e & %1.2f \\\\ \n', E(i,2),E(i,1),E(i,3),E(i,4),E(i,5),E(i,6),E(i,7),E(i,8));
   %fprintf('n_half_dofs = %4d, n_nofs = %4d, h = %1.6e, Einf = %1.6e, Rinf = %1.6f, EL2 = %1.6e, RL2 = %1.6f, EH1 = %1.6e, RH1 = %1.6f \n', E(i,2),(E(i,2)-1)/2+1,E(i,1),E(i,3),E(i,4),E(i,5),E(i,6),E(i,7),E(i,8));

end

%figure
i=1:1:4;
lmbd=figure;
%plot(E(i,1),E(i,3),'-+',E(i,1),E(i,5),'-o',E(i,1),E(i,7),'-*',E(i,1),(1/2).^i*E(1,3),':');

loglog(E(i,1),E(i,3),'-+',E(i,1),E(i,5),'-o',E(i,1),E(i,7),'-*',E(i,1),(1/8).^(i-1)*(E(1,5)+E(1,3))/2,':',E(i,1),(1/4).^(i-1)*(E(1,3)+E(1,7))/2,'--');
xlabel('mesh');
ylabel('error');
legend('Inf','L2','H1','R=8','R=4');
title({['homework ',num2str(problem_case)],['Stiffness qua=',',',num2str(choose(1)),',','mass qua=',num2str(choose(2))]});
%title({['homework ',num2str(problem_case)],['Stiffness qua=',',',num2str(choose(1)),',','mass qua=',num2str(choose(2))]});
%    saveas(lmbd, ['\\blender\homes\l\u\lusun8825\nt\AccountSettings\Desktop\610\FEM1D_For_Students\sunlu\result\','EX1_hm',num2str(problem_case),...
%        '_S',num2str(choose(1)),'_M',num2str(choose(2)),...%'_N',num2str(DoFHandler.n_dofs),
%        '.png'],'png');

end


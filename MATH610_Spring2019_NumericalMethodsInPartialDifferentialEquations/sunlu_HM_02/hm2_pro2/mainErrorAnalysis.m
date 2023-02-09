function [mesh,Inf_error] =mainErrorAnalysis(method,epsillion,numIterations,numElements)


E = zeros(numIterations,5);
for i = 1:numIterations
    
   [method, mesh_size, n_dofs, Linferror,uh] = runFEM1D( numElements(i),method,epsillion);
   E(i,:) = [method , mesh_size, n_dofs, Linferror, 0 ];
   
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
       E(i,5) = log(E(i,4)/E(i-1,4))/log(E(i,2)/E(i-1,2)); %Linf error
   end
   fprintf('%1d & %4d & %4d & %1.2e & %1.6e & %1.6e \\\\ \n', E(i,1),E(i,3)-1,E(i,3),E(i,2),E(i,4),E(i,5));
   %fprintf('method=%1d, n_mesh = %4d, n_dofs = %4d, h = %1.2e, Einf = %1.6e, Rinf = %1.6f \n', E(i,1),E(i,3)-1,E(i,3),E(i,2),E(i,4),E(i,5));

end
mesh = E(:,2);
Inf_error = E(:,4);
end




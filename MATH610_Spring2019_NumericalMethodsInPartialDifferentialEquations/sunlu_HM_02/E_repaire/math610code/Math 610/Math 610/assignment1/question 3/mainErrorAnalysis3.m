
p = 1; % fe polynomial degree
quad_n_points = 4;
numIterations = 6; 
baseNumElements = 40; %num elements on first solution
numElements = baseNumElements*2.^(0:(numIterations-1));

E = zeros(numIterations,8);
for i = 1:numIterations
    
   [mesh_size, n_dofs, L2error, H1error, Linferror] = runFEM1D2(numElements(i), p,quad_n_points, false,u_exact,grad_u_exact);
   E(i,:) = [mesh_size, n_dofs, Linferror, 0, L2error, 0, H1error, 0];
   
   % calculate the one step convergence rate using the mesh size, h
   if (i>1)
       E(i,4) = log(E(i,3)/E(i-1,3))/log(E(i,1)/E(i-1,1)); %Linf error
       E(i,6) = log(E(i,5)/E(i-1,5))/log(E(i,1)/E(i-1,1)); %L2 error
       E(i,8) = log(E(i,7)/E(i-1,7))/log(E(i,1)/E(i-1,1)); %H1 error
   end
    
%       % calculate the one step convergence rate using the N_dofs
%     if (i>1)
%        dim = 1; % space dimension
%        E(i,4) = -1/dim*log(E(i,3)/E(i-1,3))/log(E(i,2)/E(i-1,2)); %Linf error
%        E(i,6) = -1/dim*log(E(i,5)/E(i-1,5))/log(E(i,2)/E(i-1,2)); %L2 error
%        E(i,8) = -1/dim*log(E(i,7)/E(i-1,7))/log(E(i,2)/E(i-1,2)); %H1 error
%     end
   
   % output table of errors and rates
   fprintf('n_dofs = %4d, h = %1.2e, Einf = %1.2e, Rinf = %1.2f, EL2 = %1.2e, RL2 = %1.2f, EH1 = %1.2e, RH1 = %1.2f \n', E(i,2),E(i,1),E(i,3),E(i,4),E(i,5),E(i,6),E(i,7),E(i,8));
   
end



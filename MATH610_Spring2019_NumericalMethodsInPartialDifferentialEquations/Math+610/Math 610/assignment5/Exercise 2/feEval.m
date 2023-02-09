function [FE_at_quad] = feEval( Quad, p )

x = Quad.xhat(:,1);
y = Quad.xhat(:,2);

if (p == 1)
    %                   phi_1, phi_2, phi_3
    % Linear elements :   1-x-y,  x,   y
    % x derivative    :   -1,   1
    FE_at_quad.hat_phi  = [1-x-y, x, y]; %nq * 3
%     FE_at_quad.hat_phix1 = [-1*ones(Quad.nq,1), -1*ones(Quad.nq,1)];
%     FE_at_quad.hat_phix2 = [1*ones(Quad.nq,1), zeros(Quad.nq,1)];
%     FE_at_quad.hat_phix3 = [zeros(Quad.nq,1), 1*ones(Quad.nq,1)];
    FE_at_quad.hat_phix = [-1*ones(Quad.nq,1), 1*ones(Quad.nq,1), zeros(Quad.nq,1)]; %nqx3
    FE_at_quad.hat_phiy = [-1*ones(Quad.nq,1), zeros(Quad.nq,1), 1*ones(Quad.nq,1)]; %nqx3

elseif (p==2)
   % construct 3 basis functions on [0,1]:  phi_0,  phi_1,  phi_{1/2}
   
   % construct 3 x-derivatives of basis functions on [0,1]:  dphi_0/dx, dphi_1/dx,  dphi_{1/2}/dx
    
elseif (p==3)
   % construct 4 basis functions on [0,1]:  phi_0,  phi_1,  phi_{1/3}, phi_{2/3}
   
   % construct 4 x-derivatives of basis functions on [0,1]:  dphi_0/dx, dphi_1/dx,  dphi_{1/3}/dx , dphi_{2/3}/dx
    
end




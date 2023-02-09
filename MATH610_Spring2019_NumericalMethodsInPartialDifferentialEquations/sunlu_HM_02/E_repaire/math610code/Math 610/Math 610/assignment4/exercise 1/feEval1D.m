
function [FE_at_quad] = feEval1D( Quad, p )
x = Quad.xhat(:,1);
if (p == 1)
    FE_at_quad.hat_phi  = [1-x, x];
    
end
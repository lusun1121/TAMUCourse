function [mean]=meanatelement(localvertices, uh_at_localdofs, Quad, FE_at_Quad)
mat_B = [(localvertices(2,:) - localvertices(1,:))',(localvertices(3,:) - localvertices(1,:))']; %[(dim)x(dim)]
det_B = det(mat_B); 

uh_at_q_points = FE_at_Quad.hat_phi*uh_at_localdofs;

mean=uh_at_q_points'*Quad.what*det_B;


function [grad]=grad_at_element(localvertices,uh_at_localdofs)
mat_B = [(localvertices(2,:) - localvertices(1,:))',(localvertices(3,:) - localvertices(1,:))']; %[(dim)x(dim)]
det_B = det(mat_B); 
inv_B = 1/det_B*[mat_B(2,2),-mat_B(1,2);-mat_B(2,1),mat_B(1,1)];
gradmatrix=[-1,1,0;-1,0,1];%2x3
grad=(gradmatrix*uh_at_localdofs)'*inv_B;%1*2
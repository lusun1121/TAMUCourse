function A_matrix = Differnece_One(A, T ,options,b)
% options    finite  defferences method
%   1              unpwind
%   2              downwind
%   3              centered
h=T.meshsize;
switch (options)
    case 1
        a_ij = [-1,2,-1]/(h^2) + [-1,1,0]*b/h;
        for cell = 2:(T.n_nodes-1)
            dofIndices =  (cell-1): (cell+1);
            A(cell,dofIndices) = a_ij;
        end
    case 2
        a_ij = [-1,2,-1]/(h^2) + [0,-1,1]*b/h;
        for cell = 2:(T.n_nodes-1)
            dofIndices =  (cell-1): (cell+1);
            A(cell,dofIndices) = a_ij;
        end
    case 3
        a_ij = [-1,2,-1]/(h^2) + [-1,0,1]*b/h/2;
        for cell = 2:(T.n_nodes-1)
            dofIndices =  (cell-1): (cell+1);
            A(cell,dofIndices) = a_ij;
        end
end

A_matrix=A;

end


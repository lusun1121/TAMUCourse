function elementMeshSizes = computeMeshSizes(elements, nodes)

n_ele = size(elements,1);

elementMeshSizes = zeros(n_ele,1);

for cell = 1:n_ele
    vertices = nodes(elements(cell,:),:);
    hK = 0.5*det([ones(1,3); vertices']);  % calculate the area of triangle
    elementMeshSizes(cell) = hK;
end
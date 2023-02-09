function T = TriangleReader(filename_nodes, filename_elements, filename_edges)
% T = TriangleReader(filename_nodes, filename_elements, filename_edges)
%
%  input:
%    filename_nodes = string filename of .node file
%    filename_elements = string filename of .ele file
%    optional:
%      filename_edges = string filename of .edge file
%  output:
%    T = struct of triangulation components
%
%  example: Suppose we have A.1.node and A.1.ele in my current directory,
%  then calling
%  
%   T = TriangleReader('A.1.node','A.1.ele');
%
%  extracts the triangulation information and returns it as T.
%
%  Note: You can access elements of A by using '.' notation
% 
% For example,
%   A.nodes = list of nodes
%   A.nodes(5,1) = x coordinate of fifth node
%   A.nodes(5,2) = y coordinate of fifth node etc
%
%   A.elements = list of elements
%
% Caution: The below code assumes that the .node and .ele files that are 
% passed in have no comments at the beginning of the files.  This is
% consistent with the way they are generated by Triangle, but you must
% ensure that this is the case or change below to accomodate comments.
%
% Possible extensions that could later be useful: 
%
% 1. Add in the possibility of a filename_edges
% 'A.1.edge' file and compile a list of all edges and of boundary edges
% which can be added to the struct below with the integer n_edges.
%
% 2. Given a domain with mixed dirichlet boundary conditions and neumann
% boundary conditions, construct two lists for boundary nodes and edges.
% That is DirichletBoundaryNodes,  NeumannBoundaryNodes,
% DirichletBoundaryEdges, NeumannBoundaryEdges and add them to struct
% below.  These can be identified by setting/reading the boundary indicator
% or material id values in the poly file.
%
% 3. Add in a list of neighbors of triangles (.neigh file) from which you 
% could concievable construct a list of neighbors across edges and degree 
% of freedom support patches around each vertex.  This would be good for
% things like the clement interpolant or an edge jump calculation say for
% an a posteriori estimator for adaptive refinement.
%

read_in_edge_file = false;
if nargin > 2
    read_in_edge_file = true;
end

%
% Read in node information
%
fid = fopen(filename_nodes);

% assume no comments at beginning of node file
n_nodes = fscanf(fid, '%d', 1);
dim = fscanf(fid,'%d', 1);
n_materialproperties = fscanf(fid,'%d', 1);
b_useboundaryindicators = fscanf(fid,'%d',1);

nodes = zeros(n_nodes,dim);

if (n_materialproperties > 0)
	materialproperties = zeros(n_nodes,n_materialproperties);
else
	materialproperties = [];
end


% define coordinate patterns by dimension
if (dim == 1)
  Pattern = '%f';
elseif (dim == 2)
  Pattern = '%f %f';
elseif (dim == 3)
  Pattern = '%f %f %f';
end

boundarynodes = zeros(n_nodes,1);
interiornodes = zeros(n_nodes,1);
	
n_boundarynodes = 0;
n_interiornodes = 0;

% read off node coordinates and properties
for n = 1:n_nodes
    fscanf(fid,'%d',1); % do nothing for first integer of line
	nodes(n,:) = fscanf(fid,Pattern,dim);
	for m = 1:n_materialproperties
		materialproperties(n,m) = fscanf(fid,'%f',1);
	end
	if (b_useboundaryindicators)
        isboundarynode = fscanf(fid,'%d',1);
        if (isboundarynode > 0.5)
            n_boundarynodes = n_boundarynodes + 1;
            boundarynodes(n_boundarynodes,1) = n;  % add pointer to list of boundary nodes;
        else
            n_interiornodes = n_interiornodes + 1;
            interiornodes(n_interiornodes,1) = n; % add pointer to list of interior nodes;
        end
	end
end

% trim extra space
boundarynodes = boundarynodes(1:n_boundarynodes,1);
interiornodes = interiornodes(1:n_interiornodes,1);

fclose(fid);  % close node file



%
% Read in element information 
%

fid = fopen(filename_elements);

% assume no comments at beginning of element file
n_elements = fscanf(fid, '%d', 1);
n_nodes_per_element = fscanf(fid,'%d', 1);
n_elementattributes = fscanf(fid,'%d', 1);

elements = zeros(n_elements,n_nodes_per_element);
if (n_elementattributes > 0)
  elementattributes = zeros(n_elements,n_elementattributes);
else
  elementattributes = [];
end

% read off node coordinates and properties
for el = 1:n_elements
    fscanf(fid,'%d',1); % do nothing for first integer of line
	for m = 1:n_nodes_per_element
		elements(el,m) = fscanf(fid,'%d', 1);
	end
	
	for m = 1:n_elementattributes
		elementattributes(el,m) = fscanf(fid,'%f',1);
	end
end
fclose(fid);



if read_in_edge_file
    %
    % Read in edge information 
    %

    fid = fopen(filename_edges);

    % assume no comments at beginning of edge file
    n_edges = fscanf(fid, '%d', 1);
    n_nodes_per_edge = 2;
    b_useboundaryedgemarkers = fscanf(fid,'%d',1);

    edges = zeros(n_edges,n_nodes_per_edge);
    boundaryedges = zeros(n_edges,1);
    interioredges = zeros(n_edges,1);
    
    n_boundaryedges = 0;
    n_interioredges = 0;
    
    % read off edges
    for e = 1:n_edges
        fscanf(fid,'%d',1); % do nothing for first integer of line
        
        for m = 1:n_nodes_per_edge
            edges(e,m) = fscanf(fid,'%d', 1);
        end
        if (b_useboundaryedgemarkers)
            isboundaryedge = fscanf(fid,'%d',1);
            if (isboundaryedge > 0.5)
                n_boundaryedges = n_boundaryedges + 1;
                boundaryedges(n_boundaryedges,1) = e;  % add pointer to list of boundary nodes;
            else
                n_interioredges = n_interioredges + 1; 
                interioredges(n_interioredges,1) = e;  % add pointer to list of interior nodes;
            end
        end
    end
    
    % trim extra space allocated
    boundaryedges = boundaryedges(1:n_boundaryedges,1);
    interioredges = interioredges(1:n_interioredges,1);
    
    fclose(fid);
    
end


%
% Construct struct of mesh information
%
T=struct('dim',dim,...
         'n_nodes',n_nodes,...
         'nodes',nodes,...
         'n_boundarynodes',n_boundarynodes,...
         'boundarynodes',boundarynodes,...
         'n_interiornodes',n_interiornodes,...
         'interiornodes',interiornodes,...
         'n_materialproperties',n_materialproperties,...
         'materialproperties',materialproperties,...
         'n_elements',n_elements,...
         'elements',elements,...
         'n_elementattributes',n_elementattributes,...
         'elementattributes',elementattributes); 
     
 if read_in_edge_file
    T.('n_edges') = n_edges;
    T.('edges') = edges;
    T.('n_boundaryedges') = n_boundaryedges;
    T.('boundaryedges') = boundaryedges;
    T.('n_interioredges') = n_interioredges;
    T.('interioredges') = interioredges;
 
 end

% if you want to add in the edges flag for dirichlet, neumann or robin bc
% you should create a new function that takes this general Triangle struct T
% and adds a new field called T.('edgeFlags') which would be a specific
% function based on geometry and checking both end points of edge to determin
% whether an edge is on boundary or not...

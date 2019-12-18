function [ind, int_func] = mytsearchn(mesh, coord)
% MYTSEARCHN carries out a spatial point location search within the MESH. 
% 
% [IND, INT_FUNC] = MYTSEARCHN(MESH, COORD) Determines which, if any,
%   MESH element a given point (COORD) falls in: 
%     MESH - FEM mesh in NIRFAST format, can be 2D (triangels) or 3D
%      (tetrahedrons)
%     COORD - n by m matrix of n cartesian coordinates to search.
%      m=2 for 2D (x,y) and m=3 for 3D (x,y,z) meshes. 
%     IND - n by 1 list of indexes of MESH elements where the n-th COORD
%      belongs. NaN if the point is outside of the MESH. 
%     INT_FUNC - n by 3 for 2D MESH or n by 4 for 3D MESH barycentric
%      coordinates of n-th COORD. NaN if the point is outside of the
%      MESH. 
% 
% First, it uses triangulation search algorithm that locates the
% triange/tetrahedron enclosing a query point as implemented in
% MYTSEARCHN_FAST. For all points that the triangulation algorithm
% determine outside the MESH, a search based on a 'volume' (3D) or a
% 'surface' (2D) criterion is used.
% 
% The 'volume' criterion is satisfied if a difference between an element
% volume and four volumes spanned on the COORD point and the tested
% tetrahedral element vertices is at leaset 10^6 times smaller than the
% volume of the element itself. As implemented in INSIDE_TETRAHEDRON.
% 
% The 'area' criterion is satisfied if a difference between an element
% surface and three surfaces spanned on the COORD point and the tested
% triangle element vertices is at leaset 10^6 times smaller than the
% surface of the element itself. As implemented in INSIDE_TRIANGLE.
% 
% The 'volume' and 'area' methods are used as might be required by other
% NIRFAST functions.
% 
% This function also compensates for numerical errors ensuring that sum of
% the barycentric coordinates (INT_FUNC) for a given MESH element is 1. 
% 
% See also MYTSEARCHN_FAST, INSIDE_TRIANGLE, INSIDE_TETRAHEDRON, TSEARCHN,
%          TRIANGULATION, POINTLOCATION. 
% 
%   Part of NIRFAST package.
%   H. Dehghani and S. Wojtkiewicz 2018

%% check in/out

narginchk(2,2);
nargoutchk(2,2);


%% TRY TO USE FAST FIRST

% Use fast first. This will return NaNs if tested point does not belong to
% the mesh.

% turn off warning that a coordinate is outside the mesh, we will handle this later 
warning('off','MATLAB:triangulation:PtsNotInTriWarnId')
% call the fast triangulation search
[ind, int_func] = mytsearchn_fast(mesh, coord);
% turn on the warning
warning('on','MATLAB:triangulation:PtsNotInTriWarnId')

% if all checked point inside the mesh
if sum(~isfinite(ind)) == 0
    % all done here
    return;
end

%% Now try the second method based on the 'volume' or 'surface' criterion

% This method is used as might be required by other NIRFAST functions.

% make copies and mask to merge this and fast methods at the end
ind_fast = ind;
int_func_fast = int_func;

mask_coords_outside = ~isfinite(ind_fast); 
coord = coord(mask_coords_outside,:);

% The 'volume' or 'surface' criterion starts here

if mesh.dimension == 2
    
    N = size(coord,1);
    ind = NaN(N,1); int_func = NaN(N,3);
    for i = 1:N
        
        % determine distance of all nodes to coord.  This applies to source
        % positions or when using this program to create regionized
        % meshes
        dist = sum((mesh.nodes(:,1:mesh.dimension) - coord(i,1:mesh.dimension)).^2,2);
%         dist = sqrt(dist);
        
        % sort nodes in order of nearest to farthest
        [~,temp_ind] = sort(dist);
        
        % Start with nearest node and test surface elements of which the node is a
        % vertex
        j = 1; true = 0;
        
        % find elements with node 'x' as a vertex
        while (true == 0) && (j <= 10)
            % temp_ind is list of node numbers in order of distance from
            % det point.  Find the elements which contain this node:
            [r,~] = find(mesh.elements == temp_ind(j));
            
            % foo is matrix of all nodes connected to node 'x', in node
            % connectivity list form.  We'll call it a 'shortened connectivity list'.
            % Contains both surface and internal nodes.
            foo = mesh.elements(r,:);
            [n,~] = size(foo);
            
            % Try each element to see if coord is inside
            k = 1; true = 0;
            while (true == 0) && (k <= n)
                % To make the syntax a little easier to read, define points P, Q, R - vertices of the surface triangle
                % which we are testing
                P = mesh.nodes(foo(k,1),1:mesh.dimension); Q = mesh.nodes(foo(k,2),1:mesh.dimension); R = mesh.nodes(foo(k,3),1:mesh.dimension);
                
                % check to see if intersection point is within element:
                true = inside_triangle(coord(i,1:mesh.dimension),P,Q,R);
                if true == 1
                    % if is in element, store element number
                    ind(i,1) = r(k);
                    % Calculate barycentric coordinate of coord in
                    % triangular element:  This is the integrating function.
                    A = [P(1) Q(1) R(1); P(2) Q(2) R(2); 1 1 1];
                    b = [coord(i,1); coord(i,2); 1];
                    int_func(i,:) = (A\b)';
                    % correction for numerical errors
                    int_func(i,end) = 1 - sum(int_func(i,1:end-1),2);
                elseif true == 0
                    k = k + 1;
                end
            end
            j = j + 1;
        end
    end
elseif mesh.dimension == 3
    
    N = size(coord,1);
    ind = NaN(N,1); int_func = NaN(N,4);
    for i = 1:N
        
        % determine distance of all nodes to coord.  This applies to source
        % positions or when using this program to create regionized
        % meshes
        dist = sum((mesh.nodes(:,1:mesh.dimension) - coord(i,1:mesh.dimension)).^2,2);
%         dist = sqrt(dist);
        
        % sort nodes in order of nearest to farthest
        [~,temp_ind] = sort(dist);
        
        % Start with nearest node and test surface elements of which the node is a
        % vertex
        j = 1; true = 0;
        
        % find elements with node 'x' as a vertex
        while true == 0 && j <= 10
            % temp_ind is list of node numbers in order of distance from
            % det point.  Find the elements which contain this node:
            [r,~] = find(mesh.elements == temp_ind(j));
            
            
            % foo is matrix of all nodes connected to node 'x', in node
            % connectivity list form.  We'll call it a 'shortened connectivity list'.
            % Contains both surface and internal nodes.
            foo = mesh.elements(r,:);
            [n,~] = size(foo);
            
            % Try each element to see if coord is inside
            k = 1; true = 0;
            while (true == 0) && (k <= n)
                % To make the syntax a little easier to read, define points P, Q, R, S
                %   - vertices of the tetrahedron which we are testing
                P = mesh.nodes(foo(k,1),1:mesh.dimension);
                Q = mesh.nodes(foo(k,2),1:mesh.dimension);
                R = mesh.nodes(foo(k,3),1:mesh.dimension);
                S = mesh.nodes(foo(k,4),1:mesh.dimension);
                
                % check to see if intersection point is within element:
                true = inside_tetrahedron(coord(i,1:mesh.dimension),P,Q,R,S);
                if true == 1
                    % if det is in element, store element number
                    ind(i,1) = r(k);
                    % Calculate barycentric coordinate of coord in
                    % triangular element:  This is the integrating function.
                    A = [P(1) Q(1) R(1) S(1); P(2) Q(2) R(2) S(2); P(3) Q(3) R(3) S(3); 1 1 1 1];
                    b = [coord(i,1); coord(i,2); coord(i,3); 1];
                    int_func(i,:) = (A\b)';
                    % correction for numerical errors
                    int_func(i,end) = 1 - sum(int_func(i,1:end-1),2);
                elseif true == 0
                    k = k + 1;
                end
            end
            j = j + 1;
        end
    end
end

%% Merge

ind_fast(mask_coords_outside) = ind;
int_func_fast(mask_coords_outside,:) = int_func;

ind = ind_fast;
int_func = int_func_fast;


end

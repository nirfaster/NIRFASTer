function [ind, int_func] = mytsearchn_fast(mesh, coord)
% MYTSEARCHN_FAST carries out a spatial point location search. It uses
%   triangulation search algorithm that locates the triange/tetrahedron
%   enclosing a query point. This is an interface for NIRFAST mesh format.
%   This function also compensates for numerical errors ensuring that sum of
%   the barycentric coordinates (INT_FUNC) for a given MESH element is 1.
% 
%   [IND, INT_FUNC] = MYTSEARCHN_FAST(MESH, COORD) Determines which, if any,
%   MESH element a given point (COORD) falls in: 
%       MESH - FEM mesh in NIRFAST format, can be 2D (triangels) or 3D
%        (tetrahedrons)
%       COORD - n by m matrix of n cartesian coordinates to search.
%        m=2 for 2D (x,y) and m=3 for 3D (x,y,z) meshes. 
%       IND - n by 1 list of indexes of MESH elements where the n-th COORD
%        belongs. NaN if the point is outside of the MESH. 
%       INT_FUNC - n by 3 for 2D MESH or n by 4 for 3D MESH barycentric
%        coordinates of n-th COORD. NaN if the point is outside of the
%        MESH. 
% 
% IMPORTANT:
%   MYTSEARCHN_FAST reurns NaNs if COORD does not belong to the mesh
%   volume.
%   If you want to try and 'snap' the COORD to the MESH elements, please use
%   the MYTSEARCHN that uses 'coarse' test conditions.
% 
% See also MYTSEARCHN, TSEARCHN, TRIANGULATION, POINTLOCATION.
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

%% check in/out

narginchk(2,2);
nargoutchk(0,2);

%% BODY

[ind, int_func] = pointLocation(triangulation(mesh.elements, mesh.nodes(:,1:mesh.dimension)), coord(:,1:mesh.dimension));

% correction for numerical errors
int_func(:,end) = 1 - sum(int_func(:,1:end-1),2);

end



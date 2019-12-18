function [edges] = nodal_distance_tetrahedron(a,b,c,d)
% NODAL_DISTANCE_TETRAHEDRON Calcultes tetrahedron edges length.
% 
% SYNTAX:
%  [EDGES] = NODAL_DISTANCE_TETRAHEDRON(A,B,C,D)
%  NODAL_DISTANCE_TETRAHEDRON(A,B,C,D)
% 
% [EDGES] = NODAL_DISTANCE_TETRAHEDRON(A,B,C,D)
%   - A, B, C and D - are m by 3 matrices of (x,y,z) cartesian coordinates
%                     of vertices of m tetrahedrons.
%   - EDGES - is m by 6 matrix of m tetrahedrons edges length. Expressed
%             in the A, B, C and D coordinates unit. 
% 
% See also volume_tetrahedron, height_tetrahedron, facets_area_tetrahedron. 
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

%% check in/out

narginchk(4,4);
nargoutchk(0,1);

%% BODY

edges = sqrt([sum((b-a).^2,2) sum((c-a).^2,2) sum((d-a).^2,2),...
              sum((c-b).^2,2) sum((b-d).^2,2) sum((c-d).^2,2)]);

end


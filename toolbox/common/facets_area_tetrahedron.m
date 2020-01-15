function [S] = facets_area_tetrahedron(a,b,c,d)
% FACETS_AREA_TETRAHEDRON Calcultes tetrahedron facets area.
% 
% SYNTAX:
%  [S] = FACETS_AREA_TETRAHEDRON(A,B,C,D)
%  FACETS_AREA_TETRAHEDRON(A,B,C,D)
% 
%  [S] = FACETS_AREA_TETRAHEDRON(A,B,C,D)
%   - A, B, C and D - are m by 3 matrices of (x,y,z) cartesian coordinates
%                     of vertices od m tetrahedrons.
%   - S - is m by 4 matrix of m tetrahedrons facets area. Expressed in the
%         A, B, C and D coordinates unit.
% 
% See also nodal_distance_tetrahedron, volume_tetrahedron, height_tetrahedron. 
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

%% check in/out

narginchk(4,4);
nargoutchk(0,1);

%% BODY

S = [sqrt(sum(cross(b - a,c - a,2).^2,2))/2,...
     sqrt(sum(cross(a - d,b - d,2).^2,2))/2,...
     sqrt(sum(cross(b - d,c - d,2).^2,2))/2,...
     sqrt(sum(cross(a - d,c - d,2).^2,2))/2];

end
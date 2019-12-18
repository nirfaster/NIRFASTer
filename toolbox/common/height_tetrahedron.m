function h = height_tetrahedron(a,b,c,d)
% HEIGHT_TETRAHEDRON Calcultes tetrahedron heights.
% 
% SYNTAX:
%  [H] = HEIGHT_TETRAHEDRON(A,B,C,D)
%  HEIGHT_TETRAHEDRON(A,B,C,D)
% 
%  [H] = HEIGHT_TETRAHEDRON(A,B,C,D)
%   - A, B, C and D - are m by 3 matrices of (x,y,z) cartesian coordinates
%                     of vertices od m tetrahedrons.
%   - H - is m by 4 matrix of m tetrahedrons heights. Expressed in the A,
%         B, C and D coordinates unit.
% 
% See also nodal_distance_tetrahedron, volume_tetrahedron, facets_area_tetrahedron. 
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

%% check in/out

narginchk(4,4);
nargoutchk(0,1);

%% BODY

% h = 3 * V / S;

V = (abs(dot(a - d, cross(b - d,c - d,2),2))/6);

h = 3*[V./(sqrt(sum(cross(b - a,c - a,2).^2,2))/2),...
       V./(sqrt(sum(cross(a - d,b - d,2).^2,2))/2),...
       V./(sqrt(sum(cross(b - d,c - d,2).^2,2))/2),...
       V./(sqrt(sum(cross(a - d,c - d,2).^2,2))/2)];
end
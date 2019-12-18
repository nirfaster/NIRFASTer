function V = volume_tetrahedron(a,b,c,d)
% VOLUME_TETRAHEDRON Calcultes tetrahedron volume.
% 
% SYNTAX:
%  [V] = VOLUME_TETRAHEDRON(A,B,C,D)
%  VOLUME_TETRAHEDRON(A,B,C,D)
% 
%  [V] = VOLUME_TETRAHEDRON(A,B,C,D)
%   - A, B, C and D - are m by 3 matrices of (x,y,z) cartesian coordinates
%                     of vertices od m tetrahedrons.
%   - V - is m elements vector of m tetrahedrons volume. Expressed in the
%         A, B, C and D coordinates unit.
% 
% See also nodal_distance_tetrahedron, height_tetrahedron, facets_area_tetrahedron. 
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

%% check in/out

narginchk(4,4);
nargoutchk(0,1);

%% BODY

V = abs(dot(a - d, cross(b - d,c - d,2),2))/6;

end
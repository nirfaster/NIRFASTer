function h = height_triangle(a,b,c)
% HEIGHT_TRIANGLE Calcultes triangles heights.
% 
% SYNTAX:
%  [H] = HEIGHT_TRIANGLE(A,B,C)
%  HEIGHT_TRIANGLE(A,B,C)
% 
%  [H] = HEIGHT_TRIANGLE(A,B,C)
%   - A, B and C - are m by 2 matrices of (x,y) cartesian coordinates
%                  of vertices of m tetrahedrons.
%   - H - is m by 3 matrix of m triangles heights. Expressed in the A,
%         B and C coordinates unit.
% 
% See also nodal_distance_triangle, surface_triangle. 
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

%% check in/out

narginchk(3,3);
nargoutchk(0,1);

%% BODY

% h = 2* area / edge;

h = 2 * surface_triangle(a,b,c) ./ nodal_distance_triangle(a,b,c);

   
end
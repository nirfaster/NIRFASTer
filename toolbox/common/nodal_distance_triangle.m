function [edges] = nodal_distance_triangle(a,b,c)
% NODAL_DISTANCE_TRIANGLE Calcultes triangles edges length.
% 
% SYNTAX:
%  [EDGES] = NODAL_DISTANCE_TRIANGLE(A,B,C)
%  NODAL_DISTANCE_TRIANGLE(A,B,C)
% 
% [EDGES] = NODAL_DISTANCE_TRIANGLE(A,B,C)
%   - A, B and C - are m by 2 matrices of (x,y) cartesian coordinates of
%                  vertices of m triangles. 
%   - EDGES - is m by 3 matrix of m ttriangles edges length. Expressed
%             in the A, B and C coordinates unit.
% 
% See also surface_triangle, height_triangle.
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

%% check in/out

narginchk(3,3);
nargoutchk(0,1);

%% BODY

edges = sqrt([sum((b-a).^2,2) sum((c-a).^2,2) sum((c-b).^2,2)]);

end


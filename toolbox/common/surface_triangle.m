function S = surface_triangle(a,b,c)
% SURFACE_TRIANGLE Calcultes riangles surface.
% 
% SYNTAX:
%  [S] = SURFACE_TRIANGLE(A,B,C)
%  SURFACE_TRIANGLE(A,B,C)
% 
%  [S] = SURFACE_TRIANGLE(A,B,C)
%   - A, B and C - are m by 2 matrices of (x,y) cartesian coordinates
%                  of vertices of m triangles.
%   - S - is m elements vector of m triangles surface. Expressed in the
%         A, B and C coordinates unit.
% 
% See also nodal_distance_triangle, height_triangle. 
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

%% check in/out

narginchk(3,3);
nargoutchk(0,1);

%% BODY

b = b - a;
c = c - a;
S = abs(b(:,1).*c(:,2) - b(:,2).*c(:,1))/2;

end
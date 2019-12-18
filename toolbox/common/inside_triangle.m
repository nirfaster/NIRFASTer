function [True] = inside_triangle(P,P1,P2,P3)
% INSIDE_TRIANGLE checks if a 2D point belongs to a 2D triangle usign the
% 'area' criterion.
% 
%   TRUE = INSIDE_TRIANGLE(P,P1,P2,P3) Tests if point P belongs to the
%   triangle spanned on P1, P2 and P3. All points are cartesian
%   coordinates.
%    P - the point to test
%    P1 - the triangel vertex 1
%    P2 - the triangel vertex 2
%    P3 - the triangel vertex 3
%    TRUE - 1 - inside, 0 - outside the triangle
%
% The 'area' criterion is satisfied if a difference between an element
% surface and three surfaces spanned on the COORD point and the tested
% triangle element vertices is at leaset 10^6 times smaller than the
% surface of the element itself.
%
% IMPORTANT:
%   This function accepts points in at least 2 dimensions. However,
%   dimensions >2 are ignored.
% 
%   Part of NIRFAST package.

%% check in/out

narginchk(4,4);
nargoutchk(1,1);

%% BODY
Area_P1P2P3 = 1/2. *abs(det([P1(1) P1(2) 1;P2(1) P2(2) 1;P3(1) P3(2) 1]));

if abs(1/2. *((abs(det([P(1) P(2) 1;P1(1) P1(2) 1;P2(1) P2(2) 1])))...
        + (abs(det([P(1) P(2) 1;P2(1) P2(2) 1;P3(1) P3(2) 1])))...
        + (abs(det([P(1) P(2) 1;P3(1) P3(2) 1;P1(1) P1(2) 1])))...
        )-Area_P1P2P3)/Area_P1P2P3 < 10^-6
    True = 1;
else
    True = 0;
end

end
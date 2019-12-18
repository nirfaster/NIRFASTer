function [True] = inside_tetrahedron(P,P1,P2,P3,P4)
% INSIDE_TETRAHEDRON checks if a 3D point belongs to a tetrahedron usign
% the 'volume' criterion.
% 
%   TRUE = INSIDE_TETRAHEDRON(P,P1,P2,P3,P4) Tests if point P belongs to the
%   tetrahedron spanned on P1, P2, P3 and P4. All points are cartesian
%   coordinates.
%    P - the point to test
%    P1 - the tetrahedron vertex 1
%    P2 - the tetrahedron vertex 2
%    P3 - the tetrahedron vertex 3
%    P4 - the tetrahedron vertex 4
%    TRUE - 1 - inside, 0 - outside the tetrahedron
%
% The 'volume' criterion is satisfied if a difference between an element
% volume and four volumes spanned on the COORD point and the tested
% tetrahedral element vertices is at leaset 10^6 times smaller than the
% volume of the element itself.
%
% IMPORTANT:
%   This function accepts points in at least 3 dimensions. However,
%   dimensions >3 are ignored.
% 
%   Part of NIRFAST package.

%% check in/out

narginchk(5,5);
nargoutchk(1,1);

%% BODY

V0 = 1/6.*abs(det([P1(1) P1(2) P1(3) 1;P2(1) P2(2) P2(3) 1;P3(1) P3(2) P3(3) 1;P4(1) P4(2) P4(3) 1]));

if abs(1/6.*((abs(det([P(1) P(2) P(3) 1;P2(1) P2(2) P2(3) 1;P3(1) P3(2) P3(3) 1;P4(1) P4(2) P4(3) 1])))...
        + (abs(det([P1(1) P1(2) P1(3) 1;P(1) P(2) P(3) 1;P3(1) P3(2) P3(3) 1;P4(1) P4(2) P4(3) 1])))...
        + (abs(det([P1(1) P1(2) P1(3) 1;P2(1) P2(2) P2(3) 1;P(1) P(2) P(3) 1;P4(1) P4(2) P4(3) 1])))...
        + (abs(det([P1(1) P1(2) P1(3) 1;P2(1) P2(2) P2(3) 1;P3(1) P3(2) P3(3) 1;P(1) P(2) P(3) 1])))...
        ) - V0)/V0 < 10^-6
    True = 1;
else
    True = 0;
end

end
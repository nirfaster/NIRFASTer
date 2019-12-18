function [solver_name] = solver_name_backslash
% SOLVER_NAME_BACKSLASH Defines MATLAB '\' solver name.
%
% [SOLVER_OUT] = SOLVER_NAME_BACKSLASH  Always returns 'backslash'.
%
% See also SOLVER_NAME_CPU, SOLVER_NAME_GPU, SOLVER_NAME_MATLAB_ITERATIVE,
%          isCUDA.
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

%% check in/out

narginchk(0,0);
nargoutchk(0,1);

%% BODY

solver_name = 'backslash';

end
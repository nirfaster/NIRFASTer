function [solver_name] = solver_name_CPU
% SOLVER_NAME_CPU Defines CPU solver name.
%
% [SOLVER_OUT] = SOLVER_NAME_CPU  Always returns 'BiCGStab_CPU'.
%
% See also SOLVER_NAME_GPU, SOLVER_NAME_MATLAB_ITERATIVE,
%          SOLVER_NAME_BACKSLASH, isCUDA.
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

%% check in/out

narginchk(0,0);
nargoutchk(0,1);

%% BODY

solver_name = 'BiCGStab_CPU';

end
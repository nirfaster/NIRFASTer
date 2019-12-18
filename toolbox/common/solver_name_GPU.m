function [solver_name] = solver_name_GPU
% SOLVER_NAME_GPU Defines GPU solver name.
%
% [SOLVER_OUT] = SOLVER_NAME_GPU  Always returns 'BiCGStab_GPU'.
%
% See also SOLVER_NAME_CPU, SOLVER_NAME_MATLAB_ITERATIVE,
%          SOLVER_NAME_BACKSLASH, isCUDA. 
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

%% check in/out

narginchk(0,0);
nargoutchk(0,1);

%% BODY

solver_name = 'BiCGStab_GPU';

end
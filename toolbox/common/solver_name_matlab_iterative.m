function [solver_name] = solver_name_matlab_iterative
% SOLVER_NAME_MATLAB_ITERATIVE Defines MATLAB Biconjugate gradients stabilized solver name.
%
% [SOLVER_OUT] = SOLVER_NAME_MATLAB_ITERATIVE  Always returns 'BiCGStab_MATLAB'.
%
% See also SOLVER_NAME_CPU, SOLVER_NAME_GPU, SOLVER_NAME_BACKSLASH, isCUDA. 
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

%% check in/out

narginchk(0,0);
nargoutchk(0,1);

%% BODY

solver_name = 'BiCGStab_MATLAB';

end
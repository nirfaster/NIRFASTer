function [OPTIONS] = solver_options
% SOLVER_OPTIONS Defines BiCGStab solver options.
%
% [OPTIONS] = SOLVER_OPTIONS Always returns structure with default
%   values of BiCGStab solver options. The OPTIONS structure have following
%   entries: 
%       - OPTIONS.no_of_iter (default 1000)
%         Maximum number of BiCGStab iterations.
%       - OPTIONS.tolerance (default 1e-12);
%         Absolute convergence tolerance (norm of solution error).
%       - OPTIONS.rel_tolerance (default 1e-8);
%         Relative convergence tolerance (norm of solution error related to
%         norm of error of first iteration).
%       - OPTIONS.divergence_tol (default 1e8);
%         Relative divergence tolerance (related to first iteration).
%       - OPTIONS.GPU (default -1);
%         Index of GPU requested to be used. By default (-1) the GPU of the
%         highest compute capability is used. Use the 'isCUDA' function to
%         get the installed GPUs info. The indexing is 0-based and the
%         indexing order is the same as in the list of GPUs returned by the
%         'isCUDA' function. Change this value here to switch the used GPU
%         globally.
%
% See also FEMDATA_FD, FEMDATA_TR, FEMDATA_DCS. 
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

%% check in/out

narginchk(0,0);
nargoutchk(0,1);

%% BODY

OPTIONS.no_of_iter = 1000;% (default 1000)
OPTIONS.tolerance = 1e-12;% absolute tolerance (norm of solution error)
% below optins are ignored by MATLAB iterative method
OPTIONS.rel_tolerance = 1e-12;% relative tolerance (norm of solution related to norm of error of first iteration)
OPTIONS.divergence_tol = 1e8;% (default 1e8 for DOUBLE and SINGLE);
% Change the GPU index here here to switch the requested GPU globaly.
OPTIONS.GPU = -1; % by default, do not set the GPU. -1 means choose the best one automatically

end
function [solver_out] = get_solver(solver_in)
% GET_SOLVER Checks and returns supported solver.
% 
% SYNTAX:
%  GET_SOLVER
%  GET_SOLVER(SOLVER_IN)
%  [SOLVER_OUT] = GET_SOLVER
%  [SOLVER_OUT] = GET_SOLVER(SOLVER_IN)
% 
% [SOLVER_OUT] = GET_SOLVER Returns the default solver. If a NVidia GPU
%   with CUDA technology is present, the BiConjugate Gradient Stabilized
%   method parallelized over GPU is used. SOLVER_OUT is set to
%   'BiCGStab_GPU'. Otherwise, the BiConjugate Gradient Stabilized method
%   parallelized over CPU is used. SOLVER_OUT is set to 'BiCGStab_CPU'.
% 
% [SOLVER_OUT] = GET_SOLVER(SOLVER_IN) Check if the user specified solver
%   SOLVER_IN is supported and copies it onto SOLVER_OUT. Issues warning
%   and returns default solver is SOLVER_IN is not supported. Suportes
%   solvers ranked by speed:
%    - 'BiCGStab_GPU' - BiConjugate Gradient Stabilized method. Iterative
%       solver that converges to the 'backslash' solution with a given
%       tolerance. Parallelized for multithreading with NVidia CUDA. Uses
%       all available power of the strongest NVidia GPU installed. This
%       function is limited to work one a single GPU only. Use this option
%       if you have NVidia GPU(s) installed.
%    - 'BiCGStab_CPU' - BiConjugate Gradient Stabilized method. The CPU
%       version is parallelized for multithreading with OpenMP. Uses all
%       available CPU power. Use this option if you do not have NVidia
%       GPU(s).
%    - 'backslash' - Matlab backslash '\'. It will use direct solver with
%       LU decomposition. Most accurate but slow and memory consuming. Use
%       this solver as a gold standard reference when comparing speed,
%       accuracy, etc.. 
%    - 'BiCGStab_MATLAB' - It will use Matlabs implementation of the
%       BiConjugate Gradient Stabilized method. Use this solver for
%       'BiCGStab_GPU' and 'BiCGStab_CPU' references only. Current (R2018b)
%       implementation is slower compared with the 'backslash' and slows
%       further for COMPLEX (FREQUENCY>0) data. Does not handle COMPLEX
%       data well....
%   SOLVER_IN can be also used in short version as:
%    - 'GPU' for 'BiCGStab_GPU'
%    - 'CPU' for 'BiCGStab_CPU'
%    - '\'   for 'backslash'
%   SOLVER_IN names are accessible using following functions:
%    - 'solver_name_GPU'
%    - 'solver_name_CPU'
%    - 'solver_name_matlab_iterative'
%    - 'solver_name_backslash'    
%
% WARNING: MacOS users, no GPU version supported.
% 
% See also isCUDA, SOLVER_NAME_GPU, SOLVER_NAME_CPU,
%          SOLVER_NAME_MATLAB_ITERATIVE, FEMDATA_FD, FEMDATA_TR, 
%          FEMDATA_TR_MOMENTS, FEMDATA_DCS. 
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

%% check in/out

narginchk(0,1);
nargoutchk(0,1);

%% BODY

% set default based on NVidia GPU presence
if isCUDA
    default = solver_name_GPU; %'BiCGStab_GPU';
else
    default = solver_name_CPU; %'BiCGStab_CPU';
end

% if user wants a specified solver
if nargin == 1
    % if GPU requested 
    if strcmp(solver_in, solver_name_GPU) || strcmp(solver_in,'GPU')
        if isCUDA
            % if CUDA GPU present
            solver_out = solver_name_GPU;
        else
            % if no CUDA GPU, step down to parallel CPU
            warning(['''' solver_in ''' requested. No CUDA GPU detected. Using ''' solver_name_CPU ''' instead.'])
            solver_out = solver_name_CPU;
        end
    elseif strcmp(solver_in,solver_name_CPU) || strcmp(solver_in,'CPU')
        % if CPU requested 
        solver_out = solver_name_CPU;
    elseif strcmp(solver_in,solver_name_matlab_iterative)
        % if MATLAB BicGStab requested 
        solver_out = solver_name_matlab_iterative;
    elseif strcmp(solver_in,solver_name_backslash) || strcmp(solver_in,'\')
        % if MATLAB backslash '\' requested
        solver_out = solver_name_backslash;
    else
        % if not supported solver requested
        warning(['Unknown ''' solver_in ''' requested. Using default solver ''' default ''' instead.'])
        solver_out = default;
    end
else
    % use default
    solver_out = default;
end

end


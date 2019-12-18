function [data, info] = get_field_FD_bicgstab_matlab(i_index, j_index, value, qvec, OPTIONS)
% GET_FIELD_FD_BICGSTAB_MATLAB Solves sparse system of linear equations using
%   MATLAB implementation of BiCGStab solver and incomplete LU or Cholesky
%   factorisations for preconditioning.
% 
% SYNTAX:
%   [PHI] = GET_FIELD_FD_BICGSTAB_MATLAB(I_INDEX, J_INDEX, VALUE, QVEC)
%   [PHI] = GET_FIELD_FD_BICGSTAB_MATLAB(I_INDEX, J_INDEX, VALUE, QVEC, OPTIONS)
%   [PHI, INFO] = GET_FIELD_FD_BICGSTAB_MATLAB(___)
% 
%   [PHI] = GET_FIELD_FD_BICGSTAB_MATLAB(I_INDEX, J_INDEX, VALUE, QVEC)
%     Solves a sparse system of linear equations: 
%           A * PHI = QVEC (Ax = b),
%       where A is n by n sparse matrix defined in COO (coordinate list) 
%       format as a list of (row, column, value) tuples (I_INDEX, J_INDEX,
%       VALUE). I_INDEX and J_INDEX should use 0-base indexing.
%       QVEC is the right-hand-side matrix of n by m size.
%       VALUE and QVEC (n by m) can be REAL or COMPLEX.
%       The entries I_INDEX and J_INDEX should be sorted first by row index
%       (I_INDEX) and then by column index (J_INDEX). Index pairs (I_INDEX,
%       J_INDEX) should not repeat. The GEN_MASS_MATRIX_... functions
%       return COO sparse matrix sorted and integrated.
% 
%   [PHI] = GET_FIELD_FD_BICGSTAB_MATLAB(I_INDEX, J_INDEX, VALUE, QVEC, OPTIONS)
%   The OPTIONS is an optional structure of parameters of BiCGStab solver.
%       - OPTIONS.no_of_iter (default 1000)
%         Maximum number of BiCGStab iterations.
%       - OPTIONS.tolerance (default 1e-12);
%         Absolute convergence tolerance (norm of solution error).
%   Entries within the structure are optional.
% 
%   [PHI, INFO] = GET_FIELD_FD_BICGSTAB_MATLAB(___) The INFO holds solver
%     information and has following entries: 
%       - INFO.convergence_iterations (m by 1 vector)
%         Number of iterations to converge to the tolerances for m column
%         of QVEC (right hand side). Rounded up (see 'bicgstab' help to see
%         why we do that).
%       - INFO.convergence_tolerances (m by 1 vector)
%         Exact absolute tolerance at convergence for m column of QVEC
%         (right hand side).
%       - INFO.status (m by 1 vector)
%         Convergence flag (from 0 to 4) as returned by MATLAB 'bicgstab'. For all m column of
%         QVEC (right hand side). '0' means success.
%  
% See also GET_SOLVER, FEMDATA_SPEC_FD, GET_FIELD_FD_CPU, GET_FIELD_FD_CUDA,
%          GEN_MASS_MATRIX_FD_CPU, GEN_MASS_MATRIX_FD_CUDA, GEN_SOURCES. 
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz and H. Dehghani 2018

%% check in/out

narginchk(4,5);
nargoutchk(1,2);

%% solver options

% BiCGStab options
tol = 1e-12; % absolute residuals norm, equivalent of OPTIONS.tolerance above
maxit = 1000; % max number of iterations, equivalent of OPTIONS.no_of_iter above

if nargin == 5
    if isfield(OPTIONS,'rel_tolerance')
        tol = OPTIONS.tolerance;
    else
        warning(['Bad relative tolerance option name. Default value of ' num2str(tol) ' will be used.'])
    end
    if isfield(OPTIONS,'no_of_iter')
        maxit = OPTIONS.no_of_iter;
    else
        warning(['Bad maximum iterations option name. Default value of ' num2str(maxit) ' will be used.'])
    end
end


%% Calculate field for all sources


%% preconditioner
% calculate lower (L) and upper (U) tiangular matrices: L*U = sparse(i_index+1, j_index+1, value)
try
    % The incomplete Cholesky (ichol) is fast but tends to fail for
    % compelex type matrices... 
    % By thefault we use incomplete LU (ilu) decomposition as this is
    % almost as fast as ichol and works for real and complex matrices. 
    [L, U] = ilu(sparse(i_index+1, j_index+1, value));
catch
    % If incomplete LU factorizatin failed, try the full Cholesky factorisation
    % Full LU docomposition worked much slower than full Cholesky on my tests
    % The full Cholesky will be ~100 time slower than the incomlete LU
    L = chol(sparse(i_index+1, j_index+1, value),'lower');
    % in Cholesky: L*L' = sparse(i_index+1, j_index+1, value), thus:
    U = L';
end

%% solver

% loop through the sources and solve
data = zeros(size(qvec)); % status per source
status = zeros(size(qvec,2),1); % status per source
iterations = zeros(size(qvec,2),1); % iterations per source
tolerance = zeros(size(qvec,2),1); % absolute tolerance per source
for ind_source = 1:size(qvec,2)    
    [data(:,ind_source),status(ind_source),~,iterations(ind_source),tolerance_all_steps] = ...
        bicgstab(sparse(i_index+1, j_index+1, value),qvec(:,ind_source),tol,maxit,L,U);
    
    tolerance(ind_source) = tolerance_all_steps(end);
end
if any(status > 0)
    warning(['''' solver_name_matlab_iterative ''' encountered some problems. ''bicgstab'' status code per source: ' ...
        num2str(reshape(status,1,numel(status))) '.'])
end

%% return the info structure if needed

if nargout == 2
    info.convergence_iterations = ceil(iterations);
    info.convergence_tolerances = tolerance;
    info.status = status;
end

end

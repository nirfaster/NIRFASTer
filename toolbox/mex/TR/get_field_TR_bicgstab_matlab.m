function [data, info] = get_field_TR_bicgstab_matlab(i_index, j_index, value_1, value_2, qvec, t, dt, OPTIONS)
% GET_FIELD_TR_BICGSTAB_MATLAB Solves sparse time-varying system of linear
%   equations using MATLAB implementation of BiCGStab solver and incomplete
%   LU or Cholesky factorisations for preconditioning. The time derivative
%   of sparse system utilizes trapezoidal Clark-Nicholson integration
%   method in the time dimension.
% 
% SYNTAX:
%   [PHI] = GET_FIELD_TR_BICGSTAB_MATLAB(I_INDEX, J_INDEX, VALUE_1, VALUE_2, QVEC, T, DT)
%   [PHI] = GET_FIELD_TR_BICGSTAB_MATLAB(I_INDEX, J_INDEX, VALUE_1, VALUE_2, QVEC, T, DT, OPTIONS)
%   [PHI, INFO] = GET_FIELD_TR_BICGSTAB_MATLAB(___)
% 
%   [PHI] = GET_FIELD_TR_BICGSTAB_MATLAB(I_INDEX, J_INDEX, VALUE_1, VALUE_2, QVEC, T, DT)
%     Solves the time (t) derivative of sparse system of linear equations
%     at time steps 0:DT:T-DT:
%         A1 * PHI(t=0) =  QVEC,
%         A1 * PHI(t+DT) = -A2 * PHI(t),
%     where A1 and A2 are n by n sparse matrices in COO (coordinate list) 
%     format as a list of (row, column, value) tuples (I_INDEX, J_INDEX,
%     VALUE_1 (for A1) and VALUE_2 (for A2)). I_INDEX and J_INDEX should
%     use 0-base indexing. The 'gen_mass_matrix_TR...' functions return the
%     indexes and values sorted and integrated.
%     QVEC is the right-hand-side matrix of n by m size.
%     T is the maximum time.
%     DT is the time step.
%     PHI is 3D matrix of n by m by floor(T/DT).
% 
%   [PHI] = GET_FIELD_TR_BICGSTAB_MATLAB(I_INDEX, J_INDEX, VALUE_1, VALUE_2, QVEC, T, DT, OPTIONS)
%     The OPTIONS is an optional structure of parameters of BiCGStab solver.
%       - OPTIONS.no_of_iter (default 1000)
%         Maximum number of BiCGStab iterations.
%       - OPTIONS.tolerance (default 1e-12);
%         Absolute convergence tolerance (norm of solution error).
%   Entries within the structure are optional.
% 
%   [PHI, INFO] = GET_FIELD_TR_BICGSTAB_MATLAB(___) The INFO holds solver
%     information for the first calculation step at t=0 and has following
%     entries:
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
% See also FEMDATA_TR, GET_SOLVER, FEMDATA_SPEC_TR, GET_FIELD_TR_CPU,
%          GET_FIELD_TR_CUDA, GEN_MASS_MATRIX_TR_CPU,
%          GEN_MASS_MATRIX_TR_CUDA, GEN_SOURCES.
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz and H. Dehghani 2018

%% check in/out

narginchk(7,8);
nargoutchk(0,2);

%% solver options

% BiCGStab options
tol = 1e-12; % absolute residuals norm, equivalent of OPTIONS.tolerance above
maxit = 1000; % max number of iterations, equivalent of OPTIONS.no_of_iter above

if nargin == 8
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
% calculate lower (L) and upper (U) tiangular matrices: L*U = sparse(i_index+1, j_index+1, value_1)
try
    % The incomplete Cholesky (ichol) is fast but tends to fail for
    % compelex type matrices... 
    % By thefault we use incomplete LU (ilu) decomposition as this is
    % almost as fast as ichol and works for real and complex matrices. 
    [L, U] = ilu(sparse(i_index+1, j_index+1, value_1));
catch
    % If incomplete LU factorizatin failed, try the full Cholesky factorisation
    % Full LU docomposition worked much slower than full Cholesky on my tests
    % The full Cholesky will be ~100 time slower than the incomlete LU
    L = chol(sparse(i_index+1, j_index+1, value_1),'lower');
    % in Cholesky: L*L' = sparse(i_index+1, j_index+1, value_1), thus:
    U = L';
end

%% solver

N_time = floor(t/dt);

% loop through the sources and solve
data = zeros(size(qvec,1),size(qvec,2),N_time); % nodal data per source per time point
status = zeros(size(qvec,2),1); % status per source
iterations = zeros(size(qvec,2),1); % iterations per source
tolerance = zeros(size(qvec,2),1); % absolute tolerance per source
for ind_time = 1:N_time
    for ind_source = 1:size(qvec,2)
        if ind_time == 1
            % A1 * PHI(t=0) =  QVEC
            [data(:,ind_source,ind_time),status(ind_source),~,iterations(ind_source),tolerance_all_steps] = ...
                bicgstab(sparse(i_index+1, j_index+1, value_1),qvec(:,ind_source),tol,maxit,L,U);
            
            tolerance(ind_source) = tolerance_all_steps(end);
        else
            % A1 * PHI(t+DT) = -A2 * PHI(t)
            [data(:,ind_source,ind_time),~,~,~,~] = ...
                bicgstab(sparse(i_index+1, j_index+1, value_1),-sparse(i_index+1, j_index+1, value_2)*squeeze(data(:,ind_source,ind_time-1)),tol,maxit,L,U);
        end
    end
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

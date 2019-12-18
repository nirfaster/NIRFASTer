% GET_FIELD_TR_CUDA Solves time derivative of sparse system of linear
%   equations using trapezoidal Clark-Nicholson integration method in the
%   time dimension and BiCGStab solver with FSAI preconditioning to solve at
%   derivative steps. Parallelized for multithreading with NVidia CUDA.
%   Uses all available power of the strongest NVidia GPU installed. This
%   function is limited to work one a single GPU only. 
%   The problem is defined as:
%              d PHI(t)
%   A*PHI(t) + -------- = QVEC
%                 dt
%             /
%    PHI(t) = |(QVEC - A*PHI(t))dt
%             /
% SYNTAX:
%   [PHI] = GET_FIELD_TR_CUDA(I_INDEX, J_INDEX, S1_VALUE, S2_VALUE, QVEC, T, DT)
%   [PHI] = GET_FIELD_TR_CUDA(I_INDEX, J_INDEX, S1_VALUE, S2_VALUE, QVEC, T, DT, OPTIONS)
%   [PHI, INFO] = GET_FIELD_TR_CUDA(___)
% 
%   [PHI] = GET_FIELD_TR_CUDA(I_INDEX, J_INDEX, S1_VALUE, S2_VALUE, QVEC, T, DT)
%       Solves the time (t) derivative of sparse system of linear equations: 
%           A1 * PHI(t=0) =  QVEC,
%           A1 * PHI(t+DT) = -A2 * PHI(t),
%       where A1 and A2 are n by n sparse matrices in COO (coordinate list) 
%       format as a list of (row, column, value) tuples (I_INDEX, J_INDEX,
%       S1(2)_VALUE). I_INDEX and J_INDEX should use 0-base indexing.
%       QVEC is the right-hand-side matrix of n by m size.
%       T is the maximum time.
%       DT is the time step.
%       PHI is 3D matrix of n by m by floor(T/DT).
%       S1_VALUE, S2_VALUE and QVEC can be of SINGLE or DOUBLE precision.
%       The input parameters precision determines the BiCGStab solver
%       implementation/calculation precision. Precision affects accuracy
%       and performance. We advise DOUBLE. 
%       The entries I_INDEX and J_INDEX should be sorted first by row index
%       (I_INDEX) and then by column index (J_INDEX). Index pairs (I_INDEX,
%       J_INDEX) should not repeat. The GEN_MASS_MATRIX_TR... functions
%       return COO sparse matrices sorted and integrated.
% 
%   [PHI] = GET_FIELD_TR_CUDA(I_INDEX, J_INDEX, S1_VALUE, S2_VALUE, QVEC, T, DT, OPTIONS)
%   The OPTIONS is an optional structure od parameters of BiCGStab solver: 
%       - OPTIONS.no_of_iter (default 1000)
%         Maximum number of BiCGStab iterations.
%       - OPTIONS.tolerance (default 1e-12 for DOUBLE and 1e-8 for SINGLE);
%         Absolute convergence tolerance.
%       - OPTIONS.rel_tolerance (default 1e-12 for DOUBLE and 1e-8 for SINGLE);
%         Relative convergence tolerance (related to first iteration).
%       - OPTIONS.divergence_tol (default 1e8 for DOUBLE and SINGLE);
%         Relative divergence tolerance (related to first iteration).
%   Entries within the structure are optional.
% 
%   [PHI, INFO] = GET_FIELD_TR_CUDA(___) The
%   INFO holds solver information at t=0 and has following entries:
%       - INFO.convergence_iterations (m by 1 vector)
%         Number of iterations to converge to the tolerances for m column
%         of QVEC (right hand side).
%       - INFO.convergence_tolerances (m by 1 vector)
%         Exact absolute tolerance at convergence for m column of QVEC
%         (right hand side).
%       - INFO.isConverged (m by 1 vector)
%         Logical flag if solution converged to tolerances for m column of
%         QVEC (right hand side).
%       - INFO.isConvergedToAbsTol (m by 1 vector)
%         Logical flag if solution converged to the absolute tolerance for
%         m column of QVEC (right hand side).
%       - INFO.convergence_iterations_limit
%         Maximum number of convergence iterations allowed.
%       - INFO.convergence_tolerance_limit
%         Absolute convergence tolerance requested.
%       - INFO.convergence_relative_tolerance_limit
%         Relative convergence tolerance (related to first iteration)
%         requested. 
%       - INFO. divergence_tolerance_limit
%         Maximum relative divergence tolerance (related to first
%         iteration) allowed. 
% 
% The FSAI preconditioner (approximate inversion of the A1 matrix) is our
% original version which keeps precision at satisfactory level and
% significantly gains speed. The factorised preconditioner matrices are
% calculated for the A1 (I_INDEX, J_INDEX, S1_VALUE) system matrix where only
% up to three values per row are significant. The significant values are
% the ones biggest in absolute value.
% 
% The A1 (I_INDEX,J_INDEX,S1_VALUE) and A2 (I_INDEX, J_INDEX, S2_VALUE)
% matrices should be symmetrical. The solver usually converges for
% non-symmetrical matrices. However, this is not guaranteed as the
% asymmetry of matrices degenerates the FSAI preconditioner and in
% consequence the solver might be unstable and fail.
%
% See also FEMDATA_TR, GEN_MASS_MATRIX_TR_CUDA, GET_FIELD_TR_CPU,
%          GEN_MASS_MATRIX_TR_CPU, GET_FIELD_TR_BICGSTAB_MATLAB.
% 
% Read: Arridge S.R. et al. A finite element approach for modeling photon
%       transport in tissue. Med. Phys., 20, 1993, s. 299–309.
% 
%   MEX File function.
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

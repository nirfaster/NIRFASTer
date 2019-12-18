% GET_FIELD_FD_CUDA Solves sparse system of linear equations using BiCGStab 
%   solver and FSAI preconditioning. Parallelized for multithreading with
%   NVidia CUDA. Uses all available power of the strongest NVidia GPU
%   installed. This function is limited to work one a single GPU only.
% 
% SYNTAX:
%   [PHI] = GET_FIELD_FD_CUDA(I_INDEX, J_INDEX, S_VALUE, QVEC)
%   [PHI] = GET_FIELD_FD_CUDA(I_INDEX,J_INDEX,S_VALUE,QVEC,OPTIONS)
%   [PHI, INFO] = GET_FIELD_FD_CUDA(___)
% 
%   [PHI] = GET_FIELD_FD_CUDA(I_INDEX, J_INDEX, S_VALUE, QVEC) Solves a sparse
%   system of linear equations: 
%           A * PHI = QVEC (Ax = b),
%       where A is n by n sparse matrix defined in COO (coordinate list) 
%       format as a list of (row, column, value) tuples (I_INDEX, J_INDEX,
%       S_VALUE). I_INDEX and J_INDEX should use 0-base indexing.
%       QVEC is the right-hand-side matrix of n by m size.
%       S_VALUE and QVEC (n by m) can be REAL or COMPLEX.
%       S_VALUE and QVEC can be of SINGLE or DOUBLE precision. The input
%       parameters precision determines the BiCGStab solver
%       implementation/calculation precision. Precision affects accuracy
%       and performance. We advise DOUBLE. 
%       The entries I_INDEX and J_INDEX should be sorted first by row index
%       (I_INDEX) and then by column index (J_INDEX). Index pairs (I_INDEX,
%       J_INDEX) should not repeat. The GEN_MASS_MATRIX_... functions
%       return COO sparse matrix sorted and integrated. 
% 
%   [PHI] = GET_FIELD_FD_CUDA(I_INDEX,J_INDEX,S_VALUE,QVEC,OPTIONS) 
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
%   [PHI, INFO] = GET_FIELD_FD_CUDA(___) The
%   INFO holds solver information and has following entries:
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
% The FSAI preconditioner (approximate inversion of the A matrix) is our
% original version which keeps precision at satisfactory level and
% significantly gains speed. The factorised preconditioner matrices are
% calculated for the A (I_INDEX, J_INDEX, S_VALUE) system matrix where only
% up to three values per row are significant. The significant values are
% the ones biggest in absolute value.
% 
% The A (I_INDEX, J_INDEX, S_VALUE) system matrix should be symmetrical. The
% solver usually converges for non-symmetrical A matrices. However, this is
% not guaranteed as the asymmetry of A matrix degenerates the FSAI
% preconditioner and in consequence the solver might be unstable and fail. 
%
% See also GET_FIELD_FD_CPU, GET_FIELD_TR_CPU, GET_FIELD_TR_CUDA,
%          GEN_MASS_MATRIX_FD_CPU, GEN_MASS_MATRIX_FD_CUDA,
%          GEN_MASS_MATRIX_TR_CPU, GEN_MASS_MATRIX_TR_CUDA. 
% 
%   MEX File function.
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

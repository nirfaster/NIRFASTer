% GET_FIELD_TR_MOMENTS_CUDA Solves sparse system of linear equations in a
%   following iterative manner:
%     A * x_{0} = b; => x_{0} = A\b;
%     q*A * x_{q} = B * x_{q-1}; => x_{q} = (qA)\(B * x_{q-1});
%   using BiCGStab solver and FSAI preconditioning. Parallelized for
%   multithreading with NVidia CUDA. Uses all available power of the
%   strongest NVidia GPU installed. This function is limited to work one a
%   single GPU only.
%   This function is used in direct calculations of statistical moments of
%   time-resolved curves.
% 
% SYNTAX:
%   [PHI] = GET_FIELD_TR_MOMENTS_CUDA(I_INDEX, J_INDEX, S1_VALUE, S2_VALUE, QVEC, MAX_ORDER)
%   [PHI] = GET_FIELD_TR_MOMENTS_CUDA(I_INDEX, J_INDEX, S1_VALUE, S2_VALUE, QVEC, MAX_ORDER, OPTIONS)
%   [PHI, INFO] = GET_FIELD_TR_MOMENTS_CUDA(___)
% 
% [PHI] = GET_FIELD_TR_MOMENTS_CUDA(I_INDEX, J_INDEX, S1_VALUE, S2_VALUE, QVEC, MAX_ORDER)
%   Solves the following sparse system of linear equations in an
%   iterative manner: 
%       A1 * PHI_{q} =  QVEC if q=0,
%       q*A1 * PHI_{q} = A2 * PHI_{q-1} if q>0,
%       q = 0, 1, 2, ... , MAX_ORDER-1
%   where A1 and A2 are n by n sparse matrices in COO (coordinate list) 
%   format as a list of (row, column, value) tuples (I_INDEX, J_INDEX,
%   S1(2)_VALUE). I_INDEX and J_INDEX should use 0-base indexing.
%   QVEC is the right-hand-side matrix of n by m size.
%   MAX_ORDER controls number of iterations.
%   PHI is 3D matrix of n by m by MAX_ORDER.
%   S1_VALUE, S2_VALUE and QVEC can be of SINGLE or DOUBLE precision.
%   The input parameters precision determines the BiCGStab solver
%   implementation/calculation precision. Precision affects accuracy
%   and performance. We advise DOUBLE. 
%   The entries I_INDEX and J_INDEX should be sorted first by row index
%   (I_INDEX) and then by column index (J_INDEX). Index pairs (I_INDEX,
%   J_INDEX) should not repeat. The GEN_MASS_MATRIX_TR_MOMENTS... functions
%   return COO sparse matrices sorted and integrated.
%   Default stopping criterai of the BiCGStab are used.
% 
% [PHI] = GET_FIELD_TR_MOMENTS_CUDA(I_INDEX, J_INDEX, S1_VALUE, S2_VALUE, QVEC, MAX_ORDER, OPTIONS)
%   The OPTIONS is an optional structure of parameters of BiCGStab solver: 
%    - OPTIONS.no_of_iter (default 1000)
%      Maximum number of BiCGStab iterations.
%    - OPTIONS.tolerance (default 1e-12 for DOUBLE and 1e-8 for SINGLE);
%      Absolute convergence tolerance.
%    - OPTIONS.rel_tolerance (default 1e-12 for DOUBLE and 1e-8 for SINGLE);
%      Relative convergence tolerance (related to first iteration).
%    - OPTIONS.divergence_tol (default 1e8 for DOUBLE and SINGLE);
%      Relative divergence tolerance (related to first iteration).
%    - OPTIONS.GPU (default -1);
%      Index of GPU requested to be used. By default (-1) the GPU of the
%      highest compute capability is used. Use the 'isCUDA' function to
%      get the installed GPUs info. The indexing is 0-based and the
%      indexing order is the same as in the list of GPUs returned by the
%      'isCUDA' function.
%   Entries within the structure are optional. The 'OPTIONS.tolerance' is
%   automatically squared at following iterations throug orders (up to the
%   MAX_ORDER) as the absolute values of PHI_{q} decrease rapidly (roughly
%   in a square manner).
% 
% [PHI, INFO] = GET_FIELD_TR_MOMENTS_CUDA(___) The
%   INFO holds solver information at the last order (MAX_ORDER) and has
%   following entries: 
%    - INFO.convergence_iterations (m by 1 vector)
%      Number of iterations to converge to the tolerances for m column
%      of QVEC (right hand side).
%    - INFO.convergence_tolerances (m by 1 vector)
%      Exact absolute tolerance at convergence for m column of QVEC
%      (right hand side).
%    - INFO.isConverged (m by 1 vector)
%      Logical flag if solution converged to tolerances for m column of
%      QVEC (right hand side).
%    - INFO.isConvergedToAbsTol (m by 1 vector)
%      Logical flag if solution converged to the absolute tolerance for
%      m column of QVEC (right hand side).
%    - INFO.convergence_iterations_limit
%      Maximum number of convergence iterations allowed.
%    - INFO.convergence_tolerance_limit
%      Absolute convergence tolerance requested.
%    - INFO.convergence_relative_tolerance_limit
%      Relative convergence tolerance (related to first iteration)
%      requested. 
%    - INFO. divergence_tolerance_limit
%      Maximum relative divergence tolerance (related to first
%      iteration) allowed. 
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
% See also FEMDATA_TR_MOMENTS, FEMDATA_STND_TR_MOMENTS,
%          FEMDATA_SPEC_TR_MOMENTS, GEN_MASS_MATRIX_TR_MOMENTS_CUDA,
%          GEN_MASS_MATRIX_TR_MOMENTS_CPU.
% 
% Read: Arridge S.R. et al. A finite element approach for modeling photon
%       transport in tissue. Med. Phys., 20, 1993, s. 299–309.
%       Arridge S.R. and M. Schweiger Photon-measurement density
%       functions. Part 2: Finite-element-method calculations. Applied
%       Optics, Vol. 34, No. 34, 1995, p. 8026-8037
% 
%   MEX File function.
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

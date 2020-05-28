% GEN_MASS_MATRIX_TR_CUDA Builds sparse matrices representation of
%   time-resolved diffusion transport equation within FEM mesh in 2D
%   (triangles) or 3D (tetrahedrons). Parallelized for multithreading with
%   NVidia CUDA. Uses all available power of the strongest NVidia GPU
%   installed. This function is limited to work one a single GPU only. 
% 
% SYNTAX:
%   [I_INDEX, J_INDEX, S1_VALUE, S2_VALUE] = GEN_MASS_MATRIX_TR_CUDA(MESH, DT)
%   [I_INDEX, J_INDEX, S1_VALUE, S2_VALUE] = GEN_MASS_MATRIX_TR_CUDA(MESH, DT, GPU)
%   [I_INDEX, J_INDEX, S1_VALUE, S2_VALUE] = GEN_MASS_MATRIX_TR_CUDA(NODES, ELEMENTS, ...
%                                               BNDVTX, MUA, KAPPA, KSI, CM, DT)
%   [I_INDEX, J_INDEX, S1_VALUE, S2_VALUE] = GEN_MASS_MATRIX_TR_CUDA(NODES, ELEMENTS, ...
%                                               BNDVTX, MUA, KAPPA, KSI, CM, DT, GPU)
% 
%   [I_INDEX, J_INDEX, S1_VALUE, S2_VALUE] = GEN_MASS_MATRIX_TR_CUDA(MESH, DT)
%       Returns the sparse system matrices A1 and A2 in COO (coordinate
%       list) format as a list of (row, column, value) tuples (I_INDEX,
%       J_INDEX, S1_VALUE) and (I_INDEX, J_INDEX, S2_VALUE) sorted and with
%       no duplicates. Sorting is carried out first by row index (I_INDEX)
%       and then by column index (J_INDEX). Row (I_INDEX) and column
%       (J_INDEX) indexes are 0-based. Matrices non-zero entries (S1_VALUE
%       and S2_VALUE) are returned as REAL and DOUBLE.
%       MESH - FEM mesh in NIRFAST structure format. This can be 2D
%        (triangles) or 3D (tetrahedrons) mesh.
%       DT - time integration step in seconds. Can be SINGLE or DOUBLE
%        precision. 
% 
%   [I_INDEX, J_INDEX, S1_VALUE, S2_VALUE] = GEN_MASS_MATRIX_TR_CUDA(MESH, DT, GPU)
%       Works as GEN_MASS_MATRIX_TR_CUDA(MESH, DT). The GPU is an optional
%       index of requested GPU, 0-based and up to the number of capable
%       GPUs installed. Please see also help for the isCUDA function.
% 
%   [I_INDEX, J_INDEX, S1_VALUE, S2_VALUE] = GEN_MASS_MATRIX_TR_CUDA(NODES, ELEMENTS, ...
%                                               BNDVTX, MUA, KAPPA, KSI, CM, DT)
%       Works as GEN_MASS_MATRIX_TR_CUDA(MESH, DT). However, the input MESH
%       structure is represented by separate values:
%       NODES - n by 2 or n by 3 matrix of FEM mesh nodes in Cartesian
%        coordinates system (x,y) or (x,y,z). Each row represents mesh node,
%        columns represent consecutive Cartesian coordinates. The coordinates
%        should be expressed in mm. Can be SINGLE or DOUBLE precision.
%       ELEMENTS - m by 3 or m by 4 nodes connectivity list. List of
%        FEM mesh elements. Each row defines an element showing which NODES
%        are connected together. Entries are positive natural numbers
%        representing row indexes (1-based) of the NODES matrix. Order of
%        ELEMENTS and NODES within ELEMENTS is arbitrary. Sorting is handled
%        internally. Can be INT32, INT64, SINGLE or DOUBLE precision.
%       BNDVTX - n by 1 vector indicating if the MESH node (NODES) belong
%        to the boundary. Can be LOGICAL, INT32, INT64, SINGLE or DOUBLE
%        precision. 
%       MUA - n by 1 vector of absorption at each node (NODES). Expressed in
%        mm^-1. Can be SINGLE or DOUBLE precision.
%       KAPPA - n by 1 vector of diffusion coefficient (KAPPA = [3(MUA(r) +
%        MUSP(r))]^-1) at each node (NODES). Expressed in mm. Can be SINGLE
%        or DOUBLE precision.
%       KSI - n by 1 vector of boundary condition at each node (NODES). KSI
%        = 1./(2*AT) where AT is air boundary attenuation of field (PHI)
%        returning into the medium. AT can be defined as in [1-2]. This
%        parameter does not have an unit. Can be SINGLE or DOUBLE precision.
%       CM - n by 1 vector of speed of light at each node (NODES).
%        Calculated as speed of light in the vacuum normalized by the node
%        refractive index. CM should be expressed in mm/s. Can be SINGLE or
%        DOUBLE precision. 
%       DT - time integration step in seconds. Can be SINGLE or DOUBLE
%        precision. 
% 
%   [I_INDEX, J_INDEX, S1_VALUE, S2_VALUE] = GEN_MASS_MATRIX_TR_CUDA(NODES, ELEMENTS, ...
%                                               BNDVTX, MUA, KAPPA, KSI, CM, DT, GPU)
%       The GPU is an optional index of requested GPU, 0-based and up to
%       the number of capable GPUs installed. Please see also help for the
%       isCUDA function. 
% 
% GEN_MASS_MATRIX_TR_CUDA was tested to give symmetric A1 (I_INDEX, J_INDEX,
% S1_VALUE) and A2 (I_INDEX, J_INDEX, S2_VALUE) matrices. 
%
% We use the following time-resolved diffusion transport equation notation:
%                                                           1       d 
%   div(KAPPA(r)) * grad(PHI(r,t)) - MUA(r) * PHI(r,t) - ------- * ----PHI(r,t) = -S(r,t), (1)
%                                                        c0/N(r)    dt
% where r - FEM mesh node position, KAPPA = [3(MUA(r) + MUSP(r))]^-1 is the
% diffusion coefficient, MUA(r) is absorption coefficient in mm^-1, MUSP(r)
% is the reduced scattering coefficient in mm^-1, PHI(r,t) is te fluence
% rate we are looking for, defined at FEM mesh nodes and at time t
% expressed in seconds, c0 is speed of light in vacuum, N(r) is the medium
% refractive index, j = sqrt(-1) and S(r,t) is a source term defined at all
% nodes and time where sum(S) = 1 and S(r,t~=0) = 0. div - divergence, grad
% - gradient.
% GEN_MASS_MATRIX_TR_CUDA transforms Eq. (1) into its matrix form of:
% 
%   A1 * PHI(r,t=0) =  S(r,t),
%   A1 * PHI(r,t+DT) = -A2 * PHI(t),
% 
% and returns sparse A1 and A2 matrices in COO (coordinate list) format
% sorted and with no duplicates.
% 
% Read: [1] Dehghani H. et al. Near infrared optical tomography using NIRFAST: algorithm
%       for numerical model and image reconstruction. Communications in Numerical
%       Methods in Engineering, 25, 2008, s. 711–732.
%       [2] Haskell R.C. et al. Boundary conditions for the diffusion equation in radiative
%       transfer. J. Opt. Soc. Am. A, 11, 1994, s. 2727–2741.
%       [3] Arridge S.R. et al. A finite element approach for modeling photon
%       transport in tissue. Med. Phys., 20, 1993, s. 299–309.
% 
% See also GEN_MASS_MATRIX_TR_CPU, GET_FIELD_TR_CUDA, GET_FIELD_TR_CPU, ISCUDA.
% 
%   MEX File function.
%   Part of NIRFAST package.
%   S. Wojtkiewicz and H. Dehghani 2018

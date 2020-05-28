% GEN_MASS_MATRIX_TR_MOMENTS_CUDA Builds sparse matrices representation of
%   time-resolved diffusion transport equation for direct calculation of
%   moments of the time-resolved curves. Matrices can represent a FEM mesh
%   in 2D (triangles) or 3D (tetrahedrons) space. Parallelized for
%   multithreading with NVidia CUDA. Uses all available power of the
%   strongest NVidia GPU installed. This function is limited to work one a
%   single GPU only.
% 
% SYNTAX:
%   [I_INDEX, J_INDEX, S1_VALUE, S2_VALUE] = GEN_MASS_MATRIX_TR_MOMENTS_CUDA(MESH)
%   [I_INDEX, J_INDEX, S1_VALUE, S2_VALUE] = GEN_MASS_MATRIX_TR_MOMENTS_CUDA(MESH, GPU)
%   [I_INDEX, J_INDEX, S1_VALUE, S2_VALUE] = GEN_MASS_MATRIX_TR_MOMENTS_CUDA(NODES, ELEMENTS, ...
%                                               BNDVTX, MUA, KAPPA, KSI, CM)
%   [I_INDEX, J_INDEX, S1_VALUE, S2_VALUE] = GEN_MASS_MATRIX_TR_MOMENTS_CUDA(NODES, ELEMENTS, ...
%                                               BNDVTX, MUA, KAPPA, KSI, CM, GPU)
% 
% [I_INDEX, J_INDEX, S1_VALUE, S2_VALUE] = GEN_MASS_MATRIX_TR_MOMENTS_CUDA(MESH)
%   Returns the sparse system matrices A1 and A2 in COO (coordinate
%   list) format as a list of (row, column, value) tuples (I_INDEX,
%   J_INDEX, S1_VALUE) and (I_INDEX, J_INDEX, S2_VALUE) sorted and with
%   no duplicates. Sorting is carried out first by row index (I_INDEX)
%   and then by column index (J_INDEX). Row (I_INDEX) and column
%   (J_INDEX) indexes are 0-based. Matrices non-zero entries (S1_VALUE
%   and S2_VALUE) are returned as REAL and DOUBLE.
%   MESH - FEM mesh in NIRFAST format. This can be 2D (triangles) or 3D
%    (tetrahedrons) mesh. 
% 
% [I_INDEX, J_INDEX, S1_VALUE, S2_VALUE] = GEN_MASS_MATRIX_TR_MOMENTS_CUDA(MESH, GPU)
%   Works as GEN_MASS_MATRIX_TR_MOMENTS_CUDA(MESH). The GPU is an optional
%   index of requested GPU, 0-based and up to the number of capable GPUs
%   installed. Please see also help for the isCUDA function.
% 
% [I_INDEX, J_INDEX, S1_VALUE, S2_VALUE] = GEN_MASS_MATRIX_TR_MOMENTS_CUDA(NODES, ELEMENTS, ...
%                                             BNDVTX, MUA, KAPPA, KSI, CM)
%   Works as the syntax with only MESH as input. However, the input MESH
%   structure is represented by separate values:
%   NODES - n by 2 or n by 3 matrix of FEM mesh nodes in Cartesian
%    coordinates system (x,y) or (x,y,z). Each row represents mesh node,
%    columns represent consecutive Cartesian coordinates. The coordinates
%    should be expressed in mm. Can be SINGLE or DOUBLE precision.
%   ELEMENTS - m by 3 or m by 4 nodes connectivity list. List of
%    FEM mesh elements. Each row defines an element showing which NODES
%    are connected together. Entries are positive natural numbers
%    representing row indexes (1-based) of the NODES matrix. Order of
%    ELEMENTS and NODES within ELEMENTS is arbitrary. Sorting is handled
%    internally. Can be INT32, INT64, SINGLE or DOUBLE precision.
%   BNDVTX - n by 1 vector indicating if the MESH node (NODES) belong
%    to the boundary. Can be LOGICAL, INT32, INT64, SINGLE or DOUBLE
%    precision. 
%   MUA - n by 1 vector of absorption at each node (NODES). Expressed in
%    mm^-1. Can be SINGLE or DOUBLE precision.
%   KAPPA - n by 1 vector of diffusion coefficient (KAPPA = [3(MUA(r) +
%    MUSP(r))]^-1) at each node (NODES). Expressed in mm. Can be SINGLE
%    or DOUBLE precision.
%   KSI - n by 1 vector of boundary condition at each node (NODES). KSI
%    = 1./(2*AT) where AT is air boundary attenuation of field (PHI)
%    returning into the medium. AT can be defined as in [1-2]. This
%    parameter does not have an unit. Can be SINGLE or DOUBLE precision.
%   CM - n by 1 vector of speed of light at each node (NODES).
%    Calculated as speed of light in the vacuum normalized by the node
%    refractive index. CM should be expressed in mm/s. Can be SINGLE or
%    DOUBLE precision.
% 
% [I_INDEX, J_INDEX, S1_VALUE, S2_VALUE] = GEN_MASS_MATRIX_TR_MOMENTS_CUDA(NODES, ELEMENTS, ...
%                                             BNDVTX, MUA, KAPPA, KSI, CM, GPU)
%   The GPU is an optional index of requested GPU, 0-based and up to the
%   number of capable GPUs installed. Please see also help for the isCUDA
%   function.
% 
% GEN_MASS_MATRIX_TR_MOMENTS_CUDA was tested to give symmetric A1 (I_INDEX, J_INDEX,
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
% GEN_MASS_MATRIX_TR_MOMENTS_CUDA transforms Eq. (1) into the following
% matrix form where the M_{q} is a q-th order statistical moment of the
% time-resolved photon fluence rate PHI(r,t):
% 
%   A1 * M_{q}(r) =  S(r,t=0) if q=0,
%   q*A1 * M_{q}(r) = A2 * M_{q-1}(r) if q>0,
%   q = 0, 1, 2, ...
% and the moments definition:
%   M_{0}(r) = int(PHI(r,t),dt,0,inf),
%   M_{q>0}(r) = int(t^q * PHI(r,t),dt,0,inf)/M_{0}(r),
% 
% GEN_MASS_MATRIX_TR_MOMENTS_CUDA returns sparse A1 and A2 matrices in COO
% (coordinate list) format sorted and with no duplicates.
% 
% Read: Arridge S.R. et al. A finite element approach for modeling photon
%       transport in tissue. Med. Phys., 20, 1993, s. 299–309.
%       Arridge S.R. and M. Schweiger Photon-measurement density
%       functions. Part 2: Finite-element-method calculations. Applied
%       Optics, Vol. 34, No. 34, 1995, p. 8026-8037
% 
% See also FEMDATA_TR_MOMENTS, FEMDATA_STND_TR_MOMENTS,
%          FEMDATA_SPEC_TR_MOMENTS, GEN_MASS_MATRIX_TR_MOMENTS_CPU,
%          GET_FIELD_TR_MOMENTS_CUDA, GET_FIELD_TR_MOMENTS_CPU.
% 
%   MEX File function.
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

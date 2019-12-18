% IntgradFgradG_tet4_CPU carries out a spatial convolution of gradients of
%   fields within a FEM mesh. This function runs on parallel on the CPU
%   threads and is parallelized over the data pairs undergoing the
%   convolution.
%
% SYNTAX:
%  [J_BORN] = IntgradFgradG_tet4_CPU(MESH,DATA_1,DATA_2)
%  [J_BORN] = IntgradFgradG_tet4_CPU(MESH,DATA_1,DATA_2,SOURCE_INDEX,DETECTOR_INDEX)
% 
% [J_BORN] = IntgradFgradG_tet4_CPU(MESH,DATA_1,DATA_2) MESH is a 2D or 3D
%   FEM mesh in NIRFAST format used to calculate DATA_1 and DATA_2.
%    DATA_1 - n by m1 matrix (real or complex) of m1 fields (as returned by
%     e.g. 'femdata_FD') defined on n nodes of the MESH. In particular it
%     can be set of photon fluence rate spatial distributions for m1 light
%     sorces witin the MESH. 
%    DATA_2 - n by m2 matrix (real or complex) of m2 fields (as returned by
%     e.g. 'femdata_FD') defined on n nodes of the MESH. In particular it
%     can be set of adjoint photon fluence rate spatial distributions for
%     m2 light sorces positioned at the detectors witin the MESH.
%    J_BORN - n by m1*m2 matrix (real or complex) of all possible
%     combinations of spatial convolutions of gradients of fields DATA_1
%     and DATA_2. Expressed in photon packets per second per millimeter
%     cubed (phot*s^-1*mm^-3) if the input values are in photon fluence
%     rate values (phot*mm^-2s^-1 (3D) and phot*mm^-1s^-1 (2D)).
% 
% [J_BORN] = IntgradFgradG_tet4_CPU(MESH,DATA_1,DATA_2,SOURCE_INDEX,DETECTOR_INDEX)
%   SOURCE_INDEX - p elements vector, permutation vector of DATA_1. It
%    holds indexes of DATA_1 columns that will be convolved. Values should
%    be within 1 and m1 (1-based indexing of DATA_1 columns). E.g. [1 1 2 2] 
%   DETECTOR_INDEX - p elements vector, permutation vector of DATA_2. It
%    holds indexes of DATA_2 columns that will be convolved. Values should
%    be within 1 and m2 (1-based indexing of DATA_2 columns). E.g. [1 2 3 5]
%   J_BORN - n by p matrix (real or complex) of spatial convolutions of
%    gradients of fields DATA_1 and DATA_2 as indexed by SOURCE_INDEX and
%    DETECTOR_INDEX, E.g. convolved pairs gradients of DATA_1 and DATA_2 
%    [1 1;1 2;2 3;2 5].
% 
% HINT
%   #1 BORN approach
%   If you want to use the perturbation theory with the Born approach,
%   please divide J_BORN node-wise by the 'MESH.kappa' (in mm) to get the
%   proper values expressed in DATA_1 and DATA_2 units.
%   #2 RYTOV approach
%   For the Rytov approach please divide J_BORN measurement_pair-wise by
%   the boundary data for DATA_1 (as returned by 'get_boundary_data').
% 
% COMMENT
%   Efficiency tests revealed that parallelization over the convolved pairs
%   is faster than paralleliation over MESH nodes/elements. Tis was also
%   tested on GPU parallel implementation which was actually slower that
%   the CPU parallel implementation.
%   The main challenge here is that the parallelization over mesh elements
%   requires frequent manipulation of nodal values at many threads at the
%   same time. This requires using an atomic operations which slow down the
%   entire algorithm. 
%   We suggest that this effect might be limited by randomizing parallel
%   execution in a way that parallel threads work on random mesh elements
%   at the same time and as such the number of atomic nodal updates sholud
%   be reduced (low probability of parallel processing of connected
%   elements).
%   However, we decided not to test this further as the current
%   CPU-over-pairs parallelization satisfies the current needs. 
% 
% See also IntFG_tet4_CPU, BUILD_JACOBIAN
% 
%   MEX File function.
%   Part of NIRFAST package.
%   S. Wojtkiewicz and H. Dehghani 2018
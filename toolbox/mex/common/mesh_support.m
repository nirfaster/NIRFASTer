% MESH_SUPPORT Calculates total area (2D) or volume (3D) of elements
%   connected to each node. (how much support is around a node).
% 
% SYNTAX:
%   [SUPPORT] = MESH_SUPPORT(NODES,ELEMENTS,ELEMENTS_AREA)
% 
% [SUPPORT] = MESH_SUPPORT(NODES,ELEMENTS,ELEMENTS_AREA) Sums all element
%   areas/volumes around a node.
%   SUPPORT - is n by 1 vector of FEM mesh elements areas (2D) or volumes
%    (3D) summed around n nodes. Expressed in the ELEMENTS_AREA units.
%   NODES - n by d matrix of cartesian coordinates of n FEM mesh nodes. d
%    is the the Euclidean space dimmension: d=2 for 2D (x,y) and d=3 for 3D
%    (x,y,z).
%   ELEMENTS - m by d+1 matrix of connectivity list of mesh elements. Where
%    m is the number of mesh elements and d is the dimension (d=2 trianles,
%    d=3 tetrahedrons). See 'help load_mesh' for more details on ELEMENTS
%    format.
%   ELEMENTS_AREA - is m by 1 vector of m FEM mesh elements area (2D) or
%    volume (3D), expressed in area (2D) or volume (3D) units.
%    ELEMENTS_AREA is returned by balling the 'ele_area_c'.
% 
% See also LOAD_MESH, ELE_AREA_C.
% 
%   MEX File function.
%   Part of NIRFAST package.
%   S. Wojtkiewicz and H. Dehghani 2018
% 
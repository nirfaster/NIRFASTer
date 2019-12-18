% ELE_AREA_C Calculates FEM mesh elements area (2D - triangles) or volume
%   (3D - tetrahedra).
% 
% SYNTAX:
%   [AREA] = ELE_AREA_C(NODES,ELEMENTS)
% 
% [AREA] = ELE_AREA_C(NODES,ELEMENTS) 
%   AREA - is m by 1 vector of m FEM mesh elements area (2D) or volume
%    (3D), expressed in the input units squared (2D) or qubed (3D).
%   NODES - n by d matrix of cartesian coordinates of n FEM mesh nodes. d
%    is the the Euclidean space dimmension: d=2 for 2D (x,y) and d=3 for 3D
%    (x,y,z).
%   ELEMENTS - m by d+1 matrix of connectivity list of mesh elements. Where
%    m is the number of mesh elements and d is the dimension (d=2 trianles,
%    d=3 tetrahedrons). See 'help load_mesh' for more details on ELEMENTS
%    format.
% 
% See also LOAD_MESH, MESH_SUPPORT.
% 
%   MEX File function.
%   Part of NIRFAST package.
%   S. Wojtkiewicz and H. Dehghani 2018
% 
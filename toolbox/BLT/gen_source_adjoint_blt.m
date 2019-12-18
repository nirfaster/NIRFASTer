function [qvec] = gen_source_adjoint_blt(mesh);

% qvec = gen_source_adjoint(mesh);
% Used by jacobian
% caclculates RHS for adjoint source.
% 
% Part of NIRFAST package
% H Dehghani 2006

% Allocate memory
[nnodes,junk]=size(mesh.nodes);
[nmeas,junk]=size(mesh.meas.coord);
qvec = spalloc(nnodes,nmeas,nmeas*5);

% Go through all measurements and integrate to get nodal values
for i = 1 : nmeas
  qvec(mesh.elements(mesh.meas.int_func(i,1),:),i) = ...
      mesh.meas.int_func(i,2:end)' .* ...
      0.01;
end

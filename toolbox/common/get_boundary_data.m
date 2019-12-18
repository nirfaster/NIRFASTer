function [data] = get_boundary_data(mesh,phi)
% GET_BOUNDARY_DATA Calculates FEM mesh boundary data at detector positions
% 
% [DATA] = GET_BOUNDARY_DATA(MESH,PHI) Calculates boundary data based on
%   detector positions 'MESH.meas.coord' and the calculated field PHI. MESH
%   is FEM mesh in NIRFAST format. PHI is n by m matrix of values at n
%   nodes for m sources as enabled in 'MESH.link'. I.e. PHI should not have
%   data for sources disabled in the 'MESH.link'. The returned boundary
%   DATA is a p by 1 vector where p is number of source-detector pairs as
%   specified in 'MESH.link' field. Data values for disabled sources are
%   set to NaN.
% 
% WARNING
%   Only 3rd column of the 'MESH.link' matrix is used to enable/disable
%   surce-detector pairs. Please make sure that this function is called
%   properly, especially in case of 'spec' type meshes, where the columns
%   size of 'MESH.link' is length(MESH.wv)+2. 'MESH.link' turns on/off
%   sorces and wavelengths.
% 
%   Part of NIRFAST package.

%% check in/out

narginchk(2,2);
nargoutchk(0,1);

%% BODY

% issue a warning if there might be a problem with mesh.link
if size(mesh.link,2) > 3
    warning('''mesh.link'' column size exceeds 3. Only the 3rd column of ''mesh.link'' is used to calculate the boundary data.');
end

% check if sources have integration functions
if ~isfield(mesh.meas,'int_func')
    disp('Calculating missing detectors integration functions.');
    % calculate if there is no int functions in the mesh.
    [ind,int_func] = mytsearchn(mesh,mesh.meas.coord);
    mesh.meas.int_func = [ind int_func];
end

% We don't want contributions from internal nodes on boundary values
det_boundary_vartex = logical(mesh.bndvtx(mesh.elements(mesh.meas.int_func(:,1),:)));
% baricentric coordinates for surface facets
% get detector tetrahedron coordinates
int_func = mesh.meas.int_func(:,2:end);
% set not-on-boundary coordinate to 0
int_func(~det_boundary_vartex) = 0;
% rescale the coordinated to get sum equal 1
int_func = int_func./sum(int_func,2);

% define boundary data
data = NaN(size(mesh.link(:,1)));

% indexes of active pairs
index_pairs = find(logical(mesh.link(:,3)));
% [labels of active sources, index of source per active source-detector pair, ~]
[index_num, index_sources, ~] = sources_active_index(mesh);
% [~, index of detector per active source-detector pair, ~]
[~, index_detectors, ~] = detectors_active_index(mesh);
% indexes of nodes of boundary elements with detectors at turned on source-detector pairs
vtx_ind = mesh.elements(mesh.meas.int_func(index_detectors,1),:);

% loop through active source-detectors pairs
for ind_active_pair = 1:length(index_pairs)

    % required as phi columns exist only for active and unique sources! 
%     mask_source = index_num == index_sources(index_pairs(ind_active_pair));
    mask_source = index_num == index_sources(ind_active_pair);

    % Get field values originating at the pair source and registered on nodes of the pair detector element.
    % Multiply the detector element values by the integration functions and add together.
%     data(index_pairs(ind_active_pair)) = int_func(index_detectors(index_pairs(ind_active_pair)),:)*...
%                                          phi(vtx_ind(index_pairs(ind_active_pair),:),mask_source);
    data(index_pairs(ind_active_pair)) = int_func(index_detectors(ind_active_pair),:)*...
                                         phi(vtx_ind(ind_active_pair,:),mask_source);
end
end

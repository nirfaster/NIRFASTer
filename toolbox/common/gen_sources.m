function qvec = gen_sources(mesh)
% GEN_SOURCES Calculates FEM sources vector for those active in mesh.link.
% 
% [QVEC] = GEN_SOURCES(MESH) MESH is NIRFAST mesh structure.
%   QVEC is a sparse, complex matrix of size N by M of initial photon
%   fluence rate at N nodes for M sources. Spatially integrated photon
%   fluence rate for each source equals 1 + 1j*eps. 
% 
% See also GEN_SOURCE_POINT, GEN_SOURCE, FEMDATA_STND_FD. 
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz and H. Dehghani 2018

%% check in/out

narginchk(1,1);
nargoutchk(0,1);

%% BODY
% [labels of active sources, ~, mask of active and unique sources]
[index_num, ~, index_num_mask] = sources_active_index(mesh);

% number of mesh nodes
nnodes = size(mesh.nodes,1);
% number of active and unique sources
nsource = numel(index_num);

% allocate some space for the sources
qvec = spalloc(nnodes,nsource,nsource*100);
% loop through the active sources
for ind_active_source = 1:nsource
    % current source mask
    s_ind = index_num_mask(:,ind_active_source);
    if mesh.source.fwhm(s_ind) == 0
        % if current source is a point
        qvec(:,ind_active_source) = gen_source_point(mesh,s_ind); %#ok<*SPRIX>
    elseif mesh.source.fwhm(s_ind) > 0
        % if current source has Gaussian spatial profile
        qvec(:,ind_active_source) = gen_source(mesh.nodes(:,1:mesh.dimension),...
            sort(mesh.elements')',...
            mesh.dimension,...
            mesh.source.coord(s_ind,1:mesh.dimension),...
            mesh.source.fwhm(s_ind));
    else
        error('The source FWHM value is negative....');
    end
end

%% Catch error in source vector
qvec_mask = sum(qvec,1) ~= 0; % mask of good sources
if sum(qvec_mask) ~= nsource % if not all sources are good
    index_num = reshape(index_num,1,length(index_num)); % prevents from message composition error.
    warning(['Check the quality and accuracy of sources: ' num2str(index_num(~qvec_mask))]);
end
function [index_num, index_link,index_num_mask] = sources_active_index(mesh)
% SOURCES_ACTIVE_INDEX Calculates permutation vectors for mesh sources.
% 
% [INDEX_NUM, INDEX_LINK, INDEX_NUM_MASK] = SOURCES_ACTIVE_INDEX(MESH) MESH
%   is NIRFAST mesh structure. 
%    INDEX_NUM - n by 1 vector where n is the number of active and unique
%      sources as enabled in 'MESH.link'. Holds indexes of active and
%      unique sources as they appear in the 'MESH.source.num'.
%    INDEX_LINK - m by 1 vector where m is the number of active
%      source-detector pairs as enabled in 'MESH.link'. Holds indexes of
%      sources as they appear in the 'MESH.source.num'. E.g. used when
%      calculating boundary data or jacobians for active source-detector
%      pairs.
%    INDEX_NUM_MASK - p by n logical mtrix where p is number of all
%      sources as they appear in the 'MESH.source.num' and n is the number
%      of active and unique sources as enabled in 'MESH.link'. Colums of
%      this matrix can be used for logic indexing of sources coordinates,
%      etc.. E.g. used when calculating source FEM vectors. 
% 
% See also DETECTORS_ACTIVE_INDEX, GEN_SOURCES, GET_BOUNDARY_DATA.
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz

%% check in/out

narginchk(1,1);
nargoutchk(0,3);

%% BODY

% mask of active sources (sources - rows, source-detector pairs - columns)
active_mesh = (mesh.source.num == mesh.link(:,1)') & logical(mesh.link(:,3))';
% indexes of active and unique sources 
index_num = find(sum(active_mesh,2)>0);

% indexes of active sources as they appear in the 'MESH.source.num' per
% active source-detector pair
[index_link,~] = find(active_mesh);

% the logical mask of the unique and active sources 
index_num_mask = mesh.source.num == index_num';

end
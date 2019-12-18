function [index_num, index_link, index_num_mask] = detectors_active_index(mesh)
% DETECTORS_ACTIVE_INDEX Calculates permutation vectors for mesh detectors.
% 
% [INDEX_NUM, INDEX_LINK, INDEX_NUM_MASK] = DETECTORS_ACTIVE_INDEX(MESH) MESH
%   is NIRFAST mesh structure. 
%    INDEX_NUM - n by 1 vector where n is the number of active and unique
%      detectors as enabled in 'MESH.link'. Holds indexes of active and
%      unique detectors as they appear in the 'MESH.meas.num'.
%    INDEX_LINK - m by 1 vector where m is the number of active
%      source-detector pairs as enabled in 'MESH.link'. Holds indexes of
%      detectors as they appear in the 'MESH.meas.num'. E.g. used when
%      calculating boundary data or jacobians for active source-detector
%      pairs.
%    INDEX_NUM_MASK - p by n logical mtrix where p is number of all
%      detectors as they appear in the 'MESH.meas.num' and n is the number
%      of active and unique detectors as enabled in 'MESH.link'. 
% 
% See also SOURCES_ACTIVE_INDEX, GET_BOUNDARY_DATA.
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz

%% check in/out

narginchk(1,1);
nargoutchk(0,3);

%% BODY

% mask of active detectors (detectors - rows, source-detector pairs - columns)
active_mesh = (mesh.meas.num == mesh.link(:,2)') & logical(mesh.link(:,3))';
% indexes of active and unique detectors 
index_num = find(sum(active_mesh,2)>0);

% indexes of active detectors as they appear in the 'MESH.meas.num' per
% active source-detector pair
[index_link,~] = find(active_mesh);

% the logical mask of the unique and active detectors 
index_num_mask = mesh.meas.num == index_num';

end
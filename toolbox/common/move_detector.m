function mesh = move_detector(mesh)
% MOVE_DETECTOR Moves MESH detectors onto the surface of the MESH
% 
%   [MESH] = MOVE_DETECTOR(MESH) Where MESH is a FEM mesh in the NIRFAST
%    format.
%  
% See also MYTSEARCHN, MYTSEARCHN_FAST, LOAD_MESH, MOVE_SOURCE.
% 
%   Part of NIRFAST package.
%   H. Dehghani and S. Wojtkiewicz 2018

%% check in/out

narginchk(1,1);
nargoutchk(1,1);


%% load mesh if needed

if ischar(mesh) || isstring(mesh)
  mesh = load_mesh(mesh);
end

%% Check if mesh has detectors to move

if ~isfield(mesh,'meas') || ~isfield(mesh.meas,'coord')
    error('No detectors present');
end

%% get list of boundary faces

% if 3D mesh
if size(mesh.elements,2) == 4
    % logical mask of elements that touch the boundary
    mask_touch_surface = sum(mesh.bndvtx(mesh.elements),2) > 0;
    % all possible faces that touch the boundary
    faces = [mesh.elements(mask_touch_surface,[1,2,3]);
              mesh.elements(mask_touch_surface,[1,2,4]);
              mesh.elements(mask_touch_surface,[1,3,4]);
              mesh.elements(mask_touch_surface,[2,3,4])];
    % sort vertex indexes to make it comparable
    faces = sort(faces,2);
    % take unique faces
    faces = unique(faces,'rows');
    % take faces where all three vertices are on the boundary
    faces = faces(sum(mesh.bndvtx(faces),2)==mesh.dimension,:);

% if 2D mesh or 3D surface mesh
elseif size(mesh.elements,2) == 3
    if mesh.dimension == 3
        % if 3D surface mesh
        % take faces where all three vertices are on the boundary
        faces = mesh.elements(sum(mesh.bndvtx(mesh.elements),2)==mesh.dimension,:);
    elseif mesh.dimension == 2
        % logical mask of elements that touch the boundary
        mask_touch_surface = sum(mesh.bndvtx(mesh.elements),2) > 0;
        faces = [mesh.elements(mask_touch_surface,[1,2]);
                  mesh.elements(mask_touch_surface,[1,3]);
                  mesh.elements(mask_touch_surface,[2,3])];
        % sort vertex indexes to make it comparable
        faces = sort(faces,2);
        % take unique edges
        faces = unique(faces,'rows');
        % take edges where all two vertices are on the boundary
        faces = faces(sum(mesh.bndvtx(faces),2)==mesh.dimension,:);
    end
end

%% loop through detectors

for i=1:size(mesh.meas.coord,1)
    
    if mesh.dimension == 2
        
        % find closest boundary node
        dist = distance(mesh.nodes,mesh.bndvtx,[mesh.meas.coord(i,1:mesh.dimension) 0]);
        [~, r0_ind] = min(dist);
        
        % find faces including the closest boundary node
        fi = faces(sum(faces==r0_ind,2)>0,:);

        % find closest face
        dist = zeros(size(fi,1),1);
        point = zeros(size(fi,1),mesh.dimension);
        for j = 1:size(fi,1)
            [dist(j),point(j,:)] = pointLineDistance(mesh.nodes(fi(j,1),:), ...
                mesh.nodes(fi(j,2),:),mesh.meas.coord(i,1:mesh.dimension));
        end
        [~, smallest] = min(dist);
        
        % move detector to the closest point on that edge
        mesh.meas.coord(i,1:mesh.dimension) = point(smallest,1:mesh.dimension);
        
    elseif mesh.dimension == 3

        % find closest boundary node
        dist = distance(mesh.nodes,mesh.bndvtx,mesh.meas.coord(i,1:mesh.dimension));
        
        [~, r0_ind] = min(dist);

        % find faces including the closest boundary node
        fi = faces(sum(faces==r0_ind,2)>0,:);

        % find closest face
        dist = zeros(size(fi,1),1);
        point = zeros(size(fi,1),mesh.dimension);
        for j = 1:size(fi,1)
            [dist(j),point(j,:)] = pointTriangleDistance([mesh.nodes(fi(j,1),:);...
                mesh.nodes(fi(j,2),:);mesh.nodes(fi(j,3),:)],mesh.meas.coord(i,1:mesh.dimension));
        end
        [~, smallest] = min(dist);
        
        % move detector to the closest point on that face
        mesh.meas.coord(i,1:mesh.dimension) = point(smallest,1:mesh.dimension);
    end
end

end

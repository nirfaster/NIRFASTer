function mesh = move_source(mesh,mus_eff)
% MOVE_SOURCE Moves MESH sources inside the MESH by 1 scattering distance.
% 
%   [MESH] = MOVE_SOURCE(MESH, MUS_EFF) Moves MESH sources
%    (MESH.source.coord) inside the MESH by 1 scattering distance
%    (1/MUS_EFF).
%  
% See also MYTSEARCHN, MYTSEARCHN_FAST, LOAD_MESH, MOVE_DETECTOR. 
% 
%   Part of NIRFAST package.
%   H. Dehghani and S. Wojtkiewicz 2018

%% check in/out

narginchk(2,2);
nargoutchk(1,1);

% set succes flag
isFailed = false;

%% load mesh if needed
if ischar(mesh) || isstring(mesh)
  mesh = load_mesh(mesh);
end

%% Check if mesh has sources to move

if ~isfield(mesh,'source') || ~isfield(mesh.source,'coord')
    error('No sources present');
end

%% check if mus_eff is unrealistic for the mesh size

scatt_dist = 1/mean(mus_eff);
display(['mus_eff: ', num2str(mean(mus_eff))])
display(['scattering_distance: ', num2str(scatt_dist)])
mesh_size = max(mesh.nodes,[],1)-min(mesh.nodes,[],1);
if scatt_dist*10 > min(mesh_size)
    warning(['Mesh is too small for the scattering coefficient given. '...
        'Minimal mesh size: ' num2str(min(mesh_size)) 'mm. Scattering distance: ' num2str(scatt_dist) 'mm. '...
        num2str(min(mesh_size)/10) ' mm will be used for scattering distance. '...
        'You might want to ensure that the scale of your mesh and the scattering coefficient are in mm.']);
    
    scatt_dist = min(mesh_size)/10;
end

%% get list of boundary faces

out_normal = 0;
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
        out_normal = 1;
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

%% loop through sources

for i = 1:size(mesh.source.coord,1)
    
    if mesh.dimension == 2
        
        % find closest boundary node
        dist = distance(mesh.nodes,mesh.bndvtx,[mesh.source.coord(i,1:mesh.dimension) 0]);
        [~, r0_ind] = min(dist);
        
        % find edges including the closest boundary node
        fi = faces(sum(faces==r0_ind,2)>0,:);

        % find closest edge
        dist = zeros(size(fi,1),1);
        point = zeros(size(fi,1),mesh.dimension);
        for j = 1:size(fi,1)
            [dist(j),point(j,:)] = pointLineDistance(mesh.nodes(fi(j,1),:), ...
                mesh.nodes(fi(j,2),:),mesh.source.coord(i,1:mesh.dimension));
        end
        [~, smallest] = min(dist);
        
        % find normal of that edge
        a = mesh.nodes(fi(smallest(1),1),:);
        b = mesh.nodes(fi(smallest(1),2),:);
        n = [a(2)-b(2) b(1)-a(1)];
        n = n/norm(n);
        
        % move source inside mesh by 1 scattering distance
        pos1 = point(smallest,1:mesh.dimension) + n*scatt_dist;
        pos2 = point(smallest,1:mesh.dimension) - n*scatt_dist;
        [ind,int_func] = mytsearchn(mesh,[pos1;pos2]);
        if ~isnan(ind(1))
            mesh.source.coord(i,:) = pos1;
            mesh.source.int_func(i,:) = [ind(1) int_func(1,:)];
        elseif ~isnan(ind(2))
            mesh.source.coord(i,:) = pos2;
            mesh.source.int_func(i,:) = [ind(2) int_func(2,:)];
        else
            isFailed = true;
        end
        
    elseif mesh.dimension == 3

        % find closest boundary node (r0_ind)
        dist = distance(mesh.nodes,mesh.bndvtx,mesh.source.coord(i,1:mesh.dimension));
        [~, r0_ind] = min(dist);
        
        % take faces including the closest boundary node
        fi = faces(sum(faces==r0_ind,2)>0,:);

        % find closest face
        dist = zeros(size(fi,1),1);
        point = zeros(size(fi,1),mesh.dimension);
        for j = 1:size(fi,1)
            [dist(j),point(j,:)] = pointTriangleDistance([mesh.nodes(fi(j,1),:);...
                mesh.nodes(fi(j,2),:);mesh.nodes(fi(j,3),:)],mesh.source.coord(i,1:mesh.dimension));
        end
        [~, smallest] = min(dist);
        
        % find normal of that face
        a = mesh.nodes(fi(smallest,1),:);
        b = mesh.nodes(fi(smallest,2),:);
        c = mesh.nodes(fi(smallest,3),:);
        n = cross(b-a,c-a);
        n = n/norm(n);
        
        % move source inside mesh by 1 scattering distance
        pos2 = point(smallest(1),:) + n*scatt_dist;
        pos1 = point(smallest(1),:) - n*scatt_dist;
        
        if ~out_normal
            % if not 3D surface mesh
            [ind,int_func] = mytsearchn(mesh,[pos1; pos2]);
        end
        if out_normal
            % if 3D surface mesh (normal outside)
            mesh.source.coord(i,:) = pos1;
            [ind,int_func] = mytsearchn(mesh,pos1);
            mesh.source.int_func(i,:) = [ind(1) int_func(1,:)];
        elseif ~isnan(ind(1))
            % first point OK
            mesh.source.coord(i,:) = pos1;
            mesh.source.int_func(i,:) = [ind(1) int_func(1,:)];
        elseif ~isnan(ind(2))
            % second point OK
            mesh.source.coord(i,:) = pos2;
            mesh.source.int_func(i,:) = [ind(2) int_func(2,:)];
        else
            isFailed = true;
        end
    end
        
end

if isFailed
    warning('Source(s) could not be moved. The mesh structure may be poor.');
end

end

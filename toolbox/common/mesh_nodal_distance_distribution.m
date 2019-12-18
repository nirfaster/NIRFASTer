function [varargout] = mesh_nodal_distance_distribution(mesh,varargin)
% MESH_NODAL_DISTANCE_DISTRIBUTION Plots FEM mesh edges length
%   distribution.
% 
% SYNTAX:
%  MESH_NODAL_DISTANCE_DISTRIBUTION(MESH)
%  MESH_NODAL_DISTANCE_DISTRIBUTION(MESH,RESOLUTION)
%  [D] = MESH_NODAL_DISTANCE_DISTRIBUTION(MESH)
%  [D] = MESH_NODAL_DISTANCE_DISTRIBUTION(MESH,RESOLUTION)
% 
%  [D] = MESH_NODAL_DISTANCE_DISTRIBUTION(MESH,RESOLUTION) MESH is the mesh in
%    NIRFAST format. MESH can be in 2D or 3D.
%    - RESOLUTION - is a scalar to set number of edges length bins in
%                   the distribution. By default RESOLUTION=128.
%    - D - (optional) is m by 6 matrix of edges length of m tetrahedrons
%          for 3D MESH or m by 3 matrix of edges length of m triangles for
%          2D MESH. Expressed in the MESH coordionates unit.
% 
% See also mesh_volume_distribution, mesh_facets_area_distribution,
%          mesh_heights_distribution.  
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

%% check in/out

narginchk(1,2);
nargoutchk(0,1);

%% get optional inputs, handle variable input

% default
resolution = 128;

if ~isempty(varargin)
    if length(varargin) == 1
        % frequency
        if ~ischar(varargin{1}) && ~isstring(varargin{1}) && (numel(varargin{1})==1)
            % sanity check
            if varargin{1} < 0
                error('Negative bins number. Please see help for details on how to use this function.')
            else
                resolution = varargin{1};
            end
        else
            error('Bad 2nd argument value. A scalar expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) > 1
        error('Bad arguments. Please see the help for details on how to use this function.')
    end
end

%% load if needed

if ~isstruct(mesh)
    mesh = load_mesh(mesh);
end

%% BODY

a = mesh.nodes(mesh.elements(:,1),:); % vertices 1
b = mesh.nodes(mesh.elements(:,2),:); % vertices 2
c = mesh.nodes(mesh.elements(:,3),:); % vertices 3
if mesh.dimension == 3
    d = mesh.nodes(mesh.elements(:,4),:); % vertices 4
    D = nodal_distance_tetrahedron(a,b,c,d);
else
    D = nodal_distance_triangle(a,b,c);
end

% set output if needed
if nargout == 1
    varargout{1} = D;
end

%% get basic statistics
D = reshape(D,numel(D),1);

dist_mean = mean(D);
dist_std = std(D);

%% plot the histogram

figure()
histogram(D,resolution);
xlabel('nodal distance /mm','Interpreter','Latex','FontSize',18)
ylabel('number of edges /-','Interpreter','Latex','FontSize',18)
legend({['$<$distance$>=' num2str(dist_mean,'%4.2f') '\pm' num2str(dist_std,'%4.2f') '$ /mm']},'Interpreter','Latex','FontSize',18)

end
function [varargout] = mesh_volume_distribution(mesh, varargin)
% MESH_VOLUME_DISTRIBUTION Plots FEM mesh elements volume(3D)/surface(2D)
%   distribution. 
% 
% SYNTAX:
%  MESH_VOLUME_DISTRIBUTION(MESH)
%  MESH_VOLUME_DISTRIBUTION(MESH,RESOLUTION)
%  [V] = MESH_VOLUME_DISTRIBUTION(MESH)
%  [V] = MESH_VOLUME_DISTRIBUTION(MESH,RESOLUTION)
% 
%  [V] = MESH_VOLUME_DISTRIBUTION(MESH,RESOLUTION) MESH is the mesh in
%    NIRFAST format. MESH can be in 2D or 3D.
%    - RESOLUTION - is a scalar to set number of volume bins in
%                   the distribution. By default RESOLUTION=128.
%    - V - (optional) is m elements vector of m tetrahedrons volume (3D
%          mesh) or m triangles surface (2D mesh). Expressed in the MESH
%          coordionates unit. 
% 
% See also mesh_nodal_distance_distribution, mesh_facets_area_distribution,
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

%% calculate volume (3D) surface (2D)

a = mesh.nodes(mesh.elements(:,1),:); % vertices 1
b = mesh.nodes(mesh.elements(:,2),:); % vertices 2
c = mesh.nodes(mesh.elements(:,3),:); % vertices 3
if mesh.dimension == 3
    d = mesh.nodes(mesh.elements(:,4),:); % vertices 4
    % elements volume
    V = volume_tetrahedron(a,b,c,d);
else
    % elements area
    V = surface_triangle(a,b,c);
end
% set output if needed
if nargout == 1
    varargout{1} = V;
end

%% get basic statistics
vol_mean = mean(V);
vol_std = std(V);

%% plot the histogram

figure()
histogram(V,resolution);
if mesh.dimension == 3
    xlabel('element volume /mm$^{3}$','Interpreter','Latex','FontSize',18)
    legend({['$<$volume$>=' num2str(vol_mean,'%4.2f') '\pm' num2str(vol_std,'%4.2f') '$ /mm$^{3}$']},'Interpreter','Latex','FontSize',18)
else
    xlabel('element surface /mm$^{2}$','Interpreter','Latex','FontSize',18)
    legend({['$<$surface$>=' num2str(vol_mean,'%4.2f') '\pm' num2str(vol_std,'%4.2f') '$ /mm$^{2}$']},'Interpreter','Latex','FontSize',18)
end
ylabel('number of elements /-','Interpreter','Latex','FontSize',18)


end
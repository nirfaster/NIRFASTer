function [varargout] = mesh_facets_area_distribution(mesh,varargin)
% MESH_FACETS_AREA_DISTRIBUTION Plots FEM mesh facets area distribution.
% 
% SYNTAX:
%  MESH_FACETS_AREA_DISTRIBUTION(MESH)
%  MESH_FACETS_AREA_DISTRIBUTION(MESH,RESOLUTION)
%  [S] = MESH_FACETS_AREA_DISTRIBUTION(MESH)
%  [S] = MESH_FACETS_AREA_DISTRIBUTION(MESH,RESOLUTION)
% 
%  [S] = MESH_FACETS_AREA_DISTRIBUTION(MESH,RESOLUTION) MESH is the mesh in
%    NIRFAST format. MESH can be in 2D or 3D.
%    - RESOLUTION - is a scalar to set number of area bins in
%                   the distribution. By default RESOLUTION=128.
%    - S - is m by 4 matrix of m facets area of 3D MESH tetrahedrons or m
%          elements vector of surfaces of 2D MESH triangles. Expressed in
%          the MESH coordionates unit.
% 
% See also mesh_nodal_distance_distribution, mesh_volume_distribution,
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
    % facets area
    S = facets_area_tetrahedron(a,b,c,d);
else
    % elements area
    S = surface_triangle(a,b,c);
end
% set output if needed
if nargout == 1
    varargout{1} = S;
end

%% get basic statistics
S = reshape(S,numel(S),1);

area_mean = mean(S);
area_std = std(S);

%% plot the histogram

figure('Name',mesh.name)
histogram(S,resolution);
xlabel('area /mm$^{2}$','Interpreter','Latex','FontSize',18)
if mesh.dimension == 3
    ylabel('number of facets /-','Interpreter','Latex','FontSize',18)
else
    ylabel('number of elements /-','Interpreter','Latex','FontSize',18)
end
legend({['$<$area$>=' num2str(area_mean,'%4.2f') '\pm' num2str(area_std,'%4.2f') '$ /mm$^{2}$']},'Interpreter','Latex','FontSize',18)

end

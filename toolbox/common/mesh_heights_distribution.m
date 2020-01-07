function [varargout] = mesh_heights_distribution(mesh,varargin)
% MESH_HEIGHTS_DISTRIBUTION Plots FEM mesh elements heights distribution. 
% 
% SYNTAX:
%  MESH_HEIGHTS_DISTRIBUTION(MESH)
%  MESH_HEIGHTS_DISTRIBUTION(MESH,RESOLUTION)
%  [H] = MESH_HEIGHTS_DISTRIBUTION(MESH)
%  [H] = MESH_HEIGHTS_DISTRIBUTION(MESH,RESOLUTION)
% 
%  [H] = MESH_HEIGHTS_DISTRIBUTION(MESH,RESOLUTION) MESH is the mesh in
%    NIRFAST format. MESH can be in 2D or 3D.
%    - RESOLUTION - is a scalar to set number of length bins in
%                   the distribution. By default RESOLUTION=128.
%   - H - is m by 4 matrix of m tetrahedrons heights for 3D MESH or m by 3
%         matrix of m triangles heights for 2D MESH. Expressed in the MESH
%         coordionates unit.   
% 
% See also mesh_nodal_distance_distribution, mesh_facets_area_distribution,
%          mesh_volume_distribution.  
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
    % tetrahedron heights
    h = height_tetrahedron(a,b,c,d);
else
    % triangle heights
    h = height_triangle(a,b,c);
end

% set output if needed
if nargout == 1
    varargout{1} = h;
end

%% get basic statistics
h = reshape(h,numel(h),1);

height_mean = mean(h);
height_std = std(h);

%% plot the histogram

figure('Name',mesh.name)
histogram(h,resolution);
xlabel('elements height /mm','Interpreter','Latex','FontSize',18)
ylabel('number of heights /-','Interpreter','Latex','FontSize',18)
legend({['$<$height$>=' num2str(height_mean,'%4.2f') '\pm' num2str(height_std,'%4.2f') '$ /mm']},'Interpreter','Latex','FontSize',18)

end
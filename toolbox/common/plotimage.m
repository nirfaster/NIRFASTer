function h = plotimage(mesh,val)
% PLOTIMAGE Fast prewiev of data within FEM mesh in the NIRFAST format.
%   Plots an image of the values on the mesh.
%
% SYNTAX:
%   H = PLOTIMAGE(MESH,VAL)
%
% H = PLOTIMAGE(MESH,VAL) MESH is the input mesh in NIRFAST format
% (workspace variable), VAL is n by 1 vector of values to plot, where n is
% the number of MESH nodes. H is the figure handle.
% 
%   Part of NIRFAST package.

%% check input
if size(val,2)==size(mesh.nodes,1)
    val = val';
end

if (size(val,1) ~= size(mesh.nodes,1)) || (size(val,2) ~= 1)
    warning('Check size of parameter to be plotted');
    h = 0;
    return
end

%% plot
figure('Name',mesh.name);
set(gcf,'Color',[0.95 0.95 0.95]);

% h = trisurf(mesh.elements,...
% 	    mesh.nodes(:,1),...
% 	    mesh.nodes(:,2),...
% 	    mesh.nodes(:,3),...
% 	    val);

h = patch('Faces',mesh.elements,'Vertices',mesh.nodes,...
    'EdgeColor','none','FaceColor','interp',...
    'FaceVertexCData',val);

shading interp;

if mesh.dimension == 2
    view(2);
else
    view(3);
end

colorbar('horiz');

axis equal; 
axis off;
colormap hot;

end
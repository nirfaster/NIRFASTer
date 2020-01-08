function h = plotimage(mesh,val)
% h = plotimage(mesh,val)
%
% Plots an image of the values on the mesh
% 
% val is the values
% mesh is the input mesh (workspace variable)
% h is the plot

if size(val,2)==size(mesh.nodes,1)
    val = val';
end

if (size(val,1) ~= size(mesh.nodes,1)) || (size(val,2) ~= 1)
    warning('Check size of parameter to be plotted');
    h = 0;
    return
end

figure('Name',mesh.name);
set(gcf,'Color',[0.95 0.95 0.95]);

% set(gca,'FontSize',28);
% h = trisurf(mesh.elements,...
% 	    mesh.nodes(:,1),...
% 	    mesh.nodes(:,2),...
% 	    mesh.nodes(:,3),...
% 	    val);
h=patch('Faces',mesh.elements,'Vertices',mesh.nodes,...
    'EdgeColor','none','FaceColor','interp','FaceVertexCData',...
    val);
shading interp;
if mesh.dimension == 2
    view(2);
end
colorbar('horiz');
axis equal; 
axis off;
colormap hot;

end
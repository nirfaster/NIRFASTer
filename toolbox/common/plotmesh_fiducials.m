function plotmesh_fiducials(mesh)
% PLOTMESH_FIDUCIALS Plots sources and detectors of meshes in NIRFAST format.
% 
% SYNTAX:
%   PLOTMESH_FIDUCIALS(MESH)
% 
%   PLOTMESH(MESH) MESH is NIRFAST mesh structure.
%
%   Part of NIRFAST package.

%% load mesh
if ischar(mesh)== 1
    mesh = load_mesh(mesh);
end

%% plot sources/detectors

figure('Name',mesh.name)

hold on
% ind = find(mesh.bndvtx==1);
if mesh.dimension == 2    
    
    %     plot(mesh.nodes(ind,1),mesh.nodes(ind,2),'c.');

    TR = triangulation(mesh.elements, mesh.nodes(:,1:mesh.dimension));
    triplot(TR,'Color','cyan')
    
    F = freeBoundary(TR);
    plot(mesh.nodes(F,1),mesh.nodes(F,2),'-', 'Color','cyan','LineWidth',2)
    
    if isfield(mesh,'source') == 1
        tmp = sort(mesh.source.num);
        if length(tmp)>1
            s1 = mesh.source.coord(mesh.source.num == tmp(1),:);
            s2 = mesh.source.coord(mesh.source.num == tmp(2),:);
            p_sour = plot(s1(:,1),s1(:,2),'go',s2(:,1),s2(:,2),'yo',...
                mesh.source.coord(3:end,1),...
                mesh.source.coord(3:end,2),'ro','LineWidth',2,'MarkerSize',8);
        else
            s1 = mesh.source.coord;
            p_sour = plot(s1(:,1),s1(:,2),'go','LineWidth',2,'MarkerSize',8);
        end 
    end
    if isfield(mesh,'meas') == 1
        p_det = plot(mesh.meas.coord(:,1),...
            mesh.meas.coord(:,2),'bx','LineWidth',2,'MarkerSize',8);
    end
    
    axis equal;
elseif mesh.dimension == 3
    
%     plot3(mesh.nodes(ind,1),...
%         mesh.nodes(ind,2),...
%         mesh.nodes(ind,3),'c.');

    TR = triangulation(mesh.elements, mesh.nodes(:,1:mesh.dimension));
    [F,P] = freeBoundary(TR);
    trisurf(F,P(:,1),P(:,2),P(:,3), 'FaceColor','cyan','FaceAlpha',0.05);
    
    if isfield(mesh,'source') == 1
        tmp = sort(mesh.source.num);
        if length(tmp)>1
            s1 = mesh.source.coord(mesh.source.num == tmp(1),:);
            s2 = mesh.source.coord(mesh.source.num == tmp(2),:);
            p_sour = plot3(s1(:,1),s1(:,2),s1(:,3),'go',s2(:,1),s2(:,2),s2(:,3),'yo',...
                mesh.source.coord(3:end,1),...
                mesh.source.coord(3:end,2),...
                mesh.source.coord(3:end,3),'ro','LineWidth',2,'MarkerSize',8);
        else
            s1 = mesh.source.coord;
            p_sour = plot3(s1(:,1),s1(:,2),s1(:,3),'go','LineWidth',2,'MarkerSize',8);
        end

    end
    if isfield(mesh,'meas') == 1
        p_det = plot3(mesh.meas.coord(:,1),...
            mesh.meas.coord(:,2),...
            mesh.meas.coord(:,3),'bx',...
            'LineWidth',2,'MarkerSize',8);
    end
    
    axis equal;
end

if isfield(mesh,'source') && size(mesh.source.coord,1)>2
    if isfield(mesh,'meas')
        legend([p_sour(1:3)' p_det(1)], 'Source 1','Source 2','Sources +','Detector');
    else
        legend(p_sour(1:3)', 'Source 1','Source 2','Sources +');
    end
elseif isfield(mesh,'source') && size(mesh.source.coord,1)==2
    if isfield(mesh,'meas')
        legend([p_sour(1:2)' p_det(1)], 'Source 1','Source 2','Detector');
    else
        legend(p_sour(1:3), 'Source 1','Source 2');
    end
elseif isfield(mesh,'source') && size(mesh.source.coord,1)==1
    if isfield(mesh,'meas')
        legend([p_sour(1) p_det(1)], 'Source 1','Detector');
    else
        legend(p_sour(1), 'Source 1');
    end
end

xlabel('x position (mm)');
ylabel('y position (mm)');
zlabel('z position (mm)');

hold off

end
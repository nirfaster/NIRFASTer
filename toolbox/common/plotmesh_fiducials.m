function plotmesh_fiducials(mesh)

%% load mesh
if ischar(mesh)== 1
    mesh = load_mesh(mesh);
end

%% plot sources/detectors

figure()


hold on;
ind = find(mesh.bndvtx==1);
if mesh.dimension == 2
    if isfield(mesh,'source') == 1
        tmp = sort(mesh.source.num);
        if length(tmp)>1
            s1 = mesh.source.coord(mesh.source.num == tmp(1),:);
            s2 = mesh.source.coord(mesh.source.num == tmp(2),:);
            plot(s1(:,1),s1(:,2),'go',s2(:,1),s2(:,2),'yo',...
                mesh.source.coord(3:end,1),...
                mesh.source.coord(3:end,2),'ro','LineWidth',2,'MarkerSize',8);
        else
            s1 = mesh.source.coord;
            plot(s1(:,1),s1(:,2),'go','LineWidth',2,'MarkerSize',8);
        end 
    end
    if isfield(mesh,'meas') == 1
        plot(mesh.meas.coord(:,1),...
            mesh.meas.coord(:,2),'bx','LineWidth',2,'MarkerSize',8);
    end
    plot(mesh.nodes(ind,1),mesh.nodes(ind,2),'c.');
    axis equal;
elseif mesh.dimension == 3
    if isfield(mesh,'source') == 1
        tmp = sort(mesh.source.num);
        if length(tmp)>1
            s1 = mesh.source.coord(mesh.source.num == tmp(1),:);
            s2 = mesh.source.coord(mesh.source.num == tmp(2),:);
            plot3(s1(:,1),s1(:,2),s1(:,3),'go',s2(:,1),s2(:,2),s2(:,3),'yo',...
                mesh.source.coord(3:end,1),...
                mesh.source.coord(3:end,2),...
                mesh.source.coord(3:end,3),'ro','LineWidth',2,'MarkerSize',8);
        else
            s1 = mesh.source.coord;
            plot3(s1(:,1),s1(:,2),s1(:,3),'go','LineWidth',2,'MarkerSize',8);
        end

    end
    if isfield(mesh,'meas') == 1
        plot3(mesh.meas.coord(:,1),...
            mesh.meas.coord(:,2),...
            mesh.meas.coord(:,3),'bx',...
            'LineWidth',2,'MarkerSize',8);
    end
    plot3(mesh.nodes(ind,1),...
        mesh.nodes(ind,2),...
        mesh.nodes(ind,3),'c.');
    axis equal;
end

if isfield(mesh,'source') && size(mesh.source.coord,1)>2
    legend('Source 1','Source 2','Sources +','Detector');
elseif isfield(mesh,'source') && size(mesh.source.coord,1)==2
    legend('Source 1','Source 2','Detector');
elseif isfield(mesh,'source') && size(mesh.source.coord,1)==1
    legend('Source 1','Detector');
else
    legend('Detector');
end

xlabel('x position (mm)');
ylabel('y position (mm)');
zlabel('z position (mm)');
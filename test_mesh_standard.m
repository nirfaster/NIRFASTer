close all
clear variables
clc

%% display properties

set(0, 'DefaultAxesFontname', 'Arial CE');
font_size = 16;
set(0, 'DefaultAxesFontsize', font_size);
set(0, 'defaulttextfontname', 'Arial CE');
set(0, 'defaulttextfontsize', font_size);
set(0, 'defaultfigurecolor', 'w');
%% BODY

%% set this script location as a current folder

where_am_i = mfilename('fullpath');
[path, ~, ~] = fileparts(where_am_i);
cd(path);

%% load meshes

mesh_2d = load_mesh('circle2000_86_stnd');
mesh_3d = load_mesh('cylinder_stnd');
disp('  *** meshes loaded')

%% save meshes

save_mesh(mesh_2d,['test' filesep 'test_2d']);
save_mesh(mesh_3d,['test' filesep 'test_3d']);
disp('  *** meshes saved')

%% re-load meshes

mesh_2d = load_mesh(['test' filesep 'test_2d']);
mesh_3d = load_mesh(['test' filesep 'test_3d']);
disp('  *** meshes re-loaded')

%% plot meshes

plotmesh(mesh_2d);
plotmesh(mesh_3d);

%% plot meshes sources and detectors

plotmesh_fiducials(mesh_2d);
plotmesh_fiducials(mesh_3d);

%% plot meshes statistics

mesh_volume_distribution(mesh_3d)

mesh_nodal_distance_distribution(mesh_2d)
mesh_nodal_distance_distribution(mesh_3d)

mesh_facets_area_distribution(mesh_2d)
mesh_facets_area_distribution(mesh_3d)

mesh_heights_distribution(mesh_2d)
mesh_heights_distribution(mesh_3d)

%% add a blob to the 2D meshes

blob.x = 0; % x position
blob.y = -10; % y position
blob.r = 10; % radius
% blob region number, max of current regions + 1
blob.region = max(unique(mesh_2d.region)) + 1;
blob.mua = 0.023; % blob absorption in mm^-1
mesh_anomaly_2d = add_blob(mesh_2d,blob);
mesh_anomaly_2d.name = [mesh_anomaly_2d.name '_anomaly'];
plotmesh(mesh_anomaly_2d)

return



mesh_3d = load_mesh('cylinder_spec');
mesh_2d = load_mesh('circle2000_86_spec');
save_mesh(mesh_2d,'test/test_2d');
save_mesh(mesh_3d,'test/test_3d');

plotmesh(mesh_2d);
plotmesh(mesh_3d);

plotmesh_fiducials(mesh_2d);
plotmesh_fiducials(mesh_3d);

mesh_volume_distribution(mesh_2d)
mesh_volume_distribution(mesh_3d)

blob.x = 0;
blob.y = -10;
blob.r = 10;
blob.region = 10;
blob.HbO = 0.023;
mesh_anom = add_blob(mesh_2d,blob);

plotmesh(mesh_anom);
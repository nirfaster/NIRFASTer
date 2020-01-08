%% Test Frequency Domain
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

%% set the 'nirfasterroot' as current folder

cd(nirfasterroot);

%% test TR on standard mesh:

%% forward data 2D
mesh = load_mesh('circle2000_86_stnd');
data_TR = femdata_TR(mesh,3e-9,25e-12);
semilogy(data_TR.time,data_TR.tpsf')
title('Log TPSF');
% time-resolved moments direct calculation
data_M = femdata_TR_moments(mesh,3);
subplot(3,1,1)
semilogy(data_M.moments(:,1));title('Log10 Intensity');
subplot(3,1,2)
plot(data_M.moments(:,2));title('mean time of flight');
subplot(3,1,3)
plot(data_M.moments(:,3));title('2nd moment (not-centralized)');

%% forward data 3D
mesh = load_mesh('cylinder_stnd');
data_TR = femdata_TR(mesh,3e-9,10e-12);
semilogy(data_TR.time,data_TR.tpsf'); title('Log TPSF');
data_M = femdata_TR_moments(mesh,3);
semilogy(data_M.moments(:,1));title('Log Intensity');
plot(data_M.moments(:,2));title('meantime');
plot(data_M.moments(:,3));title('Variance');


%% test TR on spectral mesh:

%% forward data 2D
mesh = load_mesh('circle2000_86_spec');
data_TR = femdata_TR(mesh,3e-9,25e-12);
semilogy(data_TR.time,squeeze(data_TR.tpsf(:,3,:))'); title('Log TPSF of 3rd wv');
data_M = femdata_TR_moments(mesh,3);
semilogy(squeeze(data_M.moments(:,3,1)));title('Log Intensity, 3rd wv');
plot(squeeze(data_M.moments(:,3,2)));title('meantime, 3rd wv');
plot(squeeze(data_M.moments(:,3,3)));title('Variance, 3rd wv');


%% forward data 3D
mesh = load_mesh('cylinder_spec');
data_TR = femdata_TR(mesh,3e-9,10e-12);
semilogy(data_TR.time,squeeze(data_TR.tpsf(:,3,:))'); title('Log TPSF of 3rd wv');
data_M = femdata_TR_moments(mesh,3);
semilogy(squeeze(data_M.moments(:,3,1)));title('Log Intensity, 3rd wv');
plot(squeeze(data_M.moments(:,3,2)));title('meantime, 3rd wv');
plot(squeeze(data_M.moments(:,3,3)));title('Variance, 3rd wv');

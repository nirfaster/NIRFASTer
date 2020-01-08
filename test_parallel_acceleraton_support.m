%% Test parallel acceleraton support (CUDA and OpenMP)
close all
clear variables
clc

%% BODY

%% set the 'nirfasterroot' as current folder
cd(nirfasterroot);

%% test for compatible GPUs

[is,info] = isCUDA;
% show GPUs
if is
    for ind_gpu = 1:numel(info)
        disp(['#' num2str(ind_gpu) ' GPU:'])
        disp(info(ind_gpu));
    end
else
    disp(' ')
    disp('No Nvidia GPUs supporting CUDA found. By default, code will be accelerated in parallel on CPU cores.')
    disp(' ')
end

disp(['Nvidia CUDA support: ' num2str(is)])

%% othet info

disp('Default solver:')
disp(['  ' get_solver])
disp(' ')
disp('Other supported solvers:')
disp(['  ' get_solver(solver_name_CPU)])
disp(['  ' get_solver(solver_name_backslash)])
disp(['  ' get_solver(solver_name_matlab_iterative)])
disp(' ')
disp('''help get_solver'' says:')
disp(' ')
help get_solver

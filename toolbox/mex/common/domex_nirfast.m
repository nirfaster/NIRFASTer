function domex_nirfast()
% Routine to compile the NIRFAST CORE mex files in a batch.
% Should only be run on a machine that has mex set up properly.
% It should be called whenever NIRFAST CORE mex C/C++ files are updated.
% A normal user is unlikely to run this script!
% Written by Hamid Dehghani September, 2018

mypwd = pwd;
% Find the folder where meshing mex files reside
meshing_mex = fileparts(which('nirfast'));
cd(meshing_mex);
cd('../toolbox/mex/')
mexcommand = 'mex -v ';
mexoptions = ' '; % use this for debugging

eval([mexcommand mexoptions ' IntFG.c'])
eval([mexcommand mexoptions ' IntFG_tet4.c'])
eval([mexcommand mexoptions ' IntgradFgradG.c'])
eval([mexcommand mexoptions ' IntgradFgradG_tet4.c'])
eval([mexcommand mexoptions ' create_L_fast_final.c'])
eval([mexcommand mexoptions ' distance.c'])
eval([mexcommand mexoptions ' ele_area_c.c'])
eval([mexcommand mexoptions ' gen_matrices_2d.c'])
eval([mexcommand mexoptions ' gen_matrices_2d_TR.c'])
eval([mexcommand mexoptions ' gen_matrices_3d.c'])
eval([mexcommand mexoptions ' gen_matrices_3d_TR.c'])
eval([mexcommand mexoptions ' gen_source.c'])
eval([mexcommand mexoptions ' gen_source_fl.c'])
eval([mexcommand mexoptions ' mesh_support.c'])
eval([mexcommand mexoptions ' num_int.c'])

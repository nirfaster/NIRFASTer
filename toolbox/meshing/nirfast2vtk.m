function nirfast2vtk(mesh,outfname)
% NIRFAST2VTK Saves FEM mesh in NIRFAST format into a *.vtk file.
% 
% NIRFAST2VTK(MESH,OUTFNAME) Saves the MESH in NIRFAST format into the
%   OUTFNAME.vtk file. Usefull for visualization with e.g. ParaView.
%    MESH - mesh structure in NIRFAST format or the mesh root filename
%     (with no extension). MESH can be 2D or 3D.
%    OUTFNAME - the root filename (with no extension). A path
%     (understandable by MATLAB) or filename if saving in the curent
%     directory. A OUTFNAME.vtk file will be created.
% 
%   Following data fields of the MESH structure will be saved if present
%   (all are optional):
%  FOR ANY USER DATA:
%    'mesh.phi'      as 'phi_1', 'phi_2' ... 'phi_n' for n columns of
%                       'mesh.phi'. Number of rows of the 'mesh.phi' should
%                       be equal number of MESH nodes. E.g. use it to save
%                       photon fluence rate at n sources.
%    'mesh.jacobian' as 'Jacobian_1' ... 'Jacobian_n' for n columns of
%                       'mesh.jacobian'. Number of rows of the
%                       'mesh.jacobian' should be equal number of MESH
%                       nodes. E.g. use it to save sensitivity matrices for
%                       source-detector pairs. 'mesh.jacobian' works
%                       exactly the same as 'mesh.phi' but has a different
%                       name giving you more flexibility.
%  MESH SPECIFIC DATA:
%    'mesh.region'   as 'region/-'
%    'mesh.mua       as 'mua/mm^-1';
%    'mesh.mus'      as 'mus/mm^-1';
%    'mesh.kappa'    as 'kappa/mm';
%    'mesh.muax'     as 'muax';
%    'mesh.musx'     as 'musx';
%    'mesh.muam'     as 'muam';
%    'mesh.musm'     as 'musm';
%    'mesh.muaf'     as 'muaf';
%    'mesh.eta'      as 'eta';
%    'mesh.etamuaf'  as 'etamuaf';
%    'mesh.tau'      as 'tau';
%    'mesh.sa'       as 'sa/-';
%    'mesh.sp'       as 'sp/-';
%    'mesh.support'  as 'support/mm^3';
%    'mesh.aDb'      as 'aDb/mm^2s^-1';
%    'mesh.conc'     as 'Hbo/mM', 'deoxyHb/mM', 'Water/percent' and other chromphores if specified in 'mesh.chromscattlist'
%    'mesh.conc'     as 'HbT/mM', ('Hbo/mM' + 'deoxyHb/mM') if both presend
%    'mesh.conc'     as 'StO2/percent', (100*'Hbo/mM'./'HbT/mM') if 'Hbo/mM' and 'deoxyHb/mM' present
%    'mesh.target'   as 'target'
%    'mesh.recon'    as 'recon'

% 
% NIRFAST2VTK(OUTFNAME,MESH) Any order of arguments allowed.
% 
% EXAMPLE
%   This code will save amplitude and phase data for visualization: 
%     mesh = load_mesh('cylinder_stnd');
%     data = femdata_stnd_FD(mesh,100e6);
%     mesh.phi = abs(data.phi);
%     mesh.jacobian = angle(data.phi);
%     nirfast2vtk(mesh,'test')
%   Now load the 'test.vtk' file with e.g. the 'ParaView' software.
% 
% See also SOURCES2VTK, DETECTORS2VTK, LOAD_MESH.
% 
%   Part of NIRFAST package.

%% check in/out

narginchk(2,2);
nargoutchk(0,0);

% check if the user reversed the inputs
if isstruct(outfname) && ischar(mesh)
    fn_copy = mesh;
    mesh = outfname;
    outfname = fn_copy;
elseif ~ischar(outfname) || ~isstruct(mesh)
    error('Bad arguments. Please see help on how to use this function.')
end


%% load mesh if needed
if ischar(mesh) || isstring(mesh)
  mesh = load_mesh(mesh);
end

%% get mesh properties as specified per node

% list of parameters names
listsolfnames = {};
% matrix of the values
soldata = [];

if isfield(mesh,'region')
    listsolfnames{end+1} = 'region/-';
    soldata(:,end+1) = mesh.region(:,1);
end
if isfield(mesh,'mua')
    if size(mesh.mua,2) > 1
        for ind_experiment = 1:size(mesh.mua,2)
            listsolfnames{end+1} = ['mua_' num2str(ind_experiment)]; %#ok<*AGROW>
            soldata(:,end+1) = double(mesh.mua(:,ind_experiment));
        end
    else
        listsolfnames{end+1} = 'mua/mm^-1';
        soldata(:,end+1) = mesh.mua(:,1);
    end
end
if isfield(mesh,'mus')
    if size(mesh.mus,2) > 1
        for ind_experiment = 1:size(mesh.mus,2)
            listsolfnames{end+1} = ['mus_' num2str(ind_experiment)];
            soldata(:,end+1) = double(mesh.mus(:,ind_experiment));
        end
    else
        listsolfnames{end+1} = 'mus/mm^-1';
        soldata(:,end+1) = mesh.mus(:,1);
    end
end
if isfield(mesh,'kappa')
    if size(mesh.kappa,2) > 1
        for ind_experiment = 1:size(mesh.kappa,2)
            listsolfnames{end+1} = ['kappa_' num2str(ind_experiment)];
            soldata(:,end+1) = double(mesh.kappa(:,ind_experiment));
        end
    else
        listsolfnames{end+1} = 'kappa/mm';
        soldata(:,end+1) = mesh.kappa(:,1);
    end
end
if isfield(mesh,'muax')
    listsolfnames{end+1} = 'muax';
    soldata(:,end+1) = mesh.muax;
end
if isfield(mesh,'musx')
    listsolfnames{end+1} = 'musx';
    soldata(:,end+1) = mesh.musx;
end
if isfield(mesh,'muam')
    listsolfnames{end+1} = 'muam';
    soldata(:,end+1) = mesh.muam;
end
if isfield(mesh,'musm')
    listsolfnames{end+1} = 'musm';
    soldata(:,end+1) = mesh.musm;
end
if isfield(mesh,'muaf')
    listsolfnames{end+1} = 'muaf';
    soldata(:,end+1) = mesh.muaf;
    if isfield(mesh,'eta') && ~isfield(mesh,'etamuaf')
        mesh.etamuaf = mesh.eta.*mesh.muaf;
    end
end
if isfield(mesh,'eta')
    listsolfnames{end+1} = 'eta';
    soldata(:,end+1) = mesh.eta;
end
if isfield(mesh,'etamuaf')
    listsolfnames{end+1} = 'etamuaf';
    soldata(:,end+1) = mesh.etamuaf;
end
if isfield(mesh,'tau')
    listsolfnames{end+1} = 'tau';
    soldata(:,end+1) = mesh.tau;
end
if isfield(mesh,'sa')
    listsolfnames{end+1} = 'sa/-';
    soldata(:,end+1) = mesh.sa;
end
if isfield(mesh,'sp')
    listsolfnames{end+1} = 'sp/-';
    soldata(:,end+1) = mesh.sp;
end
if isfield(mesh,'support')
    % volume of elements surrounding a node, sort of mesh resolution spatial distribution
    listsolfnames{end+1} = 'support/mm^3';
    soldata(:,end+1) = mesh.support;
end
if isfield(mesh,'aDb')
    % DCS flow, brownian motion coefficient (mean particles displacement)
    listsolfnames{end+1} = 'aDb/mm^2s^-1';
    soldata(:,end+1) = mesh.aDb;
end
if isfield(mesh,'target')
    % target
    listsolfnames{end+1} = 'target/AU';
    soldata(:,end+1) = mesh.target;
end
if isfield(mesh,'recon')
    % target
    listsolfnames{end+1} = 'recon/AU';
    soldata(:,end+1) = mesh.recon;
end
if isfield(mesh,'phi')
    % any field value, usfull to image photon fluence rate
    % check to be sure if all will end up nice
    if size(mesh.phi,1) == size(mesh.nodes,1)
        for ind_source = 1:size(mesh.phi,2)
            listsolfnames{end+1} = ['phi_' num2str(ind_source)];
            soldata(:,end+1) = double(mesh.phi(:,ind_source));
        end
    else
        warning('''mesh.phi'' size does not match number of mesh nodes. ''mesh.phi'' not saved.');
    end
end
if isfield(mesh,'jacobian')
    % any field value, usfull to image sensitivity matrices
    if size(mesh.jacobian,1) == size(mesh.nodes,1)
        for ind_pair = 1:size(mesh.jacobian,2)
            listsolfnames{end+1} = ['Jacobian_' num2str(ind_pair)];
            soldata(:,end+1) = double(mesh.jacobian(:,ind_pair));
        end
    else
        warning('''mesh.jacobian'' size does not match number of mesh nodes. ''mesh.jacobian'' not saved.');
    end
end

% index of oxy and deohy haemoglobins
hbo_loc = -1;
deoxyhb_loc = -1;

if isfield(mesh,'conc') && isfield (mesh,'chromscattlist')
    
    mask_chromophores = ~strcmp(mesh.chromscattlist,'S-Amplitude') & ...
                        ~strcmp(mesh.chromscattlist,'S-Power');
    
    for i = 1:length(mesh.chromscattlist)
        if mask_chromophores(i) && strcmp(mesh.chromscattlist{i},'HbO')
            mesh.conc(:,i) = 1000*mesh.conc(:,i);
            hbo_loc = i;
%             disp('HbO converted to millimolar');
            listsolfnames{end+1} = [mesh.chromscattlist{i} '/mM'];
            soldata(:,end+1) = mesh.conc(:,i);
        elseif mask_chromophores(i) && strcmp(mesh.chromscattlist{i},'deoxyHb')
            mesh.conc(:,i) = 1000*mesh.conc(:,i);
            deoxyhb_loc = i;
%             disp('deoxyHb converted to millimolar');
            listsolfnames{end+1} = [mesh.chromscattlist{i} '/mM'];
            soldata(:,end+1) = mesh.conc(:,i);
        elseif mask_chromophores(i) && strcmp(mesh.chromscattlist{i},'Water')
            mesh.conc(:,i) = 100*mesh.conc(:,i);
%             disp('Water converted to %');
            listsolfnames{end+1} = [mesh.chromscattlist{i} '/percent'];
            soldata(:,end+1) = mesh.conc(:,i);
        elseif mask_chromophores(i)
            listsolfnames{end+1} = mesh.chromscattlist{i};
            soldata(:,end+1) = mesh.conc(:,i);
        end
    end
end

% if oxy and deoxy haemoglobins are present
if (hbo_loc > 0) && (deoxyhb_loc > 0)
    listsolfnames{end+1} = 'HbT/mM';
    soldata(:,end+1) = mesh.conc(:,hbo_loc) + mesh.conc(:,deoxyhb_loc);
    listsolfnames{end+1} = 'StO2/percent';
    soldata(:,end+1) = 100*mesh.conc(:,hbo_loc)./soldata(:,end);
%     disp('HbT is in millimolar');
%     disp('StO2 is in %');
end

%% write to vtk

numnodes = size(mesh.nodes,1);
numelems = size(mesh.elements,1);

outfname = add_extension(outfname,'.vtk');

% write to system temporary folder if we do not have write acces to the requested location
if ~canwrite(outfname)
    [path, fn, ext] = fileparts(outfname);
    outfname = fullfile(tempdir, strcat(fn,ext));
    warning(['No write access to ''' path ''', writing here instead: ''' outfname '''']);
end

fid = fopen(outfname,'w');

%define an VTK header for FEM mesh representation
line0 = '# vtk DataFile Version 2.0';
line1 = 'NIRFAST mesh with solutions';
line2 = 'ASCII';
line3 = 'DATASET UNSTRUCTURED_GRID';
fprintf(fid,'%s\n%s\n%s\n%s\n',line0,line1,line2,line3);

% nodes
line4 = ['POINTS ', num2str(numnodes), ' double'];
fprintf(fid,'%s\n',line4);

%elements
if mesh.dimension == 2
    fprintf(fid, '%f %f %f\n', double([mesh.nodes(:,1:mesh.dimension) zeros(numnodes,1)]'));
    line5 = ['CELLS ',num2str(numelems),' ',num2str(numelems*3+numelems)]; %connectivity maps
    fprintf(fid,'%s\n',line5);
    fprintf(fid,'%d %d %d %d\n',[3*ones(numelems,1) mesh.elements-1]');
    line6 = ['CELL_TYPES ', num2str(numelems)];
    fprintf(fid,'%s\n',line6);
    fprintf(fid,'%d\n', ones(numelems,1)*5); %specify the mesh basis 5->triangle for all connectivity maps 
else
    fprintf(fid, '%f %f %f\n', mesh.nodes');
    line5 = ['CELLS ',num2str(numelems),' ',num2str(numelems*4+numelems)]; %connectivity maps
    fprintf(fid,'%s\n',line5);
    fprintf(fid,'%d %d %d %d %d\n',[4*ones(numelems,1) mesh.elements-1]');
    line6 = ['CELL_TYPES ', num2str(numelems)];
    fprintf(fid,'%s\n',line6);
    fprintf(fid,'%d\n', ones(numelems,1)*10); %specify the mesh basis 10->tetrahedral for all connectivity maps 
end

% the data
fprintf(fid,'%s\n',['POINT_DATA ',num2str(numnodes)]); %specify the data that follows is defined on the nodes
for i = 1:size(soldata,2)
%     fprintf(fid,'%s\n',['SCALARS ', listsolfnames{i}, ' double 1']);
    fprintf(fid,'%s\n',['SCALARS ', listsolfnames{i}, ' double']);
    fprintf(fid,'%s\n','LOOKUP_TABLE default');
    fprintf(fid,'%e\n', double(soldata(:,i)));
end

% close the file
fclose(fid);

end

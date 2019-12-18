function sources2vtk(mesh,output_fn,varargin)
% SOURCES2VTK Saves sources of a FEM mesh in NIRFAST format into a *.vtk file.
% 
% SOURCES2VTK(MESH,OUTPUT_FN) Saves sources (MESH.source.coord) of a FEM
%   MESH in NIRFAST format into a *.vtk file. Usefull for visualization
%   with e.g. ParaView.
%    MESH - mesh structure in NIRFAST format or the mesh root filename
%     (with no extension). MESH can be 2D or 3D.
%    OUTPUT_FN - the root filename (with no extension). A path
%     (understandable by MATLAB) or filename if saving in the curent
%     directory. A OUTPUT_FN.vtk file will be created.
%   By default the sources coordinates are snapped to the MESH surface.
% 
% SOURCES2VTK(MESH,OUTPUT_FN,SNAP) The SNAP is an optional parameter
%   allowing to control if the sources coordinates should be snapped to the
%   MESH surface. SNAP can have following values:
%    'snap' - (default) project sources onto the MESH surface
%    'nosnap' - use unchanged sources coordinates
%   All other values of SNAP are ignored and the default 'snap' is used.
% 
% SOURCES2VTK(OUTPUT_FN,MESH) Any order of arguments allowed.
% SOURCES2VTK(OUTPUT_FN,MESH,SNAP) Any order of required arguments allowed.
% 
% See also DETECTORS2VTK, LOAD_MESH.
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

%% check in/out

narginchk(2,3);
nargoutchk(0,0);

% check if the user reversed the inputs
if isstruct(output_fn) && ischar(mesh)
    fn_copy = mesh;
    mesh = output_fn;
    output_fn = fn_copy;
elseif ~ischar(output_fn) || ~isstruct(mesh)
    error('Bad arguments. Please see help on how to use this function.')
end

% if we want to snap point to the mesh surface, the default value
snap2surface = true;

if ~isempty(varargin)
    if strcmp(varargin{1}, 'snap')
        snap2surface = true;
    elseif strcmp(varargin{1}, 'nosnap')
        snap2surface = false;
    else
        warning(['Unrecognized function parameter: ''' num2str(varargin{1}) '''. Default ''snap'' is used.'])
    end
end


%% BODY

%% open file

outfname = add_extension(output_fn,'.vtk');

% write to system temporary folder if we do not have write acces to the requested location
if ~canwrite(outfname)
    [path, fn, ext] = fileparts(outfname);
    outfname = fullfile(tempdir, strcat(fn,ext));
    warning(['No write access to ''' path ''', writing here instead: ''' outfname '''']);
end

fid = fopen(outfname,'w');

%% write header

%define an VTK header for FEM mesh representation
line0 = '# vtk DataFile Version 2.0';
line1 = '# NIRFAST mesh sources';
line2 = 'ASCII';
line3 = 'DATASET POLYDATA';
fprintf(fid,'%s\n%s\n%s\n%s\n',line0,line1,line2,line3);

%% snap sources to the surface

if snap2surface
    disp('Moving sources');
    % copy coordinates of sources to detectors
    mesh.meas.coord = mesh.source.coord;
    % this will move sources onto surface (source points are detector points now)
    mesh = move_detector(mesh);
    % sources! coordinates to save
    coords = mesh.meas.coord;
else
    coords = mesh.source.coord;
end

%% set on XY plane for 2D meshes
if mesh.dimension == 2
    coords = [coords(:,1:mesh.dimension) zeros(size(coords,1),1)];
end

%% write points data in the VTK format

line8 = ['POINTS ' num2str(size(coords,1)) ' float'];
fprintf(fid,'%s\n',line8);
fprintf(fid, '%f %f %f\n', coords');
line7 = ['VERTICES 1 ' num2str(size(coords,1) + 1)];
fprintf(fid,'%s\n',line7);
fprintf(fid, '%d ', [size(coords,1) (1:size(coords,1))-1]);
fprintf(fid,'%s\n','');

% close the file
fclose(fid);

end


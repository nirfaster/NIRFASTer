function save_mesh(mesh,fn)
% SAVE_MESH Saves a FEM mesh in NIRFAST format into ASCII files.
% 
% SYNTAX:
%  SAVE_MESH(MESH,FN)
%  SAVE_MESH(FN,MESH)
% 
% SAVE_MESH(MESH,FN) Saves the MESH structure of FEM mesh in NIRFAST
%   format into MATLAB *.mat file or ASCII files.
%    FN - the root filename of the mesh (with no extension). A path
%     (understandable by MATLAB) or filename if saveing in the curent
%     directory. FN will contain following ASCII files: required *.nodes,
%     *. elem, *.param and optional *.source, *.meas, *.link, *.region,
%     *.excoeff. 
%     E.g. FN = '/dir1/dir2/my_mesh' or FN = 'my_mesh'.
%   FN - the root file name of the mesh. A path (understandable by MATLAB)
%     or filename if saveing in the curent directory. MESH will be
%     saved as follows: 
%        - As a *.mat file - if FN is given with the '.mat' extension.
%        - As list of ASCII filles if FN.mat not present. FN must contain
%          following ASCII files: FN.node, FN. elem, FN.param. Optional
%          ASCII files: FN.source, FN.meas, FN.link, FN.region, FN.excoeff.
%          E.g. FN='/dir1/dir2/my_mesh', FN='cylinder_stnd',
%          FN='cylinder_stnd.mat'.
%   See help for 'load_mesh' for details on MESH fields and format of the
%   ASCII files.
% 
% SAVE_MESH(FN,MESH) Any order of arguments allowed.
% 
%   TO_DO: 'blt'   - bioluminescence (TO_DO: not implemented yet)
% 
% WARNING:
%   As of this NIRFAST verison support for SPN and BEM meshes has been
%   removed. The last NIRFAST version that supports SPN and BEM is NIRFAST
%   v9.1.
% 
% See also LOAD_MESH. 
% 
%   Part of NIRFAST package.

%% check in/out

narginchk(2,2);
nargoutchk(0,0);

% check if the user reversed the inputs
if isstruct(fn) && ischar(mesh)
    fn_copy = mesh;
    mesh = fn;
    fn = fn_copy;
elseif ~ischar(fn) || ~isstruct(mesh)
    error('Bad arguments. Please see help on how to use this function.')
end

%% save as MATLAB fn.mat if requested

[~,name,ext] = fileparts(fn);

% if MATLAB *.mat file in use, just save and return
if strcmp(ext,'.mat')
    save([name '.mat'], 'mesh');
    return;
end

%% saving fn.node file

% mysave([fn '.node'],[mesh.bndvtx mesh.nodes]);
dlmwrite([fn '.node'],[mesh.bndvtx mesh.nodes],'delimiter','\t','precision','%g');

%% saving fn.elem file

% mysave([fn '.elem'],mesh.elements);
dlmwrite([fn '.elem'],mesh.elements,'delimiter','\t','precision','%g');

%% saving fn.param file

if strcmp(mesh.type,'stnd')
    mesh.kappa = 1./(3.*(mesh.mua+mesh.mus));
    data = [mesh.mua mesh.kappa mesh.ri];
elseif strcmp(mesh.type,'fluor')
    mesh.kappax = 1./(3.*(mesh.muax+mesh.musx));
    mesh.kappam = 1./(3.*(mesh.muam+mesh.musm));
    data = [mesh.muax mesh.kappax mesh.ri mesh.muam ...
        mesh.kappam mesh.muaf mesh.eta mesh.tau];
elseif strcmp(mesh.type,'spec')
    % get chromophores indexes
    mask_chromophores = ~strcmpi(mesh.chromscattlist,'S-Amplitude') & ...
                        ~strcmpi(mesh.chromscattlist,'S-Power');
    % get chromophores concentrations
    data = mesh.conc(:,mask_chromophores);
    % append scattering amplituse and power
    data = [data mesh.sa mesh.sp];
elseif strcmp(mesh.type,'dcs')
    mesh.kappa = 1./(3.*(mesh.mua+mesh.mus));
    if ~isfield(mesh,'Db') || ~isfield(mesh,'alpha') || ~isfield(mesh,'wv_DCS')
        % check if we have the DCS fields, issue warning and save as for 'stnd' if not
        warning(['The mesh is of ''dcs'' type. However, ''Db'' (brownian motion displacement) and/or ''alpha'' (concentration) '...
            'and/or ''wv_DCS'' (DCS wavelength) fields are missing. No flow parameters will be saved'])
        data = [mesh.mua mesh.kappa mesh.ri];
    else
        data = [mesh.mua mesh.kappa mesh.ri mesh.alpha mesh.Db];
    end
elseif strcmp(mesh.type,'blt')
    error('''blt'' mesh type not implemented yet');
end

% open .param file
fid = fopen([fn '.param'],'w');
% write mesh type
fprintf(fid,'%s\n',mesh.type);
if strcmp(mesh.type,'spec')
    % write chromophores and scattering parameters names if needed
    for i = 1:numel(mesh.chromscattlist)
        fprintf(fid,'%s\n',mesh.chromscattlist{i});
    end
end
if strcmp(mesh.type,'dcs')
    if ~isfield(mesh,'wv_DCS')
        % check if we have the DCS wavelength field
        warning('The mesh is of ''dcs'' type. However, ''wv_DCS'' (DCS wavelength) field is missing. Default 785nm will be saved.')
        % write DCS wavelength
        fprintf(fid,'%s\n',['wv_DCS:' num2str(785)]);
    else
        % write DCS wavelength
        fprintf(fid,'%s\n',['wv_DCS:' num2str(mesh.wv_DCS)]);
    end
end

% close the file
fclose(fid);
% now append the data matrix
dlmwrite([fn '.param'],data,'delimiter',' ','precision','%g','-append');

%% save extinction file for spec mesh type

if strcmp(mesh.type,'spec')
    
    % check if mesh.wv is the right format
    if size(mesh.wv,1) < size(mesh.wv,2)
        mesh.wv = mesh.wv';
    end
    
    data = [mesh.wv mesh.excoef];
    
    % open .param file
    fid = fopen([fn '.excoef'],'w');
    % get chromophores indexes
    mask_chromophores = ~strcmpi(mesh.chromscattlist,'S-Amplitude') & ...
                        ~strcmpi(mesh.chromscattlist,'S-Power');
    % write chromophores and scattering parameters names if needed
    for i = 1:numel(mesh.chromscattlist)
        if mask_chromophores(i)
            fprintf(fid,'%s\n',mesh.chromscattlist{i});
        end
    end
    % close the file
    fclose(fid);
    % now append the data matrix
    dlmwrite([fn '.excoef'],data,'delimiter',' ','precision','%g','-append');
else
    % delete if file exists and no mesh data to save
    if exist([fn '.excoef'],'file') == 2
        delete([fn '.excoef']);
    end    
end

%% saving fn.region file

if isfield(mesh,'region')
    %mysave([fn '.region'],mesh.region);
    dlmwrite([fn '.region'],mesh.region,'delimiter','\t','precision','%g');
else
    % delete if file exists and no mesh data to save
    if exist([fn '.region'],'file') == 2
        delete([fn '.region']);
    end
end

%% saving fn.source file

if isfield(mesh,'source')
    
    %check formatting
    if ~isfield(mesh.source,'fixed') || ~isfield(mesh.source,'num') || ...
            ~isfield(mesh.source,'coord') || ~isfield(mesh.source,'fwhm')
        error('Bad ''mesh.source'' format. Requirds fields not present.')
    end
    
    % open .source file
    fid = fopen([fn '.source'],'w');
    % write if fixed
    if logical(mesh.source.fixed)
        fprintf(fid,'%s\n','fixed');
    end
    % write if source integration functions present
    if isfield(mesh.source,'int_func') && (size(mesh.source.int_func,1) == size(mesh.source.coord,1))
        % prepare column names
        if mesh.dimension == 2
            fprintf(fid, '%s %s %s %s %s %s %s %s\n',...
                         'num','x','y','fwhm','ele','ip1','ip2','ip3');
        elseif mesh.dimension == 3
            fprintf(fid, '%s %s %s %s %s %s %s %s %s %s\n',...
                         'num','x','y','z','fwhm','ele','ip1','ip2','ip3','ip4');
        end
        % prepare data
        data = [mesh.source.num, mesh.source.coord(:,1:mesh.dimension), mesh.source.fwhm mesh.source.int_func];
    else
        % if no integration functions
        % prepare column names
        if mesh.dimension == 2
            fprintf(fid, '%s %s %s %s\n',...
                         'num','x','y','fwhm');
        elseif mesh.dimension == 3
            fprintf(fid, '%s %s %s %s %s\n',...
                         'num','x','y','z','fwhm');
        end
        % prepare data
        data = [mesh.source.num, mesh.source.coord(:,1:mesh.dimension), mesh.source.fwhm];
    end
    % close the file
    fclose(fid);
    % now append the data matrix
    dlmwrite([fn '.source'],data,'delimiter',' ','precision','%.16g','-append');
else
    % delete if file exists and no mesh data to save
    if exist([fn '.source'],'file') == 2
        delete([fn '.source']);
    end
end

%% saving fn.meas file

if isfield(mesh,'meas')

    %check formatting
    if ~isfield(mesh.meas,'fixed') || ~isfield(mesh.meas,'num') || ~isfield(mesh.meas,'coord')
        error('Bad ''mesh.meas'' format. Requirds fields not present.')
    end
       
    % open .meas file
    fid = fopen([fn '.meas'],'w');
    % write if fixed
    if logical(mesh.meas.fixed)
        fprintf(fid,'%s\n','fixed');
    end
    % write if detector integration functions present
    if isfield(mesh.meas,'int_func') && (size(mesh.meas.int_func,1) == size(mesh.meas.coord,1))
        % prepare column names
        if mesh.dimension == 2
            fprintf(fid, '%s %s %s %s %s %s %s\n',...
                         'num','x','y','ele','ip1','ip2','ip3');
        elseif mesh.dimension == 3
            fprintf(fid, '%s %s %s %s %s %s %s %s %s\n',...
                         'num','x','y','z','ele','ip1','ip2','ip3','ip4');
        end
        % prepare data
        data = [mesh.meas.num, mesh.meas.coord(:,1:mesh.dimension), mesh.meas.int_func];
    else
        % if no integration functions
        % prepare column names
        if mesh.dimension == 2
            fprintf(fid, '%s %s %s\n',...
                         'num','x','y');
        elseif mesh.dimension == 3
            fprintf(fid, '%s %s %s %s\n',...
                         'num','x','y','z');
        end
        % prepare data
        data = [mesh.meas.num, mesh.meas.coord(:,1:mesh.dimension)];
    end
    % close the file
    fclose(fid);
    % now append the data matrix
    dlmwrite([fn '.meas'],data,'delimiter',' ','precision','%.16g','-append');
else
    % delete if file exists and no mesh data to save
    if exist([fn '.meas'],'file') == 2
        delete([fn '.meas']);
    end
end

%% saving fn.link file

if isfield(mesh,'link')
    % open .link file
    fid = fopen([fn '.link'],'w');
    % write mesh type
    fprintf(fid,'%s\n','source detector active');
    % close the file
    fclose(fid);
    % now append the data matrix
    dlmwrite([fn '.link'],[mesh.link],'delimiter',' ','precision','%g','-append');
else
    % delete if file exists and no mesh data to save
    if exist([fn '.link'],'file') == 2
        delete([fn '.link']);
    end
end

%% write *.vtk files for visualization

% nirfast2vtk(mesh,[fn,'.vtk']);
% sources2vtk(mesh,[fn,'.source.vtk']);
% detectors2vtk(mesh,[fn,'.meas.vtk']);

end

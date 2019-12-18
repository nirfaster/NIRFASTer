function mesh = load_mesh(fn)
% LOAD_MESH Loads meshes from specified root file name.
%
% SYNTAX:
%  MESH = LOAD_MESH(FN)
%  LOAD_MESH(FN)
%
% [MESH] = LOAD_MESH(FN) Creates the MESH structure of FEM mesh in NIRFAST
%   format. Data is loaded from the root file name FN.
%   FN - the root file name of the mesh. The filename should be available
%        in the MATLAB path. Otherwise, specify the full path. FN will be
%        interpreted as follows:
%        - As *.mat file - if FN.mat exists or FN given with *.mat
%          extension. The 'FN.mat' file is required to have the 'mesh'
%          structure stored inside.
%        - As list of ASCII filles if FN.mat not present. FN must contain
%          following ASCII files: *.node, *. elem, *.param. Optional ASCII
%          files: *.source, *.meas, *.link, *.region, *.excoeff.
%
%          FN='/dir1/dir2/my_mesh', FN='cylinder_stnd' or
%          FN='cylinder_stnd.mat'.
%
% ##### MESH structure #####
%
% ----common for all MESH types----
%
%  MESH.type - mesh type as specified in the *.param file can be either:
% 	            'stnd'  - standard (single wavelength light sources)
%               'spec'  - spectral (multi wavelength light sources)
%               'fluor' - fluorescence
%               'dcs'   - diffuse correlation spectroscopy (single
%                         wavelength, long-coherence length light sources)
%               'blt'   - bioluminescence
%  MESH.name - the mesh file name
%  MESH.dimension - a scalar equals 2 or 3, Euclidean space dimension
%  MESH.nodes - in millimeters, n by d matrix of Cartesian coordinates of n
%               nodes in 2D (x,y) or 3D (x,y,z)
%  MESH.elements - m by d+1 matrix of connectivity list of MESH elements.
%                  Where m is the number of elements and d is the
%                  'MESH.dimension' (d=2 triangles, d=3 tetrahedrons). A row
%                  entry of 'MESH.elements' defines an element as d+1 long
%                  list of indexes of nodes as stored in the 'MESH.nodes'.
%                  e.g. entry [7 13 3 19] defines a linear tetrahedron where
%                  (x,y,z) coordinates of vertices are stored at
%                  MESH.nodes([7 13 3 19],:).
%  MESH.bndvtx - n by 1 vector of flags (0 or 1) indicating if a node
%                belongs to the external boundary (1), n is number of nodes.
%  MESH.region - n by 1 vector of nodes labels showing where nodes belong,
%                e.g. 1 is scalp, 2 is skull, etc., n is number of nodes
%  MESH.ri - unitless, n by 1 vector of refractive index values, n is
%            number of nodes
%  MESH.c - in mm/s, n by 1 vector of speed of light in medium (speed of
%           light in vaccum normalized by MESH.ri, n is number of nodes
%  MESH.ksi - unitless, n by 1 vector of the photon fluence rate scale
%             factor on the mesh-outside_mesh boundary as derived from
%             Fresenel's law, n is number of nodes
%  MESH.element_area - in mm^3 (3D) or mm^2 (2D), m by 1 vector of
%                      'MESH.elements' volumes or area if in 2D, m is
%                      number of elements
%  MESH.support - in mm^3 (3D) or mm^2 (2D), n by 1 vector of volume or
%                 area if in 2D of all MESH elements connected ta a given
%                 node, n is number of nodes
%  MESH.source - structure of light sources
%     MESH.source.fixed - a scalar, 0 or 1 flag. Sources are fixed in space
%                         (1) or allowed to be checked and moved (0) onto
%                         the MESH surface and then one scattering distance
%                         inside the MESH.
%     MESH.source.num - s by 1 vector of sources identification numbers,
%                       where s is the number of sources. Typicaly,
%                       'MESH.source.num' equals 1:s.
%     MESH.source.coord - in millimeters, s by d matrix of Cartesian
%                         coordinates of s sources in 2D (d=2: x,y) or 3D
%                         (d=3: x,y,z), d equals 'MESH.dimension'
%     MESH.source.fwhm - in millimeters, s by 1 vector of s sources spatial
%                        profile FWHM. Controls spatial size of sources. If
%                        0, source is a point (pencl beam). If >0, source
%                        intensity has a spatial Gaussian distribution with
%                        the FWHM as specified.
%     MESH.source.int_func - s by d+2 matrix of barycentric coordinates of
%                            s sources, where d equals 'MESH.dimension'. A
%                            row entry will have a following structure:
%                            [elem_index, b_coordiantes], where elem_index
%                            is index of the MESH element where the source
%                            belongs (row index of the 'MESH.elements') and
%                            b_coordiantes holds barycentric coordinates of
%                            the source within the element. Number of the
%                            barycentric coordinates id d+1. Sum of the
%                            coordinates is 1.
%  MESH.meas - structure of light detectors
%     MESH.meas.fixed - see 'MESH.source.fixed'
%     MESH.meas.num - k by 1 vector of detectors identification numbers,
%                     where k is the number of detectors. Typicaly,
%                     'MESH.meas.num' equals 1:k.
%     MESH.meas.coord - in millimeters, k by d matrix of Cartesian
%                       coordinates of k detectors in 2D (d=2: x,y) or 3D
%                       (d=3: x,y,z), d equals 'MESH.dimension'
%     MESH.meas.int_func - k by d+2 matrix of barycentric coordinates of k
%                          sources, where d equals 'MESH.dimension'. See
%                          'MESH.source.int_func' for details on structure.
%
% ----'stnd' MESH types specific----
%
%  MESH.mua - in mm^-1, n by 1 vector of the absorption coefficient at n
%             MESH nodes
%  MESH.mus - in mm^-1, n by 1 vector of the reduced scattering coefficient
%             at n MESH nodes
%  MESH.kappa - in mm, n by 1 vector of the diffusion coefficient at n MESH
%               nodes, defined as (3*(MESH.mua + MESH.mus))^-1
%  MESH.link - p by 3 matrix showing how to link sources and detectors as p
%              pairs. A row entry of the MESH.link will look as follows:
%              [source_number detector_number is_pair], where source_number
%              is the identification number as in 'MESH.source.num',
%              detector_number is the identification number as in
%              'MESH.meas.num' and is_pair is 0 or 1 showing if a given
%              pair is used (1) or disabled (0).
%
% ----'spec' MESH types specific----
%
%  MESH.chromscattlist - cell array of names of wavelength dependent MESH
%                        parameters as saved in the '*.param' file.
%  MESH.conc - n by c matrix of cosituents concentrations for at n MESH
%              nodes for c constituents. If not user-redefined, colums
%              meaning is as follows: 1st - oxygenated haemoglobin
%              concentration in mM (milli molar), 2nd - deoxygenated
%              haemoglobin concentration in mM (milli molar) and 3rd -
%              water volume fraction from 0 (no water) to 1 (100% water).
%  MESH.sa - n by 1 vector of scattering amplitude at n nodes, unitless, as
%            defined by Mie theory, reference wavelength 1000nm
%  MESH.sp - n by 1 vector of scattering power at n nodes, unitless, as
%            defined by Mie theory, reference wavelength 1000nm
%  MESH.wv - in nm (nano meters), w elements vector of light sources
%            wavelengts
%  MESH.excoef - w by c matrix of extinction/attenuation coefficients at w
%                wavelengths as defined in 'MESH.wv' for c constituents at
%                concetrations as specified in MESH.conc. If not
%                user-redefined, colums meaning is as follows: 1st -
%                oxygenated haemoglobin extinction coefficient expressed in
%                ln(10)*mm^-1mM^-1 (mM - milli molar (molar -
%                moles/liter)), 2nd - deoxygenated haemoglobin extinction
%                coefficient exoressed in the same units as the oxygenated
%                haemoglobin and 3rd - water absorption coefficient
%                expressed in mm^-1.
%  MESH.link - p by 2+w matrix showing how to link sources and detectors as
%              p pairs. A row entry of the MESH.link will look as follows:
%              [source_number detector_number is_pair@lambda_1 is_pair@lambda_2 ... is_pair@lambda_w],
%              where source_number is the identification number as in
%              'MESH.source.num', detector_number is the identification
%              number as in 'MESH.meas.num' and is_pair@lambda_1 to
%              is_pair@lambda_n is 0 or 1 showing if a given pair is used
%              (1) or disabled (0) at a given wavelength. Wavelengt as
%              specified in 'MESH.wv'.
%
% ----'dcs' MESH types specific----
%
%  MESH.Db - n by k matrix of nodal values of k brownian motion
%            displacement components expressed in mm^2*s^-1.
%  MESH.alpha - n by k matrix of nodal values of concentration of the
%               MESH.Db components expressed as fraction between 0 and 1.
%               Sum of concentrations along the components (k) dimmension
%               should be 1.
%  MESH.aDb - n by 1 matrix of nodal values of the 'flow index' expressed
%             in mm^2*s^-1. It representing the brownian motion
%             displacement (Db) weighted by the concentration (alpha). It
%             equals linear combination of the 'Db' and 'alpha': MESH.aDb =
%             sum(MESH.Db.*MESH.alpha,2).
%  MESH.wv_DCS - a scalar, the DCS wavelength in nm, e.g. 785
%  MESH.link - p by 3 matrix showing how to link sources and detectors as p
%              pairs. Structure as for the 'stnd' MESH type.
%
% ----'fluo' MESH types specific----
%  MESH.muax - n by 1 matrix of nodal values of the absorption at excitation
%             wavelength
%  MESH.musx - n by 1 matrix of nodal values of the reduced scatter at excitation
%             wavelength
%  MESH.kappax - n by 1 matrix of nodal values of the diffusion coefficient at excitation
%             wavelength
%  MESH.muam - n by 1 matrix of nodal values of the absorption at emission
%             wavelength
%  MESH.musm - n by 1 matrix of nodal values of the reduced scatter at
%             emission wavelength
%  MESH.kappam - n by 1 matrix of nodal values of the diffusion coefficient
%               at emission wavelength
%  MESH.muaf - n by 1 matrix of nodal values of the absorption of the
%              fluorophore
%  MESH.eta - n by 1 matrix of nodal values of the fluorescence quantum
%  yield
%  MESH.tau - n by 1 matrix of nodal values of the fluorescence lifetime
%
%  MESH.link - How to link sources and detectors as pairs ([source_number detector_number is_pair])
%
% ----'blt' MESH types specific----
%
% same as standard mesh, except we do not have a source file, and these
% need to be created for internal sources. Measurements are t all
% measurement points
%
% ##### *.node file #####
% MESH nodes coordinates.
% Each row has a following structure:
%    is_boundary x_coordinate y_coordinate z_coordinate
% Where is_boundary is 0 or 1 to indicate if a given node is the boundary
% node, x to z coordinates are cartezian coordiantes expressed in mm.
% z_coordinate is not present if MESH is a 2D MESH. is_boundary goes into
% 'MESH.bndvtx' and [x_coordinate y_coordinate z_coordinate] into
% 'MESH.nodes'.
%
% ##### *.elem file #####
% MESH nodes connectivity list - MESH elements
% Each row has a following structure:
%    vertex_1 vertex_2 vertex_3 vertex_4
% Where vertices create the MESH element. vertex_4 is not present if MESH
% is a 2D MESH as the elements are triangels. [vertex_1 vertex_2 vertex_3
% vertex_4] is saved in 'MESH.elements'. Vertices are indexes of the
% 'MESH.nodes', 1-based indexing.
%
% ##### *.param file #####
% MESH nodal parameters.
% First line specifies the 'MESH.type'. If the first row is
%    'stnd'
% the following rowns specify nodal parameters as in following format:
%    MESH.mua MESH.kappa MESH.ri
% See above for the parameters meaning and units. If the first row is
%    'spec'
% the following few lines define names of constituents, e.g.:
%    HbO
%    deoxyHb
%    Water
% However, this list is not limited as it can be user redefined. The last
% two names are always:
%    S-Amplitude
%    S-Power
% Where, S-Amplitude corresponds to the scattering amplitude and power as
% defined by the Mie theory. The list of the constituents and the scattering
% amplitude and power are saved in the 'MESH.chromscattlist'.
% Following rows specify nodal values of the spectral parameters, e.g.:
%    HbO deoxyHb Water S-Amplitude S-Power
% Where HbO, deoxyHb and Water will be saved directly into 'MESH.conc' and
% S-Amplitude, S-Power into 'MESH.sa' and 'MESH.sp' respectively.
% First line specifies the 'MESH.type'. If the first row is
%    'flour'
% the following rowns specify nodal parameters as in following format:
%    MESH.muax MESH.kappax MESH.ri MESH.muam MESH.kappam MESH.muaf MESH.eta MESH.tau
%
% ##### *.source file #####
% MESH sources definition
% If the first line says
%    fixed
% the 'MESH.source.fixed' is set to 1. Otherwise, the first line shows the
% comumn names. If the first line is 'fixed'. This will be the second line.
%    num x y z fwhm
% or
%    num x y z fwhm ele ip1 ip2 ip3 ip4
% The 'ip4' is not present for a 2D MESH. Following lines define sources,
% one per line, e.g.:
%    source_num x_coord y_coord z_coord FWHM
% or
%    source_num x_coord y_coord z_coord FWHM elem_index b_coord_1 b_coord_2 b_coord_3 b_coord_4
% The 'z_coord' and 'b_coord_4' are not present for a 2D MESH. The
% source_num will go to 'MESH.source.num', [x_coord y_coord z_coord] into
% 'MESH.source.coord', FWHM into 'MESH.source.fwhm' and [elem_index
% b_coord_1 b_coord_2 b_coord_3 b_coord_4] into 'MESH.source.int_func'. See
% above for the MESH fields meaning.
%
% ##### *.meas file #####
% MESH detectors definition
% If the first line says
%    fixed
% the 'MESH.meas.fixed' is set to 1. Otherwise, the first line shows the
% comumn names. If the first line is 'fixed'. This will be the second line.
%    num x y z
% or
%    num x y z ele ip1 ip2 ip3 ip4
% The 'ip4' is not present for a 2D MESH. Following lines define detectors,
% one per line, e.g.:
%    detector_num x_coord y_coord z_coord
% or
%    detector_num x_coord y_coord z_coord elem_index b_coord_1 b_coord_2 b_coord_3 b_coord_4
% The 'z_coord' and 'b_coord_4' are not present for a 2D MESH. The
% detector_num will go to 'MESH.meas.num', [x_coord y_coord z_coord] into
% 'MESH.meas.coord' and [elem_index b_coord_1 b_coord_2 b_coord_3 b_coord_4]
% into 'MESH.meas.int_func'. See above for the MESH fields meaning.
%
% ##### *.link file #####
% How MESH sources and detectors are organized into pairs
% The first line provides column names
%    source detector active
% Following lines define source-detector pairs, one per line, e.g. for
% 'MESH.type' equal 'stnd':
%    source_num detector_num is_active
% or 'MESH.type' equal 'spec'
%    source_num detector_num is_active@lambda1 is_active@lambda2 is_active@lambda3
% The 'source_num' is as in 'MESH.source.num', detector_num corresponds to
% 'MESH.meas.num' and is_active or is_active@lambda1 is 0 or 1 showing if a
% given pair or a given pair at a given wavelength is used. The entries
% will be stored directly into 'MESH.link'. See above for the MESH fields
% meaning.
%
% ##### *.region file #####
% MESH regions per node
% Each line defines a MESH node region. Stored directly into 'MESH.nodes'.
% See above for the MESH fields meaning.
%
% ##### *.excoef file #####
% MESH constituents spectra.
% First lines show constituents names, e.g.:
%    HbO
%    deoxyHb
%    Water
% However, this list is not limited as it can be user redefined. Following
% lines define the spectra, one wavelength per line, e.g.:
%    wavelength ex_coeff_1 ex_coeff_2 ex_coeff_3
% The 'wavelength' is saved into 'MESH.wv' and [ex_coeff_1 ex_coeff_2
% ex_coeff_3] to 'MESH.excoef'. See above for the MESH fields meaning.
%
%
% WARNING:
%   As of this NIRFAST verison support for SPN and BEM meshes has been
%   removed. The last NIRFAST version that supports SPN and BEM is NIRFAST
%   v9.1.
%
% See also SAVE_MESH, CALC_MUA_MUS
%
%   Part of NIRFAST package.

%% check in/out

narginchk(1,1);
nargoutchk(0,1);

if isstring(fn)
    fn = char(fn);
end
if ~ischar(fn)
    error('Bad mesh name. Text input expected.');
end


%% Load *.mat file or set mesh name

[~,name,ext] = fileparts(fn);

% if MATLAB *.mat file in use, just load and return
if strcmp(ext,'.mat') || (exist([name '.mat'],'file') == 2)
    load([name '.mat'], 'mesh');
    return;
end

% if ASCII files, set the mesh name
mesh.name = name;
%% Read mesh nodes

% if *.node file exist
if exist([fn '.node'],'file') == 2
    mesh.nodes = load(strcat(fn, '.node'));
    % sets 1 if boundary node, 0 if internal
    mesh.bndvtx = (mesh.nodes(:,1)); % TO_DO: change to logical
    mesh.nodes = mesh.nodes(:,2:end);
else
    error(['''' fn '.node'' file is not present']);
end

%% Read appriopriate mesh parameters

% Loading up mesh parameters from fn.param
if exist([fn '.param'],'file') == 2
    param = importdata([fn '.param']);
    % Convert from obsolete Version 1
    if ~isfield(param,'textdata')
        warning(['The mesh ''' fn ''' was saved using an obsolete version of NIRFAST. ' ...
            'Some functionalities might be limited. Saving the mesh using ''save_mesh'' will update mesh files.']);
        [~,pa] = size(param);
        % Convert old standard
        if pa == 3
            p.data = param;
            p.textdata = cellstr('stnd');
            % Convert old Fluorescence
        elseif pa == 8
            p.data = param;
            p.textdata = cellstr('fluor');
        else
            error(['''' fn '.param'' file is incorrectly formatted']);
        end
        param = p;
    end
    
    % new file format
    if isfield(param,'textdata')
        % Load standard Nirfast Mesh
        if strcmp(param.textdata(1,1),'stnd')
            mesh.type = 'stnd';
            mesh.mua = param.data(:,1);
            mesh.kappa = param.data(:,2);
            mesh.ri = param.data(:,3);
            mesh.mus = ((1./mesh.kappa)./3)-mesh.mua;
            % Load standard Fluorfast Mesh
        elseif strcmp(param.textdata(1,1),'fluor') == 1
            mesh.type = 'fluor';
            mesh.muax = param.data(:,1);
            mesh.kappax = param.data(:,2);
            mesh.musx = ((1./mesh.kappax)./3)-mesh.muax;
            mesh.ri = param.data(:,3);
            mesh.muam = param.data(:,4);
            mesh.kappam = param.data(:,5);
            mesh.musm = ((1./mesh.kappam)./3)-mesh.muam;
            mesh.muaf =  param.data(:,6);
            mesh.eta =  param.data(:,7);
            mesh.tau =  param.data(:,8);
        elseif strcmp(param.textdata(1,1),'stnd_spn')
            error('SPN meshes not supported in this NIRFAST verson. The last NIRFAST version that supports SPN is NIRFAST v9.1');
        elseif strcmp(param.textdata(1,1),'stnd_bem') || strcmp(param.textdata(1,1),'fluor_bem') || strcmp(param.textdata(1,1),'spec_bem')
            error('BEM meshes not supported in this NIRFAST verson. The last NIRFAST version that supports BEM is NIRFAST v9.1');
            % Load spectral Mesh
        elseif strcmp(param.textdata(1,1),'spec')
            mesh.type = 'spec';
            mesh.chromscattlist = param.textdata(2:end,1);
            % Get Scatter Amplitude
            ind = find(strcmpi(mesh.chromscattlist,'S-Amplitude'));
            if isempty(ind)
                error(['S-Amplitude not present in ''' fn '*.param'' file']);
            else
                mesh.sa = param.data(:,ind);
            end
            % Get Scatter Power
            ind = find(strcmpi(mesh.chromscattlist,'S-Power'));
            if isempty(ind)
                error(['S-Power not present in ''' fn '*.param'' file']);
            else
                mesh.sp = param.data(:,ind);
            end
            % Get Chromophore concentrations
            mesh.conc = [];
            k = 1;
            for i = 1 : size(mesh.chromscattlist,1)
                if ~strcmpi(mesh.chromscattlist(i),'S-Amplitude') && ~strcmpi(mesh.chromscattlist(i),'S-Power')
                    mesh.conc(:,k) = param.data(:,i);
                    k = k + 1;
                end
            end
            % Get extintion coefficient values
            if exist([fn '.excoef'],'file') == 2
                excoef = importdata([fn '.excoef']);
            else
                excoef = importdata('excoef.txt');
            end
            mesh.wv = excoef.data(:,1);
            k = 1;
            for i = 1 : size(mesh.chromscattlist,1)
                if ~strcmpi(mesh.chromscattlist(i),'S-Amplitude') && ~strcmpi(mesh.chromscattlist(i),'S-Power')
                    ind = find(strcmpi(excoef.textdata,mesh.chromscattlist(i,1)));
                    if isempty(ind)
                        error(['The Chromophore ''' char(mesh.chromscattlist(i,1)) ''' is not defined in extinction coefficient file']);
                    else
                        mesh.excoef(:,k) = excoef.data(:,ind+1);
                        k = k + 1;
                    end
                end
            end
            if ~isfield(mesh,'excoef')
                error(['Extinction coefficient not present in the ''' fn '.excoef'' file.']);
            end
            % Load standard DCS mesh
        elseif strcmp(param.textdata(1,1),'dcs')
            mesh.type = 'dcs';
            mesh.mua = param.data(:,1);
            mesh.kappa = param.data(:,2);
            mesh.ri = param.data(:,3);
            mesh.mus = ((1./mesh.kappa)./3)-mesh.mua;
            
            % DCS specific
            
            % number of colums loaded up to now
            columns_consumed = 3;
            
            % DCS wavelength
            if strncmp(param.textdata(2,1),'wv_DCS',length('wv_DCS'))
                mesh.wv_DCS = str2double(param.textdata{2,1}(length('wv_DCS')+2:end));
            else
                warning(['''dcs'' mesh type specified. However, no DCS wavelength given. '...
                    'Default wavelength of 785nm used. Change the ''MESH.wv_DCS'' accordingly.']);
                mesh.wv_DCS = 785;
            end
            
            % number of flow components (speed + concentration create flow, thus doubled)
            number_of_flows_doubled = size(param.data,2) - columns_consumed;
            
            if mod(number_of_flows_doubled,2) ~= 0
                error(['''' fn '.param'' file is incorrectly formatted. '...
                    'Flow values should be defined as pair of concentration (0 to 1) '...
                    'and brownian motion displacement coefficeint (mm^2*s^-1).']);
            end
            
            % number of flow components (speed + concentration create flow, thus divide by 2)
            number_of_flows = number_of_flows_doubled/2;
            
            % read concentrations
            mesh.alpha = param.data(:,columns_consumed+(1:number_of_flows));
            % read brownian motion displacements
            mesh.Db = param.data(:,columns_consumed+(number_of_flows+1:2*number_of_flows));
            
            mesh.aDb = sum(mesh.Db.*mesh.alpha,2);
            
        end
    end
else
    error(['''' fn '.param'' file is not present']);
end

%% Read mesh element
% if elements file exist
if exist([fn '.elem'],'file') == 2
    mesh.elements = load(strcat(fn, '.elem'));
    mesh.dimension = size(mesh.elements,2) - 1;
    % just in case as some previous NIRFAST verions saved zeroes at the z
    % coordinate for 2D meshes.
    if (mesh.dimension == 2) && (size(mesh.nodes,2) > 2)
        mesh.nodes = mesh.nodes(:,1:2);
    end
    % sort in columns (needed for the solver)
    mesh.elements = sort(mesh.elements,2);
else
    error(['''' fn '.elem'' file is not present']);
end

%% Region file

if exist([fn '.region'],'file') == 2
    mesh.region = load(strcat(fn, '.region'));
else
    mesh.region = zeros(size(mesh.nodes,1),1);
end

%% fix element orientation for 2D triangular meshes

if mesh.dimension == 2
    mesh.elements = check_element_orientation_2d(mesh.elements,mesh.nodes);
end

%% Load source locations

if exist([fn '.source'],'file') == 2
    source = importdata([fn '.source']);
    if ~isfield(source,'textdata')
        % No text at top of source file (old format)
        mesh.source.fixed = 0;
        [ns,nc] = size(source);
        mesh.source.num = (1:ns)';
        if nc == mesh.dimension
            mesh.source.fwhm = zeros(ns,1);
            mesh.source.coord = source;
        elseif nc == (mesh.dimension+1)
            mesh.source.fwhm = source(:,mesh.dimension+1);
            mesh.source.coord = source(:,1:mesh.dimension);
        else
            mesh.source.fwhm = [];
            mesh.source.coord = [];
        end
    elseif isfield(source,'textdata') && (sum(sum(strcmp(source.textdata,'num'))) == 0)
        % If text at top of source file, but columns are not labeled ('num', 'x', 'y', etc.)
        % it should only say 'fixed' at the top of the file (old format)
        mesh.source.fixed = 1;
        [ns,nc] = size(source.data);
        mesh.source.num = (1:ns)';
        if nc == mesh.dimension
            mesh.source.fwhm = zeros(ns,1);
            mesh.source.coord = source.data;
        elseif nc == mesh.dimension+1
            mesh.source.fwhm = source.data(:,mesh.dimension+1);
            mesh.source.coord = source.data(:,1:mesh.dimension);
        else
            mesh.source.fwhm = [];
            mesh.source.coord = [];
        end
    elseif isfield(source,'textdata') && (sum(sum(strcmp(source.textdata,'num'))) == 1)
        % Text flags at top of source file with column headings (new format)
        ns = size(source.data,1);
        
        mesh.source.fixed = 0;
        if sum(sum(strcmp(source.textdata,'fixed'))) == 1
            mesh.source.fixed = 1;
            source.textdata = source.textdata(2:end,:);
        end
        
        mesh.source.num = source.data(:,logical(strcmp(source.textdata,'num')));
        mesh.source.coord(:,1) = source.data(:,logical(strcmp(source.textdata,'x')));
        mesh.source.coord(:,2) = source.data(:,logical(strcmp(source.textdata,'y')));
        % only for 3D meshes
        if sum(strcmp(source.textdata,'z')) == 1
            if mesh.dimension == 2
                warning('Sources are 3D, mesh is 2D.');
            end
            mesh.source.coord(:,3) = source.data(:,logical(strcmp(source.textdata,'z')));
        end
        % spatial FWHM of sources
        if sum(strcmp(source.textdata,'fwhm')) == 1
            mesh.source.fwhm = source.data(:,logical(strcmp(source.textdata,'fwhm')));
        else
            mesh.source.fwhm = zeros(ns,1);
        end
        % what mesh element the source does belong to
        if sum(strcmp(source.textdata,'ele')) == 1
            mesh.source.int_func(:,1) = source.data(:,logical(strcmp(source.textdata,'ele')));
            % load integration functions for the source (barycentric coordinates)
            if sum(strcmp(source.textdata,'ip1')) == 1
                mesh.source.int_func(:,2) = source.data(:,logical(strcmp(source.textdata,'ip1')));
            end
            if sum(strcmp(source.textdata,'ip2')) == 1
                mesh.source.int_func(:,3) = source.data(:,logical(strcmp(source.textdata,'ip2')));
            end
            if sum(strcmp(source.textdata,'ip3')) == 1
                mesh.source.int_func(:,4) = source.data(:,logical(strcmp(source.textdata,'ip3')));
            end
            if sum(strcmp(source.textdata,'ip4')) == 1
                if mesh.dimension == 2
                    warning('Sources ''int_func'' are 3D, mesh is 2D.');
                end
                mesh.source.int_func(:,5) = source.data(:,logical(strcmp(source.textdata,'ip4')));
            end
        end
    else
        disp(['''' fn '.source'' file has a bad format.']);
    end
    
    % Check and position sources
    if logical(mesh.source.fixed) || isfield(mesh.source,'int_func')
        % if fixed or barycentric coordinates (integration functions) loaded
        if logical(mesh.source.fixed)
            disp('Fixed sources');
        end
        if isfield(mesh.source,'int_func')
            disp('Sources integration functions loaded');
        end
    else
        % if not fixed and no barycentric coordinates (integration functions) loaded
        % based on mesh type, get the scattering coefficients
        if strcmp(mesh.type,'stnd') || strcmp(mesh.type,'dcs')
            mus_eff = mesh.mus;
        elseif strcmp(mesh.type,'fluor')
            mus_eff = mesh.musx;
        elseif strcmp(mesh.type,'spec')
            [~,mus_eff] = calc_mua_mus(mesh,mesh.wv(1));
        end
        % TO_DO 'blt' ?
        disp('Moving Sources');
        [mesh] = move_source(mesh,mus_eff);
    end
    
    if isfield(mesh.source,'int_func')
        % if integration fnctions loaded or calculated while moving sources
        % check if loaded functions still valid fot the mesh
        mesh_tmp = mesh;
        % get elements where sources fall in, based on the loaded values
        mesh_tmp.elements = mesh.elements(mesh.source.int_func(:,1),:);
        % recalculate integration functions for the loaded sources coordinates
        [~,int_func] = mytsearchn(mesh_tmp,mesh.source.coord);
        % check if differences between loaded int functions and
        % calculated ones, within some precision (0.0001%)
        if any(any(abs(mesh.source.int_func(:,2:end) - int_func) > 1e-6))
            disp('Recalculating integration functions as there is some mismatch');
            % recalculate if the loaded int functions are not for the loaded source coordinates.
            [ind,int_func] = mytsearchn(mesh,mesh.source.coord);
            mesh.source.int_func = [ind int_func];
        end
    else
        % if integration fnctions are not present in the '*.source' file
        disp('Calculating sources integration functions');
        % calculate the integration functions
        [ind,int_func] = mytsearchn(mesh,mesh.source.coord);
        mesh.source.int_func = [ind int_func];
    end
    
    % check if any of the source is outside the mesh and issue a warning
    if any(~isfinite(mesh.source.int_func(:,1)))
        warning(['Source(s) outside the mesh. '...
            'Either move them manually or remove ''fixed'' from the source file. '...
            'Check sources: ' strcat(num2str(find(~isfinite(mesh.source.int_func(:,1)'))))]);
    end
else
    disp(['''' fn '.source'' file is not present']);
end

%% Load detector locations

if exist([fn '.meas'],'file') == 2
    meas = importdata([fn '.meas']);
    
    if ~isfield(meas,'textdata')
        % No text at top of meas file (old format)
        mesh.meas.fixed = 0;
        nm = size(meas,1);
        mesh.meas.num = (1:nm)';
        mesh.meas.coord = meas;
    elseif isfield(meas,'textdata')  && (sum(sum(strcmp(meas.textdata,'num'))) == 0)
        % If text at top of source file, but columns are not labeled ('num', 'x', 'y', etc.)
        % it should only say 'fixed' at the top of the file (old format)
        mesh.meas.fixed = 1;
        nm = size(meas.data,1);
        mesh.meas.num = (1:nm)';
        mesh.meas.coord = meas.data;
    elseif isfield(meas,'textdata') && (sum(sum(strcmp(meas.textdata,'num'))) == 1)
        % Text flags at top of source file with column headings (new format)
        mesh.meas.fixed = 0;
        if sum(sum(strcmp(meas.textdata,'fixed'))) == 1
            mesh.meas.fixed = 1;
            meas.textdata = meas.textdata(2,:);
        end
        
        mesh.meas.num = meas.data(:,logical(strcmp(meas.textdata,'num')));
        mesh.meas.coord(:,1) = meas.data(:,logical(strcmp(meas.textdata,'x')));
        mesh.meas.coord(:,2) = meas.data(:,logical(strcmp(meas.textdata,'y')));
        % only for 3D meshes
        if sum(strcmp(meas.textdata,'z')) == 1
            mesh.meas.coord(:,3) = meas.data(:,logical(strcmp(meas.textdata,'z')));
        end
        % what mesh element the source does belong to
        if sum(strcmp(meas.textdata,'ele')) == 1
            mesh.meas.int_func(:,1) = meas.data(:,logical(strcmp(meas.textdata,'ele')));
            % load integration functions for the source (barycentric coordinates)
            if sum(strcmp(meas.textdata,'ip1')) == 1
                mesh.meas.int_func(:,2) = meas.data(:,logical(strcmp(meas.textdata,'ip1')));
            end
            if sum(strcmp(meas.textdata,'ip2')) == 1
                mesh.meas.int_func(:,3) = meas.data(:,logical(strcmp(meas.textdata,'ip2')));
            end
            if sum(strcmp(meas.textdata,'ip3')) == 1
                mesh.meas.int_func(:,4) = meas.data(:,logical(strcmp(meas.textdata,'ip3')));
            end
            if sum(strcmp(meas.textdata,'ip4')) == 1
                if mesh.dimension == 2
                    warning('Detectors ''int_func'' are 3D, mesh is 2D.');
                end
                mesh.meas.int_func(:,5) = meas.data(:,logical(strcmp(meas.textdata,'ip4')));
            end
        end
    else
        disp(['''' fn '.meas'' file has a bad format.']);
    end
    
    % Check and/or move detectors
    if logical(mesh.meas.fixed) || isfield(mesh.meas,'int_func')
        % if fixed or barycentric coordinates (integration functions) loaded
        if logical(mesh.meas.fixed)
            disp('Fixed detectors');
        end
        if isfield(mesh.meas,'int_func')
            disp('Detectors integration functions loaded');
        end
    else
        % if not fixed and no barycentric coordinates (integration functions) loaded
        disp('Moving detectors');
        [mesh] = move_detector(mesh);
    end
    
    if isfield(mesh.meas,'int_func') == 1
        % if integration fnctions loaded or calculated while moving detectors
        % check if loaded functions still valid fot the mesh
        mesh_tmp = mesh;
        % get elements where detectors fall in, based on the loaded values
        mesh_tmp.elements = mesh.elements(mesh.meas.int_func(:,1),:);
        % recalculate integration functions for the loaded detectors coordinates
        [~,int_func] = mytsearchn(mesh_tmp,mesh.meas.coord);
        % check if differences between loaded int functions and
        % calculated ones, within some precision (0.0001%)
        if any(any(abs(mesh.meas.int_func(:,2:end) - int_func) > 1e-6))
            disp('Recalculating integration functions as there is some mismatch');
            % recalculate if the loaded int functions are not for the loaded detector coordinates.
            [ind,int_func] = mytsearchn(mesh,mesh.meas.coord);
            mesh.meas.int_func = [ind int_func];
        end
    else
        % if integration fnctions are not present in the '*.meas' file
        disp('Calculating detectors integration functions');
        % calculate the integration functions
        [ind,int_func] = mytsearchn(mesh,mesh.meas.coord);
        mesh.meas.int_func = [ind int_func];
    end
    
    % check if any of the detector is outside the mesh and issue a warning
    if any(~isfinite(mesh.meas.int_func(:,1)))
        warning(['Detector(s) outside the mesh. '...
            'Either move them manually or remove ''fixed'' from the detector file. '...
            'Check detectors: ' strcat(num2str(find(~isfinite(mesh.meas.int_func(:,1)'))))]);
    end
else
    disp(['''' fn '.meas'' file is not present']);
end

%% Load link list for source and detector

if exist([fn '.link'],'file') == 2
    % determine if link file is legacy format
    fid = fopen([fn '.link']);
    first_line = fgetl(fid);
    fclose(fid);
    if ~strcmp(first_line(1),'s')
        % convert to new
        link = load([fn '.link']);
        [n,m] = size(link);
        for i = 1:n
            for j = 1:m
                if link(i,j) ~= 0
                    if strcmp(mesh.type,'spec')
                        mesh.link(i,:) = [i link(i,j) ones(1,length(mesh.wv))];
                    else
                        mesh.link(i,:) = [i link(i,j), 1];
                    end
                else
                    if strcmp(mesh.type,'spec')
                        mesh.link(i,:) = [i link(i,j) zeros(1,length(mesh.wv))];
                    else
                        mesh.link(i,:) = [i link(i,j), 0];
                    end
                end
            end
        end
    else
        link = importdata([fn '.link']);
        mesh.link = link.data;
        if isfield(mesh,'wv')
            if size(mesh.link,2) ~= length(mesh.wv)+2
                mesh.link = [mesh.link(:,1:2) ones(size(mesh.link,1),length(mesh.wv))];
                warning(['''mesh.wv'' length (number of wavelengths) does not match the '''...
                    fn '.link'' file. ''mesh.link'' changed to match ''mesh.wv''']);
            end
        end
    end
else
    disp(['''' fn '.link'' file is not present']);
end

%% speed of light in medium
% If a spectral mesh, assume Refractive index = 1.33
if strcmp(mesh.type,'spec')
    mesh.ri = ones(size(mesh.nodes,1),1).*1.33;
end
% speed of light
c0 = 299792458000; % / mm/s
mesh.c = c0./mesh.ri;

%% Set boundary coefficient using definition of baundary attenuation A using the Fresenel's law:

% outside mesh refractive index
n_air = 1;
% Robin boundary condition
% attenuation on the mesh-outside_mesh boundary
A = boundary_attenuation(mesh.ri,n_air,'robin');
% scale factor on the mesh-outside_mesh boundary
mesh.ksi = 1./(2*A);

%% area and support of each element
% mesh elements area (2D) or volume (3D)
mesh.element_area = ele_area_c(mesh.nodes(:,1:mesh.dimension),mesh.elements);
% total area (2D) or volume (3D) of elements connected to each node. (how much support is around a node)
mesh.support = mesh_support(mesh.nodes(:,1:mesh.dimension), mesh.elements, mesh.element_area);

end
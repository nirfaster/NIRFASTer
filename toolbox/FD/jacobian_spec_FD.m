function [J,data] = jacobian_spec_FD(mesh,varargin)
% JACOBIAN_SPEC_FD Calculates spatial distributions of sensitivity of
%   field registerd on the mesh boundary to changes of chromophores
%   concentrations and Mie scattering amplitudea nd power per mesh node.
%   This function is specialized to be called for NIRFAST meshes of 'spec'
%   type (multi-wavelength).It deliveres spectrally-resolved distributions
%   for all chromophores (e.g. haemoglobins (HbO2 and Hb), water, etc.) and
%   scattering amplitude and power. E.g.: 
%   If data are real (CW - continuous wave):
%    - J_{CHbO2}^{A} - expressesd in mM^-1 (one over milli molar (molar -
%                    moles/liter)), sensitivity of boundary attenuation (A)
%                    to oxygenated haemoglobin concentration (CHbO2 [mM]).
%    - J_{CHb}^{A}   - expressesd in mM^-1 (one over milli molar (molar -
%                    moles/liter)), sensitivity of boundary attenuation (A)
%                    to deoxygenated haemoglobin concentration (CHb [mM]).
%    - J_{W}^{A}     - unitless, sensitivity of boundary attenuation (A) to
%                    water volume fraction (W between 0 and 1). 
%    - J_{sa}^{A}    - (OPTIONAL) unitless, sensitivity of boundary
%                    attenuation (A) to Mie scatterin amplitude (sa,
%                    unitless). Returned if the optional parameter is set 
%                    to 'all'.
%    - J_{sa}^{A}    - (OPTIONAL) unitless, sensitivity of boundary
%                    attenuation (A) to Mie scatterin power (sp, unitless).
%                    Returned if the optional parameter is set to 'all'.
%   If data are complex (FD - frequency domain):
%    - J_{CHbO2}^{A}, J_{CHb}^{A}, J_{W}^{A}, J_{sa}^{A} and J_{sp}^{A} are
%                   returned af for CW. Howeverm J_{sa}^{A} and J_{sp}^{A}
%                   are always returned.
%    - J_{CHbO2}^{phi}, J_{CHb}^{phi}, J_{W}^{phi}, J_{sa}^{phi} and
%                   J_{sp}^{phi} are returned where 'phi' is the phase
%                   shift in radians. Jacobians hacve the same units as for
%                   in CW case.
% 
% SYNTAX:
%  [J,DATA] = JACOBIAN_SPEC_FD(MESH)
%  [J,DATA] = JACOBIAN_SPEC_FD(MESH,FREQUENCY)
%  [J,DATA] = JACOBIAN_SPEC_FD(MESH,FREQUENCY,SECOND_MESH_BASIS)
%  [J,DATA] = JACOBIAN_SPEC_FD(MESH,FREQUENCY,SECOND_MESH_BASIS,SOLVER)
%  [J,DATA] = JACOBIAN_SPEC_FD(MESH,FREQUENCY,SECOND_MESH_BASIS,SOLVER,OPTIONS)
%  [J,DATA] = JACOBIAN_SPEC_FD(MESH,FREQUENCY,SECOND_MESH_BASIS,SOLVER,OPTIONS,'all') 
%  [J,DATA] = JACOBIAN_SPEC_FD(MESH,FREQUENCY,[],[],[],[])
%  [J] = JACOBIAN_SPEC_FD(___)
%  
% [J,DATA] = JACOBIAN_SPEC_FD(MESH) MESH is NIRFAST mesh structure where
%   'MESH.type' is set to 'spec'. J is the jacobian for CW (continuous
%   wave) type sources. J is a matrix of size P*L by C*N where P is
%   number of active source-detectors pairs as enabled in 'MESH.link', L is
%   number of active wavelengths as enabled in the 'MESH.link', C is number
%   of chromophores in the MESH structure (columns od 'MESH.conc') and N is
%   number ot the MESH nodes. Returned J is a matrix build in the following
%   mannere for wavelengths lam1, lam2 to lamL:
%    |J_{CHbO2}^{A}@lam1 J_{CHb}^{A}@lam1 J_{W}^{A}@lam1|
%    |J_{CHbO2}^{A}@lam2 J_{CHb}^{A}@lam2 J_{W}^{A}@lam2|
%    |                         ...                      |
%    |J_{CHbO2}^{A}@lamL J_{CHb}^{A}@lamL J_{W}^{A}@lamL|
%   Where each sub-matrix has P by N size.
%   DATA is the structure as returned by 'femdata_FD' with additional field
%   'aphi'. 'DATA.phi' is the forward field 'DATA.aphi' represents the
%   adjoint field (detectors are new sources) and the boundary data are
%   present in 'DATA.complex'.
%    - DATA.PHI - (direct field) N by M by L matrix of photon fluence rate
%       at N nodes for M sources at L wavelengths as used for
%       source-detector pairs enabled in 'MESH.link'. Expressed in
%       mm^-2s^-1 (3D) and mm^-1s^-1 (2D). 
%    - DATA.APHI - (adjoint field) N by K by L matrix of photon fluence
%       rate at N nodes for K detectors used as sources and at L
%       wavelengths as used for source-detector pairs enabled in
%       'MESH.link'. Detectors as used for source-detector pairs enabled in
%       'MESH.link' field. Expressed in mm^-2s^-1 (3D) and mm^-1s^-1 (2D).
%    - DATA.COMPLEX - (boundary data) P by W matrix of photon fluence
%       rate for P source-detectors pairs and W wavelength as specified in
%       'MESH.link'. Expressed in mm^-2s^-1 (3D) and mm^-1s^-1 (2D). Values
%       for turned off pairs and wavelengths are set to NaN.
% 
% [J,DATA] = JACOBIAN_SPEC_FD(MESH,FREQUENCY) Optional FREQUENCY is
%   modulation frequency of sources in Hz. J is a matrix of size 2*P*L by
%   (C+S)*N where P is number of active source-detectors pairs as enabled
%   in 'MESH.link', L is number of active wavelengths as enabled in the
%   'MESH.link', C is number of chromophores in the MESH structure (columns
%   od 'MESH.conc'), S is equal 2 (scattering amplitude and power) and N is 
%   number ot the MESH nodes. Returned J is a matrix build in the following
%   mannere for wavelengths lam1, lam2 to lamL:
%    |J_{CHbO2}@lam1 J_{CHb}@lam1 J_{W}@lam1 J_{sa}@lam1 J_{sp}@lam1|
%    |J_{CHbO2}@lam2 J_{CHb}@lam2 J_{W}@lam2 J_{sa}@lam2 J_{sp}@lam2|
%    |                           ...                                |
%    |J_{CHbO2}@lamL J_{CHb}@lamL J_{W}@lamL J_{sa}@lamL J_{sp}@lamL|
%   Where each sub-matrix has 2*P by N size. E.g:
%       J_{CHbO2}@lam1(1:2:end,:) = J_{CHbO2}^{A};
%       J_{CHbO2}@lam1(2:2:end,:) = J_{CHbO2}^{phi};
% 
% [J,DATA] = JACOBIAN_SPEC_FD(MESH,FREQUENCY,SECOND_MESH_BASIS)
%   Optional SECOND_MESH_BASIS is a mesh in NIRFAST format that will be
%   used to calculate the Jacobians. This mesh is usually a coarse mesh as
%   compared with the MESH. See e.g. 'pixel_basis' or 'second_mesh_basis'
%   for details how to create such mesh.
% 
% [J,DATA] = JACOBIAN_SPEC_FD(MESH,FREQUENCY,SECOND_MESH_BASIS,SOLVER) Optional,
%   user specified SOLVER. Type 'help get_solver' for how to use solvers.
% 
% [J,DATA] = JACOBIAN_SPEC_FD(MESH,FREQUENCY,SECOND_MESH_BASIS,SOLVER,OPTIONS)
%   Optional OPTIONS is a structure that allows to control solvers
%   parameters. Default OPTIONS structure is returned by 'solver_options'
%   function. See 'help solver_options' on how to use this structure.
% 
% [J,DATA] = JACOBIAN_SPEC_FD(MESH,FREQUENCY,SECOND_MESH_BASIS,SOLVER,OPTIONS,'all') 
%   The 'all' option has effect only if data are real (CW - continuous
%   wave). If specified for the CW data, the attenuation sensitivity to
%   scattering amplitude (sa) and power (sp) are returned as well. In this
%   mode, J is a matrix of size P*L by (C+S)*N where P is number of active
%   source-detectors pairs as enabled in 'MESH.link', L is number of active
%   wavelengths as enabled in the 'MESH.link', C is number of chromophores
%   in the MESH structure (columns od 'MESH.conc'), S is equal 2
%   (scattering amplitude and power) and N is number ot the MESH nodes.
%   Returned J is a matrix build in the following mannere for wavelengths
%   lam1, lam2 to lamL: 
%    |J_{CHbO2}^{A}@lam1 J_{CHb}^{A}@lam1 J_{W}^{A}@lam1 J_{sa}^{A}@lam1 J_{sp}^{A}@lam1|
%    |J_{CHbO2}^{A}@lam2 J_{CHb}^{A}@lam2 J_{W}^{A}@lam2 J_{sa}^{A}@lam2 J_{sp}^{A}@lam2|
%    |                                          ...                                     |
%    |J_{CHbO2}^{A}@lamL J_{CHb}^{A}@lamL J_{W}^{A}@lamL J_{sa}^{A}@lamL J_{sp}^{A}@lamL|
%   Where each sub-matrix has P by N size.
% 
% [J,DATA] = JACOBIAN_SPEC_FD(MESH,FREQUENCY,[],[],[],'all') Any of the
%   parameters: SECOND_MESH_BASIS, SOLVER and OPTIONS can be set to empty.
%   This will set a given parameter to its default value or no 'coarse'
%   mesh is used.
% 
% [J] = JACOBIAN_SPEC_FD(___) Does not retur the forward data. Jacobians only.
% 
% See also JACOBIAN_FD, JACOBIAN_STND_FD, BUILD_JACOBIAN_FD.
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz and H. Dehghani 2018


%% check in/out

narginchk(1,6);
nargoutchk(0,2);

%% get optional inputs, handle variable input

% default
frequency = 0;
second_mesh_basis = [];
solver = get_solver;
OPTIONS = solver_options;
isAll = false;

if ~isempty(varargin)
    if length(varargin) >= 1
        % frequency
        if ~ischar(varargin{1}) && ~isstring(varargin{1}) && (numel(varargin{1})==1)
            % sanity check
            if varargin{1} < 0
                error('Negative frequency value. Please see help for details on how to use this function.')
            else
                frequency = varargin{1};
            end
        else
            error('Bad 2nd argument value. A scalar expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) >= 2
        % second mesh basis
        if isstruct(varargin{2}) || isempty(varargin{2})
            second_mesh_basis = varargin{2};
        else
            error('Bad 3nd argument value. A mesh structure expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) >= 3
        % solver
        if ischar(varargin{3}) || isstring(varargin{3})
            % user specified, sanity check
            solver = get_solver(varargin{3});
        elseif isstruct(varargin{3})
            OPTIONS = varargin{3};
        elseif isempty(varargin{3})
            solver = get_solver;
        else
            error('Bad 4th argument value. Solver name or solver settings structure expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) >= 4
        % solver options
        if isstruct(varargin{4})
            OPTIONS = varargin{4};
        elseif isempty(varargin{4})
            OPTIONS = solver_options;
        else
            error('Bad 5th argument value. Solver settings structure expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) == 5
        % if scattering for CW data as well
        if ischar(varargin{5}) || isstring(varargin{5})
            if strcmp(varargin{5},'all')
                isAll = true;
            end
        else
            error('Bad 6th argument value. Text expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) > 5
        error('Bad arguments. Please see the help for details on how to use this function.')
    end
end

%% If not a workspace variable, load mesh
if ~isstruct(mesh)
    mesh = load_mesh(mesh);
end

%% Calculate field and boundary data for all sources

data = femdata_spec_FD(mesh, frequency, solver, OPTIONS);

%% Now calculate adjoint source vector (detectors are new sources)

% make sources copy and swap with detectors
sources_copy = mesh.source;
detectors_copy = mesh.meas;
link_copy = mesh.link;

% swap source
mesh.source.fixed = 0;
mesh.source.num = detectors_copy.num;
mesh.source.coord = detectors_copy.coord;
mesh.source.fwhm = zeros(size(detectors_copy.num));
mesh.source.int_func = detectors_copy.int_func;
%swap detector
mesh.meas.fixed = 0;
mesh.meas.num = sources_copy.num;
mesh.meas.coord = sources_copy.coord;
%swap link
mesh.link(:,1) = link_copy(:,2);
mesh.link(:,2) = link_copy(:,1);
% a missing detectors 'mesh.meas.int_func' will calculate itself in 'get_boundary_data' if needed

% calculate adjoint field for all sources (detectors)
data_detector = femdata_spec_FD(mesh, frequency, solver, OPTIONS);
data.aphi = data_detector.phi;

% swap back sources, detectors and link
mesh.source = sources_copy;
mesh.meas = detectors_copy;
mesh.link = link_copy;


%% LOOP THROUGH WAVELENGTHS

% number of chromophores
chromophores_num = size(mesh.excoef,2);
% number of nodes
if ~isempty(second_mesh_basis)
    nodes_num = size(second_mesh_basis.nodes,1);
else
    nodes_num = size(mesh.nodes,1);
end

% number of scattering phase and amplitude jacobians per wv
scattering_jacobians = 0;

% number of active measurements
measurements_active = sum(isfinite(data.complex(:,1)));
if (~isreal(data.phi) && ~isreal(data.aphi))
    measurements_active = measurements_active * 2; % amplitude and phase
    scattering_jacobians = 2; % number scattering phase and amplitude jacobians per wv
end

% check if a wavelength is turned off for all sources.
mask_wv_off = sum(mesh.link(:,3:end),1) == 0;
% number of turned on wavelengths
wavelength_on_num = sum(~mask_wv_off);

% if scattering for CW measurements as well
if isAll
    % append scattering phase and amplitude to chromophores
    scattering_jacobians = 2; % number scattering phase and amplitude jacobians per wv
end

% initialize spectral Jacobian
J = zeros(measurements_active * wavelength_on_num, (chromophores_num + scattering_jacobians) * nodes_num);

for ind_wv = 1:length(mesh.wv)
    % calculate only if wavelength requested.
    if ~mask_wv_off(ind_wv)
        
        %%% spatially convolve fileds for sources and detectors

        % copy single wavelength fields
        data_single_wv.phi = squeeze(data.phi(:,:,ind_wv));
        data_single_wv.aphi = squeeze(data.aphi(:,:,ind_wv));
        % jacobian build requires data for enabled pairs only
        data_single_wv.complex = data.complex(logical(data.link(:,3)),ind_wv);

        % Calculate Jacobian
        if ~isempty(second_mesh_basis) % use second mesh basis for Jacobian
            % make a copy of spectral link
            link_copy = second_mesh_basis.link;
            % put a single wavelength link into the mesh
            second_mesh_basis.link = second_mesh_basis.link(:,[1 2 ind_wv+2]);

            % get optical properties at nodes using spectra
            [second_mesh_basis.mua, second_mesh_basis.mus, second_mesh_basis.kappa, Ex_coefs] = calc_mua_mus(second_mesh_basis, second_mesh_basis.wv(ind_wv));
            
            % if sources are not fixed, move sources depending on mus
            if ~logical(second_mesh_basis.source.fixed)
                second_mesh_basis = move_source(second_mesh_basis,second_mesh_basis.mus);
            end
            
            % build jacobian on data interpoladed onto the coarse mesh
            if isAll
                [J_single_wv] = build_jacobian_FD(second_mesh_basis,interpolatef2r(mesh,second_mesh_basis,data_single_wv),'all');
            else
                [J_single_wv] = build_jacobian_FD(second_mesh_basis,interpolatef2r(mesh,second_mesh_basis,data_single_wv));
            end
            
            % put back the spectral link into the mesh
            second_mesh_basis.link = link_copy;
            
            % see Dehghani2008 (page 7) for the factors meaning (in short, they convert the kappa jacobian into Mie scattering amplitude and power jacobians)
            sa_factor = -3*(second_mesh_basis.kappa.^2).*((second_mesh_basis.wv(ind_wv)/1000).^(-second_mesh_basis.sp));
            sp_factor = (-3*(second_mesh_basis.kappa.^2).*second_mesh_basis.mus).*(-log(second_mesh_basis.wv(ind_wv)/1000));
        else
            % make a copy of spectral link
            link_copy = mesh.link;
            % put a single wavelength link into the mesh
            mesh.link = mesh.link(:,[1 2 ind_wv+2]);

            % get optical properties at nodes using spectra
            [mesh.mua, mesh.mus, mesh.kappa, Ex_coefs] = calc_mua_mus(mesh, mesh.wv(ind_wv));
            
            % if sources are not fixed, move sources depending on mus
            if ~logical(mesh.source.fixed)
                mesh = move_source(mesh,mesh.mus);
            end
            
            % build jacobian on data interpoladed onto the coarse mesh
            if isAll
                [J_single_wv] = build_jacobian_FD(mesh,data_single_wv,'all');
            else
                [J_single_wv] = build_jacobian_FD(mesh,data_single_wv);
            end
            
            % put back the spectral link into the mesh
            mesh.link = link_copy;
            
            % see Dehghani2008 (page 7) for the factors meaning (in short, they convert the kappa jacobian into Mie scattering amplitude and power jacobians)
            sa_factor = -3*(mesh.kappa.^2).*((mesh.wv(ind_wv)/1000).^(-mesh.sp));
            sp_factor = (-3*(mesh.kappa.^2).*mesh.mus).*(-log(mesh.wv(ind_wv)/1000));
        end
        
        
        
        %%% Assign outputs        
        if size(J_single_wv.complete,2) == (2*nodes_num)
            % jacobians for chromophores
            for ind_chromo = 1:chromophores_num
                J((ind_wv-1)*measurements_active+1:(ind_wv)*measurements_active,(ind_chromo-1)*nodes_num+1:(ind_chromo)*nodes_num) =...
                    J_single_wv.complete(:,nodes_num+1:end)*Ex_coefs(ind_chromo);
            end
            % increment column index to save scattering amplitude
            ind_chromo = ind_chromo + 1;
            J((ind_wv-1)*measurements_active+1:(ind_wv)*measurements_active,(ind_chromo-1)*nodes_num+1:(ind_chromo)*nodes_num) =...
                J_single_wv.complete(:,1:nodes_num).*sa_factor';
            % increment column index to save scattering power
            ind_chromo = ind_chromo + 1;
            J((ind_wv-1)*measurements_active+1:(ind_wv)*measurements_active,(ind_chromo-1)*nodes_num+1:(ind_chromo)*nodes_num) =...
                J_single_wv.complete(:,1:nodes_num).*sp_factor';
        else
            % jacobians for chromophores
            for ind_chromo = 1:chromophores_num
                J((ind_wv-1)*measurements_active+1:(ind_wv)*measurements_active,(ind_chromo-1)*nodes_num+1:(ind_chromo)*nodes_num) =...
                    J_single_wv.complete(:,1:nodes_num)*Ex_coefs(ind_chromo);
            end
        end
    end
end

% remove fields
data = rmfield(data,'aphi');
data = rmfield(data,'phi');

end

function [data_recon] = interpolatef2r(fwd_mesh,recon_mesh,data_fwd)
% This function interpolates fwd_mesh data onto recon_mesh
% Used to calculate the Jacobian on second mesh

% copy boundary data
data_recon.complex = data_fwd.complex;

% calculate interpolation functions if needed
if ~isfield(fwd_mesh,'fine2coarse')
    fwd_mesh.fine2coarse = second_mesh_basis(fwd_mesh,recon_mesh);
end

% loop through nodes of recon_mesh
for i = 1:size(recon_mesh.nodes,1)
    if fwd_mesh.fine2coarse(i,1) > 0
        data_recon.phi(i,:) = (fwd_mesh.fine2coarse(i,2:end) * ...
            data_fwd.phi(fwd_mesh.elements(fwd_mesh.fine2coarse(i,1),:),:));
        data_recon.aphi(i,:) = (fwd_mesh.fine2coarse(i,2:end) * ...
            data_fwd.aphi(fwd_mesh.elements(fwd_mesh.fine2coarse(i,1),:),:));
        %         else
        %             error(['The interpolated mesh (fwd_mesh) has bad endtry in the integration functions ''fine2coarse''. ' ...
        %                 'Check row: ' num2str(i) '.'])
    end
end

end

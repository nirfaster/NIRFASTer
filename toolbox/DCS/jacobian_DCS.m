function [J,data] = jacobian_DCS(mesh,tau_DCS,varargin)
% JACOBIAN_DCS Calculates spatial distributions of sensitivity of
%   normalized autocorrelation registerd on the mesh boundary (detector) to
%   changes of flow per mesh node.
%
% SYNTAX:
%  [J,DATA] = JACOBIAN_DCS(MESH,TAU_DCS)
%  [J,DATA] = JACOBIAN_DCS(MESH,TAU_DCS,WV_DCS)
%  [J,DATA] = JACOBIAN_DCS(MESH,TAU_DCS,WV_DCS,SECOND_MESH_BASIS)
%  [J,DATA] = JACOBIAN_DCS(MESH,TAU_DCS,WV_DCS,SECOND_MESH_BASIS,SOLVER)
%  [J,DATA] = JACOBIAN_DCS(MESH,TAU_DCS,WV_DCS,SECOND_MESH_BASIS,SOLVER,OPTIONS)
%  [J,DATA] = JACOBIAN_DCS(MESH,TAU_DCS,[],[],[],[])
%  [J] = JACOBIAN_DCS(___)
% 
% [J,DATA] = JACOBIAN_DCS(MESH,TAU_DCS) Calculates J (Jacobian), the
%   sensitivity of normalized autocorrelation registerd on the mesh
%   boundary (detector) to changes of flow per mesh node. It is required
%   the MESH have following fields:
%    - MESH.aDb - nodal values of the 'flow index' expressed in mm^2*s^-1.
%    - MESH.wv_DCS - a scalar, the DCS wavelength in nm, e.g. 785
%   TAU_DCS expressed in seconds, is a vector of k requested
%    autocorrelation delay time values. 
%   J is n by p by k matrix of the sensitivities expressed in s*mm^-2,
%    where n is number of MESH nodes, p is numer of active source-detector
%    pairs as specified in 'MESH.link' and k is number of requested
%    autocorrelation times TAU_DCS.
%   DATA structure hods the boundary data and has the following entries:
%    - DATA.amplitude - p elements vector of photon fluence rate for p
%       source-detectors pairs as specified in MESH.link. Expressed in
%       mm^-2s^-1. The electric field temporal autocorrelation function for
%       the autocorrelation delay time equal zero.
%    - DATA.link - copy of MESH.link showing how sources and detectors are
%       linked together into pairs.
%    - DATA.tau_DCS - copy of the TAU_DCS, vector of k requested
%       autocorrelation delay times in seconds.
%    - DATA.G1_DCS - p by k matrix of the electric field temporal
%       autocorrelation function at the mesh boundary (at detectors) for p
%       source-detctor pairs as specified in MESH.link and sampled at k
%       autocorrelation delay times as requested in TAU_DCS. Expressed in
%       mm^-2s^-1. 
% 
% [J,DATA] = JACOBIAN_DCS(MESH,TAU_DCS,WV_DCS) WV_DCS is an override of
%   'MESH.wv_DCS'. It specifies the DCS long-coherence light source
%   wavelengths. Expressed in nm.
% 
% [J,DATA] = JACOBIAN_DCS(MESH,TAU_DCS,WV_DCS,SECOND_MESH_BASIS)
%   Optional SECOND_MESH_BASIS is a mesh in NIRFAST format that will be
%   used to calculate the Jacobians. This mesh is usually a coarse mesh as
%   compared with the MESH. See e.g. 'pixel_basis' or 'second_mesh_basis'
%   for details how to create such mesh.
% 
% [J,DATA] = JACOBIAN_DCS(MESH,TAU_DCS,WV_DCS,SECOND_MESH_BASIS,SOLVER)
%   User specified SOLVER is used. See help of the 'get_solvers' function
%   for how to use solvers. 
% 
% [J,DATA] = JACOBIAN_DCS(MESH,TAU_DCS,WV_DCS,SECOND_MESH_BASIS,SOLVER,OPTIONS)
%   Optional OPTIONS is a structure that allows to control solvers
%   parameters. Default OPTIONS structure is returned by 'solver_options'
%   function. See 'help solver_options' on how to use this structure.
% 
% [J,DATA] = JACOBIAN_DCS(MESH,TAU_DCS,[],[],[],[]) Any of the parameters:
%   WV_DCS, SECOND_MESH_BASIS, SOLVER and OPTIONS can be set to empty. This
%   will set a given parameter to its default value or no 'coarse' mesh is
%   used.
% 
% [J] = JACOBIAN_DCS(___) Does not retur the forward data. Jacobians only.
% 
% See also FEMDATA_DCS, BUILD_JACOBIAN_FD.
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

%% check in/out

narginchk(2,6);
nargoutchk(0,2);

%% get optional inputs, handle variable input

% default
wv_DCS = [];
second_mesh_basis = [];
solver = get_solver;
OPTIONS = solver_options;

if ~isempty(varargin)
    if length(varargin) >= 1
        if ~ischar(varargin{1}) && ~isstring(varargin{1}) && (numel(varargin{1})==1)
            % sanity check
            if varargin{1} < 0
                error('Negative wavelength value. Please see help for details on how to use this function.')
            else
                wv_DCS = varargin{1};
            end
        elseif ~isempty(varargin{1})
            warning('Bad 3rd argument value. A scalar or empty value expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) >= 2
        % second mesh basis
        if isstruct(varargin{2}) || isempty(varargin{2})
            second_mesh_basis = varargin{2};
        else
            error('Bad 4th argument value. A mesh structure expected. Please see help for details on how to use this function.')
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
            error('Bad 5th argument value. Solver name or solver settings structure expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) == 4
        % solver options
        if isstruct(varargin{4})
            OPTIONS = varargin{4};
        elseif isempty(varargin{4})
            OPTIONS = solver_options;
        else
            error('Bad 6th argument value. Solver settings structure expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) > 4
        error('Bad arguments. Please see the help for details on how to use this function.')
    end
end

%% If not a workspace variable, load mesh
if ~isstruct(mesh)
    mesh = load_mesh(mesh);
end

%% get the DCS wavelength from the mesh if not specified as the function argument
if isempty(wv_DCS)
    if ~isfield(mesh,'wv_DCS')
        error('The ''wv_DCS'' field is missing. Check your mesh structure.');
    end
    wv_DCS = mesh.wv_DCS;
end

%% Calculate field and boundary data for all sources

data = femdata_DCS(mesh, tau_DCS, wv_DCS, solver, OPTIONS);

%% Now calculate adjoint field (detectors are new sources)

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
data_detector = femdata_DCS(mesh, tau_DCS, wv_DCS, solver, OPTIONS);

% swap back sources, detectors and link
mesh.source = sources_copy;
mesh.meas = detectors_copy;
mesh.link = link_copy;

%% build the Jacobian

% copy data with possible NaN values for disabled pairs
boundary_data_copy = data.G1_DCS;

% jacobian build requires data for enabled pairs only
data.G1_DCS = data.G1_DCS(logical(data.link(:,3)),:);

% declare the output Jacobian
if ~isempty(second_mesh_basis) % use second mesh basis for Jacobian
    J = zeros(size(second_mesh_basis.nodes,1),size(data.G1_DCS,1),numel(tau_DCS));
    second_mesh_basis = interpolate_mesh2mesh(mesh,second_mesh_basis,'ri','mua','mus');
    % wave number
    k0 = 2*pi*second_mesh_basis.ri./(wv_DCS * 1e-6); % <- DCS wavelength converted from nm to mm
else
    J = zeros(size(mesh.nodes,1),size(data.G1_DCS,1),numel(tau_DCS));
    % wave number
    k0 = 2*pi*mesh.ri./(wv_DCS * 1e-6); % <- DCS wavelength converted from nm to mm
end

% loop through the autocorrelation delay times
for ind_tau = 1:numel(tau_DCS)
   
    % copy DCS data to match the FD structure to simply use the FD jacobian claculation procedure
    data.phi = squeeze(data.phi_DCS(:,:,ind_tau));
    data.aphi = squeeze(data_detector.phi_DCS(:,:,ind_tau));
    data.complex = data.G1_DCS(:,ind_tau);
    % Calculate Jacobian
    if ~isempty(second_mesh_basis) % use second mesh basis for Jacobian
        % build jacobian on data interpoladed onto the coarse mesh
        J_tau = build_jacobian_FD(second_mesh_basis,interpolatef2r(mesh,second_mesh_basis,data));
    else
        % Calculate Jacobian
        J_tau = build_jacobian_FD(mesh,data);
    end
    
    % assigne the current tau Jacobian
    J(:,:,ind_tau) = J_tau.complete';
    
    
    % make DCS jacobian for aDb (scale as required)
    if ~isempty(second_mesh_basis) % use second mesh basis for Jacobian
        J(:,:,ind_tau) = J(:,:,ind_tau) * 6 .* (k0).^2 .* (second_mesh_basis.mus + second_mesh_basis.mua).^2 * tau_DCS(ind_tau);
    else
        J(:,:,ind_tau) = J(:,:,ind_tau) * 6 .* (k0).^2 .* (mesh.mus + mesh.mua).^2 * tau_DCS(ind_tau);
    end
    
end

% if data output requestes, format the data
if nargout == 2
    data = rmfield(data,{'phi','phi_DCS','aphi','complex'});
    data.G1_DCS = boundary_data_copy;
end


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
        else
            error(['The interpolated mesh (fwd_mesh) has bad endtry in the integration functions ''fine2coarse''. ' ...
                'Check row: ' num2str(i) '.'])
        end
    end

end

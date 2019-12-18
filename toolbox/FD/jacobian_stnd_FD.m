function [J,data] = jacobian_stnd_FD(mesh,varargin)
% JACOBIAN_STND_FD Calculates spatial distributions of sensitivity of
%   field registerd on the mesh boundary to changes of optical properties
%   per mesh node. This function is specialized to be called for NIRFAST
%   meshes of 'stnd' type (single wavelength).It deliveres up to 4
%   distributions:
%   If data are real (CW - continuous wave):
%    - J_{mua}^{A}     - expressesd in mm, sensitivity of boundary
%                      attenuation (A) to absorption (mua [mm^-1]). This is
%                      the mean partial pathlength photons travel through
%                      the volume unit. The sum of this Jacobian gives mean
%                      photons pathlength in mm. 
%    - J_{kappa}^{A}   - (OPTIONAL) expressed in mm^-1, sensitivity of
%                      boundary attenuation (A) to diffusion coefficient
%                      'kappa' ((3*(mua+mus))^-1 [mm]). Returned if the
%                      optional parameter is set to 'all'.
%   If data are complex (FD - frequency domain):
%    - J_{mua}^{A}     - af for CW.
%    - J_{kappa}^{A}   - af for CW but always returned.
%    - J_{mua}^{phi}   -  expressed in mm, sensitivity of boundary phase
%                      shift (phi), as related to the source phase, to
%                      absorption (mua [mm]).
%    - J_{kappa}^{phi} - expressed in mm^-1, sensitivity of boundary phase
%                      shift (phi) to diffusion coefficient 'kappa'
%                      ((3*(mua+mus))^-1 [mm]).
%  
% SYNTAX:
%  [J,DATA] = JACOBIAN_STND_FD(MESH)
%  [J,DATA] = JACOBIAN_STND_FD(MESH,FREQUENCY)
%  [J,DATA] = JACOBIAN_STND_FD(MESH,FREQUENCY,SECOND_MESH_BASIS)
%  [J,DATA] = JACOBIAN_STND_FD(MESH,FREQUENCY,SECOND_MESH_BASIS,SOLVER)
%  [J,DATA] = JACOBIAN_STND_FD(MESH,FREQUENCY,SECOND_MESH_BASIS,SOLVER,OPTIONS)
%  [J,DATA] = JACOBIAN_STND_FD(MESH,FREQUENCY,SECOND_MESH_BASIS,SOLVER,OPTIONS,'all') 
%  [J,DATA] = JACOBIAN_STND_FD(MESH,FREQUENCY,[],[],[],[])
%  [J] = JACOBIAN_STND_FD(___)
% 
% [J,DATA] = JACOBIAN_STND_FD(MESH) MESH is NIRFAST mesh structure where
%   'MESH.type' is set to 'stnd'. J is the jacobian for CW (continuous
%   wave) type sources. J is structure of:
%    - J.complex is P (pairs - as enabled in 'MESH.link') by N (nodes)
%       matrix of spatial convolutions of fields (DATA.PHI and DATA.APHI)
%       expressed in photons (photon packets to be exact) per second per
%       millimeter (phot*s^-1*mm^-1). 
%    - J.complete is P (pairs) by N (nodes) matrix of the J_{mua}^{A}
%       sensitivity expressed in millimeters.
%   DATA is the structure as returned by 'femdata_FD' with additional field
%   'aphi'. 'DATA.phi' is the forward field 'DATA.aphi' represents the
%   adjoint field (detectors are new sources) and the boundary data are
%   present in 'DATA.complex'.
%    - DATA.PHI - (direct field) N by M matrix of photon fluence rate at N
%       nodes for M sources as used for source-detector pairs enabled in
%       'MESH.link'. Expressed in mm^-2s^-1 (3D) and mm^-1s^-1 (2D).
%    - DATA.APHI - (adjoint field) N by K matrix of photon fluence rate at
%       N nodes for K detectors used as sources. Detectors as used for
%       source-detector pairs enabled in 'MESH.link' field. Expressed in
%       mm^-2s^-1 (3D) and mm^-1s^-1 (2D).
%    - DATA.COMPLEX - (boundary data) P elements vector of photon fluence
%       rate for P source-detectors pairs as enabled in 'MESH.link'.
%       Expressed in mm^-2s^-1 (3D) and mm^-1s^-1 (2D).
% 
% [J,DATA] = JACOBIAN_STND_FD(MESH,FREQUENCY) Optional FREQUENCY is modulation
%   frequency of sources in Hz.
%    - J.complex is P (pairs) by 2*N (2*nodes) matrix composed of two P by
%       N matrices [J_mus J_mua]. The J_mus is spatial convolution of
%       gradients of fields (DATA.PHI and DATA.APHI) expressed in photon
%       packets per second per millimeter cubed (phot*s^-1*mm^-3). J_mua is
%       the spatial convolution of fields (DATA.PHI and DATA.APHI) as in the
%       real numbers case but carried out using complex numbers arithmetics.
%    - J.complete is 2*P (pairs) by 2*N (2*nodes) matrix of: J_{mua}^{A}
%       (in mm, P by N), J_{kappa}^{A} (in mm^-1, P by N), J_{mua}^{phi}
%       (in mm, P by N) and J_{kappa}^{phi} (in mm^-1, P by N). J.complete
%       is organized as follows:
%       J.complete(1:2:end,:) = [J_{kappa}^{A} J_{mua}^{A}];
%       J.complete(2:2:end,:) = [J_{kappa}^{phi} J_{mua}^{phi}];
% 
% [J,DATA] = JACOBIAN_STND_FD(MESH,FREQUENCY,SECOND_MESH_BASIS)
%   Optional SECOND_MESH_BASIS is a mesh in NIRFAST format that will be
%   used to calculate the Jacobians. This mesh is usually a coarse mesh as
%   compared with the MESH. See e.g. 'pixel_basis' or 'second_mesh_basis'
%   for details how to create such mesh.
% 
% [J,DATA] = JACOBIAN_STND_FD(MESH,FREQUENCY,SECOND_MESH_BASIS,SOLVER) Optional,
%   user specified SOLVER. Type 'help get_solver' for how to use solvers.
% 
% [J,DATA] = JACOBIAN_STND_FD(MESH,FREQUENCY,SECOND_MESH_BASIS,SOLVER,OPTIONS)
%   Optional OPTIONS is a structure that allows to control solvers
%   parameters. Default OPTIONS structure is returned by 'solver_options'
%   function. See 'help solver_options' on how to use this structure.
% 
% [J,DATA] = JACOBIAN_STND_FD(MESH,FREQUENCY,SECOND_MESH_BASIS,SOLVER,OPTIONS,'all') 
%   The 'all' option has effect only if data are real (CW - continuous
%   wave). If specified for the CW data, the attenuation sensitivity to
%   scattering (the 'kappa') is returned as well.
%    - J.complex is P (pairs) by 2*N (2*nodes) matrix composed of two P by
%      N matrices [J_mus J_mua]. The J_mus and J_mua are the same as in the
%      complex numbers case. However, calculations are carried out using
%      real numbers arithmetic.
%    - J.complete is P (pairs) by 2*N (2*nodes) matrix of: J_{kappa}^{A}
%      (in mm, P by N) and J_{mua}^{A} (in mm^-1, P by N). J.complete is
%      organized as follows: J.complete = [J_{kappa}^{A} J_{mua}^{A}]; 
% 
% [J,DATA] = JACOBIAN_STND_FD(MESH,FREQUENCY,[],[],[],'all') Any of the
%   parameters: SECOND_MESH_BASIS, SOLVER and OPTIONS can be set to empty.
%   This will set a given parameter to its default value or no 'coarse'
%   mesh is used.
% 
% [J] = JACOBIAN_STND_FD(___) Does not retur the forward data. Jacobians only.
% 
% See also JACOBIAN_FD, JACOBIAN_SPEC_FD, BUILD_JACOBIAN_FD.
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

data = femdata_FD(mesh, frequency, solver, OPTIONS);

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
data_detector = femdata_FD(mesh, frequency, solver, OPTIONS);
data.aphi = data_detector.phi;

% swap back sources, detectors and link
mesh.source = sources_copy;
mesh.meas = detectors_copy;
mesh.link = link_copy;

%% spatially convolve fileds for sources and detectors

% copy data with possible NaN values for disabled pairs
boundary_data_copy = data.complex;

% jacobian build requires data for enabled pairs only
data.complex = data.complex(logical(data.link(:,3)),:);

% Calculate Jacobian
if ~isempty(second_mesh_basis) % use second mesh basis for Jacobian
    % build jacobian on data interpoladed onto the coarse mesh
    if isAll
        [J] = build_jacobian_FD(second_mesh_basis,interpolatef2r(mesh,second_mesh_basis,data),'all');
    else
        [J] = build_jacobian_FD(second_mesh_basis,interpolatef2r(mesh,second_mesh_basis,data));
    end
else
    % Calculate Jacobian
%     tic;
    if isAll
        [J] = build_jacobian_FD(mesh,data,'all');
    else
        [J] = build_jacobian_FD(mesh,data);
    end
%     data2.aphi = conj(data2.aphi);
%     [J] = build_jacobian_safe_copy(mesh,data2);
%     toc;
end

% put back the boundary data with NaN for disabled pairs
data.complex = boundary_data_copy;

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

function [J,data] = jacobian_stnd_TR(mesh,varargin)
% JACOBIAN_STND_TR Calculates spatial distributions of sensitivity of
%    registerd on the mesh boundary to changes of optical properties
%   per mesh node and moment of the TPSF up to the 2nd order. 
%   This function is specialized to be called for NIRFAST
%   meshes of 'stnd' type (single wavelength).It deliveres up to 4
%   distributions:
%   If data are real (Moments - continuous wave):
%    - J0_{mua}^{A}     - expressesd in mm, 0th moment sensitivity of boundary
%                      attenuation (A) to absorption (mua [mm^-1]). This is
%                      the mean partial pathlength photons travel through
%                      the volume unit. The sum of this Jacobian gives mean
%                      photons pathlength in mm. 
%    - J0_{kappa}^{A}   - (OPTIONAL) expressed in mm^-1, 0th moment sensitivity of
%                      boundary attenuation (A) to diffusion coefficient
%                      'kappa' ((3*(mua+mus))^-1 [mm]). Returned if the
%                      optional parameter is set to 'all'.
%    - J1_{mua}^{A}    1st moment...
%    - J1_{kappa}^{A}  
%    - J2_{mua}^{A}    2nd moment
%    - J2_{kappa}^{A}  
%  
% SYNTAX:
%  [J,DATA] = JACOBIAN_STND_TR(MESH)
%  [J,DATA] = JACOBIAN_STND_TR(MESH,ORDER)
%  [J,DATA] = JACOBIAN_STND_TR(MESH,ORDER,'central')
%  [J,DATA] = JACOBIAN_STND_TR(MESH,ORDER,'central',SECOND_MESH_BASIS)
%  [J,DATA] = JACOBIAN_STND_TR(MESH,ORDER,'central',SECOND_MESH_BASIS,SOLVER)
%  [J,DATA] = JACOBIAN_STND_TR(MESH,ORDER,'central',SECOND_MESH_BASIS,SOLVER,OPTIONS)
%  [J,DATA] = JACOBIAN_STND_TR(MESH,ORDER,'central',SECOND_MESH_BASIS,SOLVER,OPTIONS,'all') 
%  [J,DATA] = JACOBIAN_STND_TR(MESH,ORDER,'central',[],[],[],[])
%
%REFS:
%  [1] Simon R. Arridge and M. Schweiger, "Photon-measurement density functions.
%      Part 2: Finite-element-method calculations," Appl. Opt. 34, 8026-8037 (1995) 
%  [2] Mozumder, M., & Tarvainen, T. (2020). Evaluation of temporal moments
%      and Fourier transformed data in time-domain diffuse optical tomography. JOSA A, 37(12), 1845-1856.
%  [3] Arridge, S. R. (1995). Photon-measurement density functions. 
%      Part I: Analytical forms. Applied Optics, 34(31), 7395-7409.
%  [4] Wabnitz et al., 2020, "Depth-selective data analysis for time-domain fNIRS: moments vs. time windows"


%% check in/out

narginchk(1,7);
nargoutchk(0,2);

%% get optional inputs, handle variable input

% default
order = 2;
central = false;
second_mesh_basis = [];
solver = get_solver;
OPTIONS = solver_options;
isAll = false;

if ~isempty(varargin)
    if length(varargin) >= 1
        % frequency
        if ~ischar(varargin{1}) && ~isstring(varargin{1}) && (numel(varargin{1})==1)
            % sanity check
            if varargin{1} > 2
                error('Maximum supported moment order is 2. Please see help for details on how to use this function.')
            else
                order = 2; 
            end
        else
            error('Bad 2nd argument value. A scalar expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) >= 2
        % do we want centralised moments
        if ischar(varargin{2}) || isstring(varargin{2})
            if strcmp(varargin{2},'central')
                central = true;
            end
        end
    end
    if length(varargin) >= 3
        % second mesh basis
        if isstruct(varargin{3}) || isempty(varargin{3})
            second_mesh_basis = varargin{3};
        else
            error('Bad 3nd argument value. A mesh structure expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) >= 4
        % solver
        if ischar(varargin{4}) || isstring(varargin{4})
            % user specified, sanity check
            solver = get_solver(varargin{4});
        elseif isstruct(varargin{4})
            OPTIONS = varargin{4};
        elseif isempty(varargin{4})
            solver = get_solver;
        else
            error('Bad 4th argument value. Solver name or solver settings structure expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) >= 5
        % solver options
        if isstruct(varargin{5})
            OPTIONS = varargin{5};
        elseif isempty(varargin{5})
            OPTIONS = solver_options;
        else
            error('Bad 5th argument value. Solver settings structure expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) == 6
        % if scattering for CW data as well
        if ischar(varargin{6}) || isstring(varargin{6})
            if strcmp(varargin{6},'all')
                isAll = true;
            end
        else
            error('Bad 6th argument value. Text expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) > 6
        error('Bad arguments. Please see the help for details on how to use this function.')
    end
end

%% If not a workspace variable, load mesh
if ~isstruct(mesh)
    mesh = load_mesh(mesh);
end

%% Calculate field and boundary data for all sources

% data = femdata_FD(mesh, frequency, solver, OPTIONS); From jacobian_stnd_FD
data = femdata_TR_moments(mesh, order, 'field', 'mellin', solver, OPTIONS); % data.phi is 3D matrix (n by m by MAX_ORDER+1).

for i = 1:order
    complex = get_boundary_data(mesh, data.phi(:,:,i));
    if i == 1
        data.complex = zeros([size(complex), order+1]);
    end
    data.complex(:,:,i) = complex;
end

%% Now calculate adjoint source vector (detectors are new sources) for all moments

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
data_detector = femdata_TR_moments(mesh, order, 'field', 'mellin', solver, OPTIONS);
data.aphi = data_detector.phi;
disp(['********** order is ', num2str(order)])
disp('********** size is ')
disp(size(data.aphi))

% swap back sources, detectors and link
mesh.source = sources_copy;
mesh.meas = detectors_copy;
mesh.link = link_copy;

%% spatially convolve fields for sources and detectors

% copy data with possible NaN values for disabled pairs
boundary_data_copy = data.complex;
% data.complex in CW contains (boundary data) P elements vector of photon fluence
%       rate for P source-detectors pairs as enabled in 'MESH.link'.
%       Expressed in mm^-2s^-1 (3D) and mm^-1s^-1 (2D)
% [TODO] this is probably a missing field in TD-NIRS when computing moments
% but we could make it equivalent to the CW NIRS case.
% WARNING:
%   Calculated statistical moments are NOT centralized moments. A q-th
%   moment M_{q} is calculated as follows: 
%     M_{0} = int(Phi(t),dt,0,inf),
%     M_{q>0} = int(t^q * Phi(t),dt,0,inf)/M_{0},
%   where t is time.
%
%   If the 'mellin' option is used, the calculated Mellin transfroms are
%   calculated as follows:
%     M_{0} = int(Phi(t),dt,0,inf),
%     M_{q>0} = int(t^q * Phi(t),dt,0,inf),
%   where t is time.

% jacobian build requires data for enabled pairs only
data.complex = data.complex(logical(data.link(:,3)),:,:); % required for scaling only in FD_NIRS

% This is potentially wrong since it is done inside the call to
% build_jacobian_TR_element, but not going to change anything since it is
% probably all ones
%[TODO] <commented> data.phi = data.phi(logical(data.link(:,3)),:,:);
%[TODO] <commented> data.aphi = data.aphi(logical(data.link(:,3)),:,:);

internal_phi_copy = data.phi;
internal_aphi_copy = data.aphi;

if order >= 0 % See Eq. 26 in REF[1]
    
    disp('YES! Computing 0th moment Jacobian...')
    
    data.phi = internal_phi_copy(:,:,1);
    data.aphi = internal_aphi_copy(:,:,1);

    % Calculate Jacobian
    if ~isempty(second_mesh_basis) % use second mesh basis for Jacobian
        % build jacobian on data interpoladed onto the coarse mesh
        if isAll
            [J0] = build_jacobian_TR_element(second_mesh_basis,interpolatef2r(mesh,second_mesh_basis,data),'all');
        else
            [J0] = build_jacobian_TR_element(second_mesh_basis,interpolatef2r(mesh,second_mesh_basis,data));
        end
    else
        % Calculate Jacobian
        if isAll
            [J0] = build_jacobian_TR_element(mesh,data,'all');
        else
            [J0] = build_jacobian_TR_element(mesh,data);
        end
    end
    
end
    
if order >= 1   % Eq. 27 in REF[1]
    
    disp('YES! Computing 1st moment Jacobian...')
    
    % Calculate Jacobian
    if ~isempty(second_mesh_basis) % use second mesh basis for Jacobian
        % build jacobian on data interpoladed onto the coarse mesh
        if isAll
            data.phi = internal_phi_copy(:,:,1);
            data.aphi = internal_aphi_copy(:,:,2);
            [J1a] = build_jacobian_TR_element(second_mesh_basis,interpolatef2r(mesh,second_mesh_basis,data),'all');
            
            data.phi = internal_phi_copy(:,:,2);
            data.aphi = internal_aphi_copy(:,:,1);
            [J1b] = build_jacobian_TR_element(second_mesh_basis,interpolatef2r(mesh,second_mesh_basis,data),'all');
        else
            data.phi = internal_phi_copy(:,:,1);
            data.aphi = internal_aphi_copy(:,:,2);
            [J1a] = build_jacobian_TR_element(second_mesh_basis,interpolatef2r(mesh,second_mesh_basis,data));
            
            data.phi = internal_phi_copy(:,:,2);
            data.aphi = internal_aphi_copy(:,:,1);
            [J1b] = build_jacobian_TR_element(second_mesh_basis,interpolatef2r(mesh,second_mesh_basis,data));
        end
    else
        % Calculate Jacobian
        if isAll
            data.phi = internal_phi_copy(:,:,1);
            data.aphi = internal_aphi_copy(:,:,2);
            [J1a] = build_jacobian_TR_element(mesh,data,'all');
            
            data.phi = internal_phi_copy(:,:,2);
            data.aphi = internal_aphi_copy(:,:,1);
            [J1b] = build_jacobian_TR_element(mesh,data,'all');
        else
            data.phi = internal_phi_copy(:,:,1);
            data.aphi = internal_aphi_copy(:,:,2);
            [J1a] = build_jacobian_TR_element(mesh,data);
            
            data.phi = internal_phi_copy(:,:,2);
            data.aphi = internal_aphi_copy(:,:,1);
            [J1b] = build_jacobian_TR_element(mesh,data);
        end
    end
    J1 = J1a;
    J1.complex = J1.complex + J1b.complex;
    
end

if order >= 2    % Eq. 27 in REF[1]
    
    disp('YES! Computing 2nd moment Jacobian...')
    
    % Calculate Jacobian
    if ~isempty(second_mesh_basis) % use second mesh basis for Jacobian
        % build jacobian on data interpoladed onto the coarse mesh
        if isAll
            data.phi = internal_phi_copy(:,:,1);
            data.aphi = internal_aphi_copy(:,:,3);
            [J2a] = build_jacobian_TR_element(second_mesh_basis,interpolatef2r(mesh,second_mesh_basis,data),'all');
            
            data.phi = internal_phi_copy(:,:,2);
            data.aphi = internal_aphi_copy(:,:,2);
            [J2b] = build_jacobian_TR_element(second_mesh_basis,interpolatef2r(mesh,second_mesh_basis,data),'all');
            
            data.phi = internal_phi_copy(:,:,3);
            data.aphi = internal_aphi_copy(:,:,1);
            [J2c] = build_jacobian_TR_element(second_mesh_basis,interpolatef2r(mesh,second_mesh_basis,data),'all');
        else
            data.phi = internal_phi_copy(:,:,1);
            data.aphi = internal_aphi_copy(:,:,3);
            [J2a] = build_jacobian_TR_element(second_mesh_basis,interpolatef2r(mesh,second_mesh_basis,data));
            
            data.phi = internal_phi_copy(:,:,2);
            data.aphi = internal_aphi_copy(:,:,2);
            [J2b] = build_jacobian_TR_element(second_mesh_basis,interpolatef2r(mesh,second_mesh_basis,data));
            
            data.phi = internal_phi_copy(:,:,3);
            data.aphi = internal_aphi_copy(:,:,1);
            [J2c] = build_jacobian_TR_element(second_mesh_basis,interpolatef2r(mesh,second_mesh_basis,data));
        end
    else
        % Calculate Jacobian
        if isAll
            data.phi = internal_phi_copy(:,:,1);
            data.aphi = internal_aphi_copy(:,:,3);
            [J2a] = build_jacobian_TR_element(mesh,data,'all');
            
            data.phi = internal_phi_copy(:,:,2);
            data.aphi = internal_aphi_copy(:,:,2);
            [J2b] = build_jacobian_TR_element(mesh,data,'all');
            
            data.phi = internal_phi_copy(:,:,3);
            data.aphi = internal_aphi_copy(:,:,1);
            [J2c] = build_jacobian_TR_element(mesh,data,'all');
        else
            data.phi = internal_phi_copy(:,:,1);
            data.aphi = internal_aphi_copy(:,:,3);
            [J2a] = build_jacobian_TR_element(mesh,data);
            
            data.phi = internal_phi_copy(:,:,2);
            data.aphi = internal_aphi_copy(:,:,2);
            [J2b] = build_jacobian_TR_element(mesh,data);
            
            data.phi = internal_phi_copy(:,:,3);
            data.aphi = internal_aphi_copy(:,:,1);
            [J2c] = build_jacobian_TR_element(mesh,data);
        end
    end
    J2 = J2a;
    J2.complex = J2.complex + 2 * J2b.complex + J2c.complex;
    
end

data.link = mesh.link;

J = [];

% Normalise to time integrated boundary data. 
% This is equivalent to: J_{\alpha, \nu} / E  
if order >= 0
    % Scale to the Rytow approximation (natural log of boundary data change)
    % For this moment this is what nirfaster usually does.
    % If we apply the other transforms it would be zero.
    J0.complete = J0.complex ./ data.complex(:,:,1);
    J.zeroth = J0;
    disp('0th moment Jacobian normalised')
end

% For higher orders, to compute the jacobians of the temporal moments 
% as oposed to computing the so far computed jaconians of the Mellin
% moments we need to apply Eq. 41 of REF[3], i.e.
% J^{<t_n>} = J^{T_n} / E - <t_n> * J^E / E ...
%           = J^{T_n} / E - T_n * J^E / E^2, as per Eq. 9 in REF[3].
% TODO - would it be worth computing J^E with the dedicated jacobian
% calculator as well as E with the corresponding function, not assumming it
% is the same than the moment of order 0?

if order >= 1
    J1.complete = J1.complex / data.complex(:,:,1)...
                  - data.complex(:,:,2) .* J0.complex / data.complex(:,:,1).^2;
    J.first = J1;
    disp('1st moment Jacobian normalised')
end

if order >= 2
    J2.complete = J2.complex / data.complex(:,:,1)...
                  - data.complex(:,:,3) .* J0.complex / data.complex(:,:,1).^2;
    J.second = J2;
    disp('2nd moment Jacobian normalised')
    if central  % REF[4] - derivative of equation (3)
       disp('Centralising second moment')
       J.second.complete = J.second.complete - 4*J.first.complete;
    end
end

% put back the boundary data with NaN for disabled pairs and for all moments
data.phi = internal_phi_copy;
data.aphi = internal_aphi_copy;

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

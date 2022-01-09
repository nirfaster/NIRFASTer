function [data,varargout] = femdata_DCS(mesh,tau_DCS,varargin)
% FEMDATA_DCS Calculates diffuse correlation spectroscopy (DCS) data
%   (electric field temporal autocorrelation function) for a given mesh at
%   a given correlation delay times.  
% 
% SYNTAX:
%  [DATA] = FEMDATA_DCS(MESH, TAU_DCS)
%  [DATA] = FEMDATA_DCS(MESH, TAU_DCS, WV_DCS)
%  [DATA] = FEMDATA_DCS(MESH, TAU_DCS, WV_DCS, SOLVER)
%  [DATA] = FEMDATA_DCS(MESH, TAU_DCS, WV_DCS, OPTIONS)
%  [DATA] = FEMDATA_DCS(MESH, TAU_DCS, WV_DCS, SOLVER, OPTIONS)
%  [DATA] = FEMDATA_DCS(MESH, TAU_DCS, [], [], [])
%  [DATA,INFO] = FEMDATA_DCS(___)
%
% [DATA] = FEMDATA_DCS(MESH, TAU_DCS) Calculates diffuse correlation
%   spectroscopy (DCS) data for a given mesh MESH in NIRFAST format and at
%   a given correlation delay times TAU_DCS. It is required the MESH have
%   following fields:
%    - MESH.aDb - nodal values of the 'flow index' expressed in mm^2*s^-1.
%    - MESH.wv_DCS - a scalar, the DCS wavelength in nm, e.g. 785
%   TAU_DCS expressed in seconds, is a vector of k requested
%    autocorrelation delay time values.
%   DATA structure has following entries:
%    - DATA.phi - n by m matrix of photon fluence rate at n nodes for m
%       sources as enabled in 'MESH.link' field. Expressed in mm^-2s^-1.
%       This is the electric field temporal autocorrelation function when
%       the autocorrelation delay time equals zero.
%    - DATA.amplitude - p elements vector of photon fluence rate for p
%       source-detectors pairs as specified in MESH.link. Expressed in
%       mm^-2s^-1. The electric field temporal autocorrelation function for
%       the autocorrelation delay time equal zero.
%    - DATA.link - copy of MESH.link showing how sources and detectors are
%       linked together into pairs.
%    - DATA.tau_DCS - copy of the TAU_DCS, vector of k requested
%       autocorrelation delay times in seconds.
%    - DATA.phi_DCS - n by m by k matrix of electric field temporal
%       autocorrelation function at n MESH nodes, originating at m sources
%       and sampled at k autocorrelation delay times. Expressed in
%       mm^-2s^-1.
%    - DATA.G1_DCS - p by k matrix of the electric field temporal
%       autocorrelation function at the mesh boundary (at detectors) for p
%       source-detctor pairs as specified in MESH.link and sampled at k
%       autocorrelation delay times. Expressed in mm^-2s^-1.
% 
% [DATA] = FEMDATA_DCS(MESH, TAU_DCS, WV_DCS) WV_DCS is an override of
%   'MESH.wv_DCS'. It specifies the DCS long-coherence light source
%   wavelengths. Expressed in nm.
% 
% [DATA] = FEMDATA_DCS(MESH, TAU_DCS, WV_DCS, SOLVER) User specified SOLVER
%   is used. See help of the 'get_solver' function for how to use solvers.
% 
% [DATA] = FEMDATA_DCS(MESH, TAU_DCS, WV_DCS, OPTIONS) Default solver is
%   used. See help of the 'get_solvers' function for how the default solver
%   is determined. OPTIONS is a structure that allows to control solvers
%   parameters. The OPTIONS structure can have following optional entries:
%       - OPTIONS.no_of_iter (default 1000)
%         Maximum number of BiCGStab iterations.
%       - OPTIONS.tolerance (default 1e-8);
%         Absolute convergence tolerance (norm of solution error).
%       - OPTIONS.rel_tolerance (default 1e-8);
%         Relative convergence tolerance (norm of solution error related to
%         norm of error of first iteration).
%       - OPTIONS.divergence_tol (default 1e8);
%         Relative divergence tolerance (related to first iteration).
%       - OPTIONS.GPU (default -1);
%         Index of GPU requested to be used. By default (-1) the GPU of the
%         highest compute capability is used. Use the 'isCUDA' function to
%         get the installed GPUs info. The indexing is 0-based and the
%         indexing order is the same as in the list of GPUs returned by the
%         'isCUDA' function.
%   Entries within the structure are optional. The default OPTIONS
%   structure is returned by 'solver_options' function.
% 
% [DATA] = FEMDATA_DCS(MESH, TAU_DCS, WV_DCS, SOLVER, OPTIONS) User
%   specified SOLVER and its OPTIONS are used. See help of the
%   'get_solvers' function for how to use solvers. The default OPTIONS
%   structure is returned by 'solver_options' function. 
% 
% [DATA] = FEMDATA_DCS(MESH, TAU_DCS, [], [], []) Any of the parameters:
%   WV_DCS, SOLVER and OPTIONS can be set to empty. This will set a given
%   parameter to its default value.
% 
% [DATA,INFO] = FEMDATA_DCS(___) Also returns solvers
%   information/statistic at TAU_DCS=0 in the INFO structure with following
%   possible entries:
%    - INFO.convergence_iterations (m by 1 vector)
%      Number of iterations to converge to the tolerances for m sources as
%      enabled in 'MESH.link'
%    - INFO.convergence_tolerances (m by 1 vector)
%      Exact absolute tolerance at convergence for m sources as enabled in
%      'MESH.link'.
%    - INFO.isConverged (m by 1 vector)
%      Logical flag if solution converged to tolerances for m sources as
%      enabled in 'MESH.link'.
%    - INFO.isConvergedToAbsTol (m by 1 vector)
%      Logical flag if solution converged to the absolute tolerance for
%      m sources as enabled in 'MESH.link'.
%    - INFO.convergence_iterations_limit
%      Maximum number of convergence iterations allowed.
%    - INFO.convergence_tolerance_limit
%      Absolute convergence tolerance requested.
%    - INFO.convergence_relative_tolerance_limit
%      Relative convergence tolerance (related to first iteration)
%      requested. 
%    - INFO.divergence_tolerance_limit
%      Maximum relative divergence tolerance (related to first
%      iteration) allowed. 
%    - INFO.status (m by 1 vector)
%      Convergence flag (from 0 to 4) as returned by MATLAB 'bicgstab'. For
%      all m sources as enabled in 'MESH.link'. '0' means success.
%   INFO is empty if the '\' solver is used.
% 
% NOTE 1:
%   The 'flow index' expressed in mm^2*s^-1 represents an effective
%   diffusion coefficient of the moving scatterers (Db) weighted by the 
%   concentration of moving particles (a). The (Db) can be identified with
%   mean-square displacement of moving particels (Brownian motion
%   coefficient) and (a) is a percentage of light scattering events from
%   moving scatterers (as related to all scattering events).
% 
% NOTE 2:
%   The normalized electric field temporal autocorrelation function can be
%   calculated by normalizing with values at the correlation time equal
%   zero. I.e, the normalized autocorrelation fields within the MESH
%   (g1_3D) for a source and  curves (g1_curve) for source-detector pairs
%   can be calculated as follows:
%     g1_3D(:,ind_source,:) = DATA.phi_DCS(:,ind_source,:)./DATA.phi(:,ind_source); 
%     g1_curve = DATA.G2_DCS./DATA.amplitude;
%   The g1 values are unitless and vary between 1 (full correlation) and 0
%   (no correlation).
% 
% WARNING:
%   Please meke sure the MESH optical properties (MESH.mua, MESH.mus,
%   MESH.kappa) are given for the wv_DCS wavelength.
% 
% READ: Durduran T. et al. Diffuse optics for tissue monitoring and
%       tomography. Rep. Prog. Phys, 73, 2010, 076701(43pp)
% 
% See also FEMDATA_STND_FD, GET_SOLVER, SOLVER_OPTIONS
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

%% check in/out

narginchk(2,5);
nargoutchk(0,2);

if isempty(tau_DCS)
    error('Empty autocorrelation time. Please see help for details on how to use this function');
end

%% get optional inputs, handle variable input

% default
wv_DCS = [];
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
        if ischar(varargin{2}) || isstring(varargin{2})
            % user specified, sanity check
            solver = get_solver(varargin{2});
        elseif isstruct(varargin{2})
            OPTIONS = varargin{2};
        elseif ~isempty(varargin{2})
            warning('Bad 4th argument value. Text, structure or empty value expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) == 3        
        if isstruct(varargin{3})
            OPTIONS = varargin{3};
        elseif ~isempty(varargin{3})
            warning('Bad 5th argument value. Structure or empty value expected. Please see help for details on how to use this function.')
        end
    end
    if length(varargin) > 3
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

%% TAU=0. Firs, get the field for tau=0 (zero autocorrelation time)

warning('off','NIRFAST:warning:meshType');
% use the continuous wave source
if nargout == 2
    [data_diffusion, varargout{1}] = femdata_stnd_FD(mesh,0,solver,OPTIONS);
else
    data_diffusion = femdata_stnd_FD(mesh,0,solver,OPTIONS);
end
warning('on','NIRFAST:warning:meshType');

% assign the tau=0 results
data.phi = data_diffusion.phi;
data.amplitude = data_diffusion.amplitude;
data.link = data_diffusion.link;

%% declare the DCS data part

% number of taus
no_of_tau = numel(tau_DCS);
% copy of taus
data.tau_DCS = tau_DCS;
% autocorrelation fields
data.phi_DCS = zeros(size(data_diffusion.phi,1),size(data_diffusion.phi,2),no_of_tau);
% autocorrelation boundary data (not normalized!)
data.G1_DCS = zeros(size(data_diffusion.amplitude,1),no_of_tau);

%% TAU>0. Calculate DCS data at tau>0

% copy mua as it will be modified for DCS
mua_tau_0 = mesh.mua;
% wave number
k0 = 2*pi*mesh.ri./(wv_DCS * 1e-6); % <- DCS wavelength converted from nm to mm

warning('off','NIRFAST:warning:meshType');
% loop through autocorrelation times
for ind_tau = 1:no_of_tau    
    % new 'autocorrelation' absorption
    mesh.mua = mua_tau_0 + 2 * mesh.aDb .* k0.^2 .* (mesh.mus + mua_tau_0) * tau_DCS(ind_tau);

	% use the continuous wave source
    data_tau = femdata_stnd_FD(mesh,0,solver,OPTIONS);
    
    % fill in the results
    data.phi_DCS(:,:,ind_tau) = data_tau.phi;
    data.G1_DCS(:,ind_tau) = data_tau.amplitude(:,1);
end

warning('on','NIRFAST:warning:meshType');

end


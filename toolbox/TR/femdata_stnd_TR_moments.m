function [data] = femdata_stnd_TR_moments(mesh, varargin)
% FEMDATA_STND_TR_MOMENTS Returns spatial distributions and boundary values
%   of statistical moments or Mellin transforms of time-resolved data for
%   the 'standard' FEM mesh in NIRFAST format. This function generates
%   sparse matrices representing the FEM problem (the time-resolved moments
%   'mass' mastrices and sparse sources vectors). Then, it calls the
%   time-resolved moments solver. 
%
% SYNTAX: 
%   [DATA] = FEMDATA_STND_TR_MOMENTS(MESH)
%   [DATA] = FEMDATA_STND_TR_MOMENTS(MESH, MAX_ORDER)
%   [DATA] = FEMDATA_STND_TR_MOMENTS(MESH, MAX_ORDER, 'field')
%   [DATA] = FEMDATA_STND_TR_MOMENTS(MESH, MAX_ORDER, 'field', 'mellin')
%   [DATA] = FEMDATA_STND_TR_MOMENTS(MESH, MAX_ORDER, 'field', 'mellin', SOLVER)
%   [DATA] = FEMDATA_STND_TR_MOMENTS(MESH, MAX_ORDER, 'field', 'mellin', OPTIONS)
%   [DATA] = FEMDATA_STND_TR_MOMENTS(MESH, MAX_ORDER, 'field', 'mellin', SOLVER, OPTIONS)
%   [DATA] = FEMDATA_STND_TR_MOMENTS(MESH, [], [], [], [], [])
%
% [DATA] = FEMDATA_STND_TR_MOMENTS(MESH)
%   Returns 0th, 1st and 2nd statistical moments of time-resolved temporal
%   point spread functions (distributions of time of flight of photons) for
%   sources and detectors as specified in the MESH structure (MESH.link).
%   MESH - FEM mesh in NIRFAST format. Standard mesh (single wavelength).
%   The output stucture DATA has following entries:
%     - DATA.phi [optional] is 3D matrix (n by m by k) of the spatial
%       distributions of the statistical moments where n - number of mesh
%       nodes, m - number of sources, k - number of statistical moments
%       requested. This value is optional and is enabled by the optional
%       input parameter 'field'. The moment unit depends ont its order q,
%       starting at the 0th order (q=0,1,2,...): phot*s^q.
%     - DATA.moments has size of (p by k) where p represents
%       source-detector combinations as specified in the 'MESH.link' and k
%       holds the requested moments. Moments units are the same as in the
%       'DATA.phi'.
%     - DATA.link is copy of the MESH.link showing available
%       source-detector combinations.
%   By default this function calculates 0th, 1st and 2nd statistical
%   moments, does not return the 'DATA.phi', uses default solver (see
%   'get_solver') and options (see 'solver_options').
% 
% WARNING: If the MESH.link specifies that a given source-detecor pair is
%          OFF, the DATA.moments values at that pair will be NaN. 
% 
% [DATA] = FEMDATA_STND_TR_MOMENTS(MESH, MAX_ORDER)
%   Returns 0th, 1st, 2nd up to MAX_ORDER statistical moments of
%   time-resolved temporal point spread functions (distributions of time of
%   flight of photons) for sources and detectors as specified  in the MESH
%   structure (MESH.link). 
%   MAX_ORDER should be >=0. As the moments order starts at 0, the returned
%   DATA fields size is related to the MAX_ORDER as follows
%     - DATA.phi [optional] is 3D matrix (n by m by MAX_ORDER+1).
%     - DATA.moments is (p by MAX_ORDER+1).
% 
% [DATA] = FEMDATA_STND_TR_MOMENTS(MESH, MAX_ORDER, 'field')
%   Returns the optional 'DATA.phi' field with spatial distributions fo the
%   statistical moments.
% 
% [DATA] = FEMDATA_STND_TR_MOMENTS(MESH, MAX_ORDER, 'field', 'mellin')
%   'mellin' stands for the Mellin transform of the time-dependent
%   intensity. If the 'mellin' option is specified, the returned 'DATA'
%   structure holds the Mellin transforms up to the degree of MAX_ORDER+1
%   (0th, 1st, 2nd up to MAX_ORDER).
% 
% [DATA] = FEMDATA_STND_TR_MOMENTS(MESH, MAX_ORDER, 'field', 'mellin', SOLVER) User
%   specified SOLVER is used. See help of the 'get_solvers' function for
%   how to use solvers.
% 
% [DATA] = FEMDATA_STND_TR_MOMENTS(MESH, MAX_ORDER, 'field', 'mellin', OPTIONS)
%   OPTIONS is a structure that allows to control solvers parameters. The
%   OPTIONS structure can have following optional entries:
%      - OPTIONS.no_of_iter (default 1000)
%        Maximum number of BiCGStab iterations.
%      - OPTIONS.tolerance (default 1e-12 and squares with moments order);
%        Absolute convergence tolerance (norm of solution error).
%      - OPTIONS.rel_tolerance (default 1e-8);
%        Relative convergence tolerance (norm of solution error related to
%        norm of error of first iteration).
%      - OPTIONS.divergence_tol (default 1e8);
%        Relative divergence tolerance (related to first iteration).
%   Entries within the structure are optional. The default OPTIONS
%   structure is returned by 'solver_options' function.
%   The absolute tolerance 'OPTIONS.tolerance' as returned by the
%   'solver_options' function is 1e-8. Here we anhance it to 1e-12 and it
%   is squared (inside mex C++ functions) as the absolute values of moments
%   decrease rapidly (roughly in a square manner).
% 
% [DATA] = FEMDATA_STND_TR_MOMENTS(MESH, MAX_ORDER, 'field', 'mellin', SOLVER, OPTIONS) 
%   Gives control on both SOLVER and OPTIONS.
% 
% [DATA] = FEMDATA_STND_TR_MOMENTS(MESH, [], [], [], [], []) Optional
%   fields can be empty, e.g. it is valid to call this function as 
%   [DATA] = FEMDATA_STND_TR_MOMENTS(MESH, [], 'field', [], [], OPTIONS)
% 
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
% 
% See also GEN_MASS_MATRIX_TR_MOMENTS_CUDA, GEN_MASS_MATRIX_TR_MOMENTS_CPU,
%          GET_FIELD_TR_MOMENTS_CUDA, FEMDATA_SPEC_TR_MOMENTS, GET_SOLVER,
%          SOLVER_OPTIONS, GEN_SOURCES. 
% 
% Read: Arridge S.R. and M. Schweiger Photon-measurement density
%       functions. Part 2: Finite-element-method calculations. Applied
%       Optics, Vol. 34, No. 34, 1995, p. 8026-8037
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

%% check in/out

narginchk(1,6);
nargoutchk(0,1);

%% serve the variable input 

% is the field at mesh nodes should be returned 
isField = false;
isMellin = false;
max_order = 2;
solver = get_solver;
OPTIONS = solver_options;

if ~isempty(varargin)
    if length(varargin) >= 1
        if ~ischar(varargin{1}) && ~isstring(varargin{1}) && ~isempty(varargin{1})
            max_order = round(varargin{1});
            if max_order < 0
                error('Bad maximum order. Should be >=0. Please see the help for details on how to use this function.')
            end
        end
    end
    if length(varargin) >= 2
        if ischar(varargin{2}) || isstring(varargin{2})
            if strcmp(varargin{2}, 'field')
                isField = true;
            else
                warning(['Unknown option ''' varargin{2} '''. Spatial distributions of time-resolved moments will not be returned. Please see the help.'])
            end
        end
    end
    if length(varargin) >= 3
        if ischar(varargin{3}) || isstring(varargin{3})
            if strcmp(varargin{3}, 'mellin')
                isMellin = true;
            else
                warning(['Unknown option ''' varargin{3} '''. Spatial distributions of time-resolved Mellin transforms will not be returned. Please see the help.'])
            end
        end
    end
    if length(varargin) >= 4
        if ischar(varargin{4}) || isstring(varargin{4})
            % user specified, sanity check
            solver = get_solver(varargin{4});
        elseif isstruct(varargin{4})
            OPTIONS = varargin{4};
        end
    end
    if length(varargin) == 5
        if isstruct(varargin{5})
            OPTIONS = varargin{5};
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

if ~strcmp(mesh.type,'stnd')
    warning(['Mesh type is ''' mesh.type '''. ''stnd'' expected. This might give unexpected results.'])
end

%% Now calculate sources vector for those active in mesh.link

% real values of sources in time domain
qvec = real(gen_sources(mesh));


%% Calculate forward data
% per mesh node for all sources at all time bins

% number of source-detector pairs
num_meas = size(mesh.link,1);


% make FEM matrices
if isCUDA
    [i_index, j_index ,value_1, value_2] = gen_mass_matrix_TR_moments_CUDA(mesh);
else
    [i_index, j_index ,value_1, value_2] = gen_mass_matrix_TR_moments_CPU(mesh);
end

% check for unimplemented solvers
if strcmp(solver,solver_name_CPU)
    % get spatial distribution of moments
    warning(['Parallel CPU (''get_field_TR_moments_CPU'') not supported yet.... ' ...
        'Using MATLAB backslash instead.'])
elseif strcmp(solver,solver_name_matlab_iterative)
    % get photon fluence rate
    warning(['MATLAB BiCGStab (''get_field_TR_moments_bicgstab_matlab'') not supported yet.... ' ...
        'Using MATLAB backslash instead.'])
end

% if GPU in use
if strcmp(solver,solver_name_GPU)
    % get spatial distribution of moments
    [phi, ~] = get_field_TR_moments_CUDA(i_index, j_index, value_1, value_2, qvec, max_order+1, OPTIONS);
% if no GPU in use, use CPU
else
    % use MATLAB backslash
    phi = zeros(size(qvec,1),size(qvec,2),max_order+1); % nodal data per source per statistical moment
    % +1 as the i_index and j_index are zero-based and MATLAB uses 1-based indexing
    A = sparse(i_index+1,j_index+1,value_1);
    B = sparse(i_index+1,j_index+1,value_2);
    % A * M_{0} =  QVEC, calculate zeroth moment M_{0} - the CW data!
    phi(:,:,1) = full(A\qvec);
    % loop through moments > 0 (mean time, variance, etc.)
    for ind_order = 1:max_order
        % QVEC_M = B * M_{ind_order-1}
        qvec_m = B * squeeze(phi(:,:,ind_order));
        % M_{ind_order} = (ind_order*A)^-1 * B * M_{ind_order-1}
        phi(:,:,ind_order + 1) =  1/ind_order * A\qvec_m;
    end
end

% finish calculations for statistical moments, normalize Mellin moments to statistical moments
if ~isMellin
    for ind_order = 2:size(phi,3)
        phi(:,:,ind_order) = phi(:,:,ind_order) ./ phi(:,:,1);
    end
end

% add raw data if needed
if isField
    data.phi = phi;
end

% get boundary moments from spatial distributions
data.moments = zeros(num_meas,size(phi,3));
for ind_order = 1:size(data.moments,2)
    data.moments(:,ind_order)= get_boundary_data(mesh,phi(:,:,ind_order));
end

% set source-detector link
data.link = mesh.link;

end
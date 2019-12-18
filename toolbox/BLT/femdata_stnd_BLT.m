function [data, varargout] = femdata_stnd_BLT(mesh, varargin)
% FEMDATA_STND_FD Calculates data (phase and amplitude) for a given standard
%   mesh at a given frequency (Hz). Outputs phase and amplitude in the
%   structure DATA.
% 
% SYNTAX:
%  [DATA] = FEMDATA_STND_FD(MESH, FREQUENCY)
%  [DATA] = FEMDATA_STND_FD(MESH, FREQUENCY, SOLVER)
%  [DATA] = FEMDATA_STND_FD(MESH, FREQUENCY, OPTIONS)
%  [DATA] = FEMDATA_STND_FD(MESH, FREQUENCY, SOLVER, OPTIONS)
%  [DATA,INFO] = FEMDATA_STND_FD(___)
% 
% [DATA] = FEMDATA_STND_FD(MESH, FREQUENCY) MESH is NIRFAST mesh
%   structure, FREQUENCY is source frequency (>=0) in Hz. A default solver
%   is used. See help of the 'get_solvers' function for default solver.
%   DATA structure has following entries:
%    - DATA.phi - N by M matrix of photon fluence rate at N nodes for M
%       sources as enabled in 'MESH.link' field. Expressed in mm^-2s^-1.
%    - DATA.complex - P elements vector of photon fluence rate for P
%       source-detectors pairs as specified in MESH.link. This vector has
%       COMPLEX entries. Expressed in mm^-2s^-1. 
%    - DATA.link - copy of MESH.link showing how sources and detectors are
%       linked together into pairs.
%    - DATA.amplitude - absolute value of the DATA.complex. Expressed in
%       mm^-2s^-1.
%    - DATA.phase - angle (polar coordinate) of DATA.complex. Expressed in
%       degrees.
%    - DATA.paa - two colums matrix [DATA.amplitude DATA.phase];
% 
% [DATA] = FEMDATA_STND_FD(MESH, FREQUENCY, SOLVER) User specified SOLVER
%   is used. See help of the 'get_solvers' function for how to use solvers.
% 
% [DATA] = FEMDATA_STND_FD(MESH, FREQUENCY, OPTIONS) Default solver is
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
%   Entries within the structure are optional. The default OPTIONS
%   structure is returned by 'solver_options' function.
% 
% [DATA] = FEMDATA_STND_FD(MESH, FREQUENCY, SOLVER, OPTIONS) User specified
%   SOLVER and its OPTIONS are used. See help of the 'get_solvers' function
%   for how to use solvers. The default OPTIONS structure is returned by
%   'solver_options' function. 
% 
% [DATA,INFO] = FEMDATA_STND_FD(___) Also returns solvers
%   information/statistic in the INFO structure with following possible
%   entries:
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
% WARNING
%   For large source-detector separations (>8cm).
%   Iterative solvers are set to use the follwoing absolute tolerance
%   'OPTIONS.tolerance=1e-8'. The absolute tlerance can be seen as a
%   maximum difference between iterative method result and direct solver
%   (MATLAB backslash '\'). 
%   As a rule of thumb, one order of the abosolute tolerance corresponds to
%   1cm of source-detector distance on a tissue. As such, the 1e-8 should
%   provide accurate results at least for 8cm source-detector distance.
%   Decrease the tolerance acordingly if you want to calculate data at
%   source-detector distance >8cm. However, the decreased tolerance
%   increases the execution time.
% 
% See also GET_SOLVER, SOLVER_OPTIONS, FEMDATA_SPEC_FD, GET_FIELD_FD_CPU,
%          GET_FIELD_FD_CUDA, GEN_MASS_MATRIX_FD_CPU,
%          GEN_MASS_MATRIX_FD_CUDA, GEN_SOURCES.
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz and H. Dehghani 2018

%% check in/out

narginchk(1,3);
nargoutchk(0,2);

%% get optional inputs, handle variable input

% default
solver = get_solver;
OPTIONS = solver_options;

if ~isempty(varargin)
    if length(varargin) == 1
        if ischar(varargin{1}) || isstring(varargin{1})
            % user specified, sanity check
            solver = get_solver(varargin{1});
        elseif isstruct(varargin{1})
            OPTIONS = varargin{1};
        else
            error('Bad 3rd argument value. Text or structure expected. Please see help for details on how to use this function.')
        end
    elseif length(varargin) == 2
        
        if ischar(varargin{1}) || isstring(varargin{1})
            % user specified, sanity check
            solver = get_solver(varargin{1});
        else
            error('Bad 2nd argument value. Text expected. Please see help for details on how to use this function.')
        end
        
        if isstruct(varargin{2})
            OPTIONS = varargin{2};
        else
            error('Bad 3rd argument value. Structure expected. Please see help for details on how to use this function.')
        end
    else
        error('Bad arguments. Please see the help for details on how to use this function.')
    end
end

%% If not a workspace variable, load mesh
if ~isstruct(mesh)
    mesh = load_mesh(mesh);
end

if ~strcmp(mesh.type,'stnd')
    warning('NIRFAST:warning:meshType',['Mesh type is ''' mesh.type '''. ''stnd'' expected. This might give unexpected results.'])
end

%% Create FEM matrix

% make FEM matrices
if isCUDA
    [i_index, j_index, value] = gen_mass_matrix_FD_CUDA(mesh,0);
else
    [i_index, j_index, value] = gen_mass_matrix_FD_CPU(mesh,0);
end

%% Now calculate sources vector for those active in mesh.link
if isfield(mesh,'blt')
    qvec = sparse(abs(mesh.blt));
else
    warning('NIRFAST:warning: BLT source not defined...exiting');
    data = [];
    return
end
    
%% Calculate field for all sources

% OPTIONS.no_of_iter = 1000;% (default 1000)
% OPTIONS.tolerance = 1e-8;% (default 1e-12 for DOUBLE and 1e-8 for SINGLE);
% % below optins ignored for MATLAB iterative method
% OPTIONS.rel_tolerance = 1e-8;% (default 1e-12 for DOUBLE and 1e-8 for SINGLE);
% OPTIONS.divergence_tol = 1e8;% (default 1e8 for DOUBLE and SINGLE);

% get photon fluence rate
% if GPU in use
if strcmp(solver,solver_name_GPU)
    % check if we need to return the info structureas well
    if nargout == 2
        [data.phi, varargout{1}] = get_field_FD_CUDA(i_index, j_index, value, qvec, OPTIONS);
    else
        data.phi = get_field_FD_CUDA(i_index, j_index, value, qvec, OPTIONS);
    end
% if no GPU in use, use CPU
elseif strcmp(solver,solver_name_CPU)
    if nargout == 2
        [data.phi, varargout{1}] = get_field_FD_CPU(i_index, j_index, value, qvec, OPTIONS);
    else
        data.phi = get_field_FD_CPU(i_index, j_index, value, qvec, OPTIONS);
    end
elseif strcmp(solver,solver_name_matlab_iterative)
    if nargout == 2
        [data.phi, varargout{1}] = get_field_FD_bicgstab_matlab(i_index, j_index, value, qvec, OPTIONS);
    else
        data.phi = get_field_FD_bicgstab_matlab(i_index, j_index, value, qvec, OPTIONS);
    end
else
    % use MATLAB backslash
    % +1 as the i_index and j_index are zero-based and MATLAB uses 1-based indexing
    data.phi = full(sparse(i_index+1, j_index+1, value)\qvec);
    if nargout == 2
        varargout{1} = [];
    end
end

%% Trap for link file.
% Measure all detectors for BLT
mesh.link = ones(length(mesh.meas.coord),3);
mesh.link(:,2) = 1 : length(mesh.meas.coord);

% Calculate boundary data
[data.complex]=get_boundary_data(mesh,data.phi);

%% Format output data
data.amplitude = abs(data.complex);
data = rmfield(data,'complex');

end

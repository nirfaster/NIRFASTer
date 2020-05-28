function [data,varargout] = femdata_FD(mesh, frequency, varargin)
% FEMDATA_FD Calculates data (phase and amplitude) for a given mesh at a
%   given frequency (Hz). Outputs phase and amplitude in the structure
%   DATA. This function is a wrapper to call mesh-type-specific femdata
%   functions.
% 
% SYNTAX:
%  [DATA] = FEMDATA_FD(MESH, FREQUENCY)
%  [DATA] = FEMDATA_FD(MESH, FREQUENCY, SOLVER)
%  [DATA] = FEMDATA_FD(MESH, FREQUENCY, OPTIONS)
%  [DATA] = FEMDATA_FD(MESH, FREQUENCY, SOLVER, OPTIONS)
%  [DATA,INFO] = FEMDATA_FD(___)
% 
% [DATA] = FEMDATA_FD(MESH, FREQUENCY) MESH is NIRFAST mesh
%   structure, FREQUENCY is source frequency (>=0) in Hz. A default solver
%   is used. See help of the 'get_solvers' function for default solver.
%   DATA structure has following entries:
%    - DATA.PHI - N by M matrix of photon fluence rate at N nodes for M
%       sources as enabled in 'MESH.link' field. Expressed in mm^-2s^-1.
%    - DATA.COMPLEX - P elements vector of photon fluence rate for P
%       source-detectors pairs as specified in MESH.LINK. This vector has
%       COMPLEX entries. Expressed in mm^-2s^-1. 
%    - DATA.LINK - copy of MESH.LINK showing how sources and detectors are
%       linked together into pairs.
%    - DATA.AMPLITUDE - absolute value of the DATA.COMPLEX. Expressed in
%       mm^-2s^-1.
%    - DATA.PHASE - angle (polar coordinate) of DATA.COMPLEX. Expressed in
%       degrees.
%    - DATA.PAA - two colums matrix [DATA.AMPLITUDE DATA.PHASE];
% 
% [DATA] = FEMDATA_FD(MESH, FREQUENCY, SOLVER) User specified SOLVER
%   is used. See help of the 'get_solvers' function for how to use solvers.
%
% [DATA] = FEMDATA_FD(MESH, FREQUENCY, OPTIONS) Default solver is
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
% [DATA] = FEMDATA_FD(MESH, FREQUENCY, SOLVER, OPTIONS) User specified SOLVER
%   is used. See help of the 'get_solvers' function for how to use solvers.
%   The default OPTIONS structure is returned by 'solver_options' function.
% 
% [DATA,INFO] = FEMDATA_FD(__) Also returns solvers
%   information/statistic in the INFO structure with following possible
%   entries:
%    - INFO.convergence_iterations (m by 1 vector)
%      Number of iterations to converge to the tolerances for m column
%      of QVEC (right hand side).
%    - INFO.convergence_tolerances (m by 1 vector)
%      Exact absolute tolerance at convergence for m column of QVEC
%      (right hand side).
%    - INFO.isConverged (m by 1 vector)
%      Logical flag if solution converged to tolerances for m column of
%      QVEC (right hand side).
%    - INFO.isConvergedToAbsTol (m by 1 vector)
%      Logical flag if solution converged to the absolute tolerance for
%      m column of QVEC (right hand side).
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
%      all m column of QVEC (right hand side). '0' means success.
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
% GUIDELINE
%   Please use the 'MESH.link' field to turn on/off wavelength at sources.
%   This field has following row entries:
%    #sour #det 0/1 0/1 ... 0/1
%   #sour - source label as in 'MESH.source.num', #det - detector label as
%   in 'MESH.meas.num', 0/1 this pair is off (0) or on (1), column position
%   of the 0/1 entries correspond to wavelengths as specified in 'MESH.wv'.
%   As such, if you want to disable a 2nd and 5th wavelengths at sources
%   labeled 13 and 16, please change the 'MESH.link' as follows:
%     sources2disable = [13 16];
%     wavelengths2disable = [2 5];
%     mask_sources = sum(mesh.link(:,1)==sources2disable,2) > 0;
%     mesh.link(mask_sources,wavelengths2disable+2) = 0;
% 
% See also FEMDATA_STND_FD, FEMDATA_SPEC_FD, GET_SOLVER, SOLVER_OPTIONS,
%          GET_FIELD_FD_CPU, GET_FIELD_FD_CUDA, GEN_MASS_MATRIX_FD_CPU,
%          GEN_MASS_MATRIX_FD_CUDA, GEN_SOURCES.
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz and H. Dehghani 2018

%% check in/out

narginchk(2,4);
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
            error('Bad 3rd argument value. Text expected. Please see help for details on how to use this function.')
        end
        
        if isstruct(varargin{2})
            OPTIONS = varargin{2};
        else
            error('Bad 4th argument value. Structure expected. Please see help for details on how to use this function.')
        end
    else
        error('Bad arguments. Please see the help for details on how to use this function.')
    end
end

%% If not a workspace variable, load mesh
if ~isstruct(mesh)
    mesh = load_mesh(mesh);
end

%% check mesh type

if strcmp(mesh.type,'stnd')
    % single wavelength
    if nargout == 2
        [data,varargout{1}] = femdata_stnd_FD(mesh,frequency,solver,OPTIONS);
    else
        data = femdata_stnd_FD(mesh,frequency,solver,OPTIONS);
    end
elseif strcmp(mesh.type,'spec')
    % spectral mesh (all mesh wavelength)
    % See 'femdata_spec_FD' help to learn how to use this function for just some of the mesh wavelengths
    if nargout == 2
        [data,varargout{1}] = femdata_spec_FD(mesh,frequency,solver,OPTIONS);
    else
        data = femdata_spec_FD(mesh,frequency,solver,OPTIONS);
    end
else
    error(['Bad mesh type: ''' mesh.type '''. This function handles ''stnd'' and ''spec'' types only.'])
end

end

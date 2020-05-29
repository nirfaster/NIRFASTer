function [data] = femdata_TR(mesh, t, dt, varargin)
% FEMDATA_TR Calculates time-resolved photon fluence rate for a FEM mesh in
%   the NIRFAST format. This function is a wrapper to call
%   mesh-type-specific femdata functions. This function handles mesh types
%   'stnd' (single wavlength) and 'spec' (spectral) only.
% 
%   [DATA] = FEMDATA_TR(MESH, T, DT)
%   [DATA] = FEMDATA_TR(MESH, T, DT, 'field')
%   [DATA] = FEMDATA_TR(MESH, T, DT, 'field', SOLVER)
%   [DATA] = FEMDATA_TR(MESH, T, DT, 'field', OPTIONS)
%   [DATA] = FEMDATA_TR(MESH, T, DT, 'field', SOLVER, OPTIONS)
%   [DATA] = FEMDATA_TR(MESH, [], [], NUM_DT, N_MEAN_TIME)
%   [DATA] = FEMDATA_TR(MESH, [], [], NUM_DT, N_MEAN_TIME, SOLVER)
%   [DATA] = FEMDATA_TR(MESH, [], [], NUM_DT, N_MEAN_TIME, SOLVER, OPTIONS)
% 
% Please refere to helf of 'femdata_stnd_TR' (single wavelength mesh) or
% 'femdata_spec_TR' (spectral mesh) for details on the function syntax.
% 
% See also FEMDATA_STND_TR, FEMDATA_SPEC_TR.
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz and H. Dehghani 2018

%% check in/out

narginchk(3,7);
nargoutchk(0,1);


%% serve the variable input 

% is the field at mesh nodes should be returned 
isField = false;
% used for variable sampling related to source-detector distance, mean
% optical properties, etc.
num_dt = 0; % number of time bins
N_mean_time = 0; % calculate up to N_mean_time times the mean time of flight
% default
solver = get_solver;
OPTIONS = solver_options;

if ~isempty(varargin)
    if length(varargin) >= 1
        if ischar(varargin{1}) || isstring(varargin{1})
            if strcmp(varargin{1}, 'field')
                isField = true;
            else
                warning(['Unknown option ''' varargin{1} '''. Spatial distributions of time-resolved field will not be returned. Please see the help.'])
            end
        else
            num_dt = varargin{1};
            if ~(num_dt >= 1)
                error('Bad number of time bins. Should be >=1. Please see the help for details on how to use this function.')
            end
            num_dt = round(num_dt);
        end
    end
    if length(varargin) >= 2
        if ischar(varargin{2}) || isstring(varargin{2})
            % user specified, sanity check
            solver = get_solver(varargin{2});
        elseif ~isstruct(varargin{2}) && ~isempty(varargin{2})
            N_mean_time = varargin{2};
            if N_mean_time <= 0
                error('Bad number of mean time of flights. Should be >0. Please see the help for details on how to use this function.')
            end
        end

    end
    if length(varargin) >= 3
        if ischar(varargin{3}) || isstring(varargin{3})
            % user specified, sanity check
            solver = get_solver(varargin{3});
        elseif isstruct(varargin{3})
            OPTIONS = varargin{3};
        end
    end
    if length(varargin) == 4
        if isstruct(varargin{4})
            OPTIONS = varargin{4};
        end
    end
    
    if length(varargin) > 4
        error('Bad arguments. Please see the help for details on how to use this function.')
    end
        
end

% set the field value to relay
if isField
    field = 'field';
else
    field = [];
end


%% If not a workspace variable, load mesh
if ~isstruct(mesh)
    mesh = load_mesh(mesh);
end

%% check mesh type

if strcmp(mesh.type,'stnd')
    % single wavelength
    if N_mean_time == 0
        data = femdata_stnd_TR(mesh, t, dt, field, solver, OPTIONS);
    else
        data = femdata_stnd_TR(mesh, [], [], num_dt, N_mean_time, solver, OPTIONS);
    end
elseif strcmp(mesh.type,'spec')
    % spectral mesh (all mesh wavelength)
    % See the GUIDELINE section in 'femdata_spec_FD' help to learn how to use just some of the mesh wavelengths
    if N_mean_time == 0
        data = femdata_spec_TR(mesh, t, dt, field, solver, OPTIONS);
    else
        data = femdata_spec_TR(mesh, [], [], num_dt, N_mean_time, solver, OPTIONS);
    end
else
    error(['Bad mesh type: ''' mesh.type '''. This function handles ''stnd'' and ''spec'' types only.'])
end

end

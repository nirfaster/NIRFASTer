function [data] = femdata_TR_moments(mesh, varargin)
% FEMDATA_TR_MOMENTS Calculates spatial distributions and boundary values
%   of statistical moments or Mellin transforms of time-resolved data for a
%   FEM mesh in the NIRFAST format. This function is a wrapper to call
%   mesh-type-specific femdata functions. This function handles mesh types
%   'stnd' (single wavlength) and 'spec' (spectral) only. 
% 
% SYNTAX: 
%   [DATA] = FEMDATA_TR_MOMENTS(MESH)
%   [DATA] = FEMDATA_TR_MOMENTS(MESH, MAX_ORDER)
%   [DATA] = FEMDATA_TR_MOMENTS(MESH, MAX_ORDER, 'field')
%   [DATA] = FEMDATA_TR_MOMENTS(MESH, MAX_ORDER, 'field', 'mellin')
%   [DATA] = FEMDATA_TR_MOMENTS(MESH, MAX_ORDER, 'field', 'mellin', SOLVER)
%   [DATA] = FEMDATA_TR_MOMENTS(MESH, MAX_ORDER, 'field', 'mellin', OPTIONS)
%   [DATA] = FEMDATA_TR_MOMENTS(MESH, MAX_ORDER, 'field', 'mellin', SOLVER, OPTIONS)
%   [DATA] = FEMDATA_TR_MOMENTS(MESH, [], [], [], [], [])
% 
% Please refere to helf of 'femdata_stnd_TR_moments' (single wavelength
% mesh) or 'femdata_spec_TR_moments' (spectral mesh) for details on the
% function syntax. 
% 
% See also FEMDATA_STND_TR_MOMENTS, FEMDATA_SPEC_TR_MOMENTS.
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

% set the 'field' value to relay
if isField
    field = 'field';
else
    field = [];
end

% set the 'mellin' value to relay
if isMellin
    mellin = 'mellin';
else
    mellin = [];
end


%% If not a workspace variable, load mesh
if ~isstruct(mesh)
    mesh = load_mesh(mesh);
end

%% check mesh type

if strcmp(mesh.type,'stnd')
    % single wavelength
    data = femdata_stnd_TR_moments(mesh, max_order, field, mellin, solver, OPTIONS);
elseif strcmp(mesh.type,'spec')
    % spectral mesh (all mesh wavelength)
    % See the GUIDELINE section in 'femdata_spec_FD' help to learn how to use just some of the mesh wavelengths
    data = femdata_spec_TR_moments(mesh, max_order, field, mellin, solver, OPTIONS);
else
    error(['Bad mesh type: ''' mesh.type '''. This function handles ''stnd'' and ''spec'' types only.'])
end

end

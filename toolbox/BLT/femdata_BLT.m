function [data,mesh]=femdata_BLT(mesh, varargin)

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
            error('Bad 2nd argument value. Text or structure expected. Please see help for details on how to use this function.')
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

%% check mesh type
if strcmp(mesh.type,'stnd')
    % single wavelength
    if nargout == 2
        [data,varargout{1}] = femdata_stnd_BLT(mesh,solver,OPTIONS);
    else
        data = femdata_stnd_BLT(mesh,solver,OPTIONS);
    end
elseif strcmp(mesh.type,'spec')
    % spectral mesh (all mesh wavelength)
    % See 'femdata_spec_FD' help to learn how to use this function for just some of the mesh wavelengths
    if nargout == 2
        [data,varargout{1}] = femdata_spec_BLT(mesh,solver,OPTIONS);
    else
        data = femdata_spec_BLT(mesh,OPTIONS);
    end
else
    error(['Bad mesh type: ''' mesh.type '''. This function handles ''stnd'' and ''spec'' types only.'])
end

end

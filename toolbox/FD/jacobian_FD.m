function [J,data] = jacobian_FD(mesh,varargin)
% JACOBIAN_FD Calculates spatial distributions of sensitivity of
%   field registerd on the mesh boundary to changes of optical properties
%   per mesh node. This function relays parameters and calls specialized
%   functions based on the mesh type.
%  
% SYNTAX:
%  [J,DATA] = JACOBIAN_FD(MESH)
%  [J,DATA] = JACOBIAN_FD(MESH,FREQUENCY)
%  [J,DATA] = JACOBIAN_FD(MESH,FREQUENCY,SECOND_MESH_BASIS)
%  [J,DATA] = JACOBIAN_FD(MESH,FREQUENCY,SECOND_MESH_BASIS,SOLVER)
%  [J,DATA] = JACOBIAN_FD(MESH,FREQUENCY,SECOND_MESH_BASIS,SOLVER,OPTIONS)
%  [J,DATA] = JACOBIAN_FD(MESH,FREQUENCY,SECOND_MESH_BASIS,SOLVER,OPTIONS,'all') 
%  [J,DATA] = JACOBIAN_FD(MESH,FREQUENCY,[],[],[],[])
%  [J] = JACOBIAN_FD(___)
% 
% [J,DATA] = JACOBIAN_FD(MESH) MESH is NIRFAST mesh structure. Based on the
%   'MESH.type' following procedures are called:
%    - 'stnd' --> 'jacobian_stnd_FD'
%    - 'spec' --> 'jacobian_spec_FD'
%   J is the jacobian for CW (continuous wave) type sources. Please see
%   help for 'jacobian_stnd_FD' and 'jacobian_spec_FD' for the J structure.
%   DATA is the structure as returned by 'femdata_FD' with additional field
%   'aphi'. 'DATA.phi' is the forward field 'DATA.aphi' represents the
%   adjoint field (detectors are new sources) and the boundary data are
%   present in 'DATA.complex'.
% 
% [J,DATA] = JACOBIAN_FD(MESH,FREQUENCY) Optional FREQUENCY is modulation
%   frequency of sources in Hz. 
% 
% [J,DATA] = JACOBIAN_FD(MESH,FREQUENCY,SECOND_MESH_BASIS)
%   Optional SECOND_MESH_BASIS is a mesh in NIRFAST format that will be
%   used to calculate the Jacobians. This mesh is usually a coarse mesh as
%   compared with the MESH. See e.g. 'pixel_basis' or 'second_mesh_basis'
%   for details how to create such mesh. 
% 
% [J,DATA] = JACOBIAN_FD(MESH,FREQUENCY,SECOND_MESH_BASIS,SOLVER) Optional,
%   user specified SOLVER. Type 'help get_solver' for how to use solvers.
% 
% [J,DATA] = JACOBIAN_FD(MESH,FREQUENCY,SECOND_MESH_BASIS,SOLVER,OPTIONS)
%   Optional OPTIONS is a structure that allows to control solvers
%   parameters. Default OPTIONS structure is returned by 'solver_options'
%   function. See 'help solver_options' on how to use this structure.
% 
% [J,DATA] = JACOBIAN_FD(MESH,FREQUENCY,SECOND_MESH_BASIS,SOLVER,OPTIONS,'all') 
%   The 'all' option has effect only if data are real (CW - continuous
%   wave). If specified for the CW data, the attenuation sensitivity to
%   scattering (the 'kappa') is returned as well.
% 
% [J,DATA] = JACOBIAN_FD(MESH,FREQUENCY,[],[],[],'all') Any of the
%   parameters: SECOND_MESH_BASIS, SOLVER and OPTIONS can be set to empty.
%   This will set a given parameter to its default value or no 'coarse'
%   mesh is used.
% 
% [J] = JACOBIAN_FD(___) Does not retur the forward data. Jacobians only.
% 
% See also JACOBIAN_STND_FD, JACOBIAN_SPEC_FD.
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

%% check mesh type

if strcmp(mesh.type,'stnd')
    % single wavelength
    if isAll
        [J, data] = jacobian_stnd_FD(mesh, frequency, second_mesh_basis, solver, OPTIONS, 'all');
    else
        [J, data] = jacobian_stnd_FD(mesh, frequency, second_mesh_basis, solver, OPTIONS);
    end
elseif strcmp(mesh.type,'spec')
    % spectral mesh (all mesh wavelength)
    % See 'jacobian_spec_FD' help to learn how to use this function for just some of the mesh wavelengths
    if isAll
        [J, data] = jacobian_spec_FD(mesh, frequency, second_mesh_basis, solver, OPTIONS, 'all');
    else
        [J, data] = jacobian_spec_FD(mesh, frequency, second_mesh_basis, solver, OPTIONS);
    end
else
    error(['Unknown mesh type: ''' mesh.type '''.'])
end

end


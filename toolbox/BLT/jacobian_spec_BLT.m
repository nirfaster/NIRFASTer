function [J, varargout] = jacobian_spec_BLT(mesh, varargin)
%  Same format as femdata_stnd_BLT
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz and H. Dehghani 2018

%% check in/out
narginchk(1,3);
nargoutchk(0,2);

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
    elseif length(varargin) == 3
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
        if numel(varargin{3}) == 1
            isField = logical(varargin{3});
        else
            error('Bad 4th argument value. Scalar expected. Please see help for details on how to use this function.')
        end
    else
        error('Bad arguments. Please see the help for details on how to use this function.')
    end
end

%% check solver

solver = get_solver(solver);

%% If not a workspace variable, load mesh
if ~isstruct(mesh)
    mesh = load_mesh(mesh);
end

% check for spectral mesh and wavelengths field
if ~strcmp(mesh.type,'spec')
    warning(['Mesh type is ''' mesh.type '''. ''spec'' expected. This might give unexpected results.'])
end

if ~isfield(mesh,'wv')
    error('The mesh is missing ''wv'' field specifying wavelengths.');
end

%% Now calculate sources vector
qvec = conj(gen_source_adjoint_blt(mesh));

%% LOOP THROUGH WAVELENGTHS

for ind_wv = 1:length(mesh.wv)
    
    % get optical properties at nodes using spectra
    [mesh.mua, mesh.mus, mesh.kappa] = calc_mua_mus(mesh, mesh.wv(ind_wv));
    
    % make FEM matrices
    if isCUDA
        [i_index, j_index, value] = gen_mass_matrix_FD_CUDA(mesh,0);
    else
        [i_index, j_index, value] = gen_mass_matrix_FD_CPU(mesh,0);
    end
    
    % need conjugate of MASS matrix value
    value = conj(value);

    % get photon fluence rate
    % if GPU in use
    if strcmp(solver,solver_name_GPU)
        % check if we need to return the info structureas well
        if nargout == 2
            [phi, varargout{1}] = get_field_FD_CUDA(i_index, j_index, value, qvec, OPTIONS);
        else
            phi = get_field_FD_CUDA(i_index, j_index, value, qvec, OPTIONS);
        end
        % if no GPU in use, use CPU
    elseif strcmp(solver,solver_name_CPU)
        if nargout == 2
            [phi, varargout{1}] = get_field_FD_CPU(i_index, j_index, value, qvec, OPTIONS);
        else
            phi = get_field_FD_CPU(i_index, j_index, value, qvec, OPTIONS);
        end
    elseif strcmp(solver,solver_name_matlab_iterative)
        if nargout == 2
            [phi, varargout{1}] = get_field_FD_bicgstab_matlab(i_index, j_index, value, qvec, OPTIONS);
        else
            phi = get_field_FD_bicgstab_matlab(i_index, j_index, value, qvec, OPTIONS);
        end
    else
        % use MATLAB backslash
        % +1 as the i_index and j_index are zero-based and MATLAB uses 1-based indexing
        phi = full(sparse(i_index+1, j_index+1, value)\qvec);
        if nargout == 2
            varargout{1} = [];
        end
    end
    eval(['J.wv' num2str(mesh.wv(ind_wv)) ' = full(abs(phi)'');']);
    
end
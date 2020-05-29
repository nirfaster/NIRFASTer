function [data] = femdata_stnd_TR(mesh,t,dt,varargin)
% FEMDATA_STND_TR Returns time-resolved photon fluence rate for the
%   'standard' FEM mesh in NIRFAST format.
%   This function generates sparse matrices representing the FEM problem
%   (the time-resolved 'mass' mastrices and sparse sources vectors). Then,
%   it calls the time-resolved solver. 
%
% SYNTAX: 
%   [DATA] = FEMDATA_STND_TR(MESH, T, DT)
%   [DATA] = FEMDATA_STND_TR(MESH, T, DT, 'field')
%   [DATA] = FEMDATA_STND_TR(MESH, T, DT, 'field', SOLVER)
%   [DATA] = FEMDATA_STND_TR(MESH, T, DT, 'field', OPTIONS)
%   [DATA] = FEMDATA_STND_TR(MESH, T, DT, 'field', SOLVER, OPTIONS)
%   [DATA] = FEMDATA_STND_TR(MESH, [], [], NUM_DT, N_MEAN_TIME)
%   [DATA] = FEMDATA_STND_TR(MESH, [], [], NUM_DT, N_MEAN_TIME, SOLVER)
%   [DATA] = FEMDATA_STND_TR(MESH, [], [], NUM_DT, N_MEAN_TIME, SOLVER, OPTIONS)
%
%   [DATA] = FEMDATA_STND_TR(MESH, T, DT)
%    Returns the time-resolved DATA for sources and detectors as specified
%    in the MESH structure (MESH.link) within the time frame 0:DT:T-DT.
%    MESH - FEM mesh in NIRFAST format. Standard mesh (single wavelength).
%    T - in s, obserwation time, e.g. 10e-9 seconds
%    DT - in s, time bin width, the T discretization step, e.g. 10e-9/128
%    seconds for 128 time steps within T.
%    The output stucture DATA has following entries:
%       - DATA.phi [optional] is 3D matrix (n by m by k) of the raw photon
%         fluence rate where n - number of mesh nodes, m - number of
%         sources, k - number of the time bins. This value is optional and 
%         is enabled by the optional input parameter 'field'.
%       - DATA.time is the time vector in seconds calculated as 0:DT:T-DT
%       - DATA.tpsf is the temporal point spread function (the
%         time-resolved curve) of size (p by k) where p represents
%         source-detector combinations as specified in the MESH.link and k
%         is the number of the time bins. 
%       - DATA.link is copy of the MESH.link showing available
%         source-detector combinations.
% 
% WARNING: If the MESH.link specifies that a given source-detecor pair is
%          OFF, the DATA.tpsf values at that pair will be NaN. 
% 
%   [DATA] = FEMDATA_STND_TR(MESH, T, DT, 'field') Also returns the raw
%    photon fluence rate (DATA.PHI) per mesh node at all sources and time
%    bins. Values other than 'field' are ignored.
% 
%   [DATA] = FEMDATA_STND_TR(MESH, T, DT, 'field', SOLVER) User specified
%    SOLVER is used. See help of the 'get_solvers' function for how to use
%    solvers.
% 
%   [DATA] = FEMDATA_STND_TR(MESH, T, DT, 'field', SOLVER, OPTIONS) OPTIONS
%    is a structure that allows to control solvers parameters. The OPTIONS
%    structure can have following optional entries:
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
%    Entries within the structure are optional. The default OPTIONS
%    structure is returned by 'solver_options' function.
%
%   [DATA] = FEMDATA_STND_TR(MESH, [], [], NUM_DT, N_MEAN_TIME) Introduces
%    variable width of the time bins. Usefull if the mesh has sources and
%    detectors at short and long separations at the same time. However, it
%    will be slower as the calculation routine is called multiple times.
%    The routine is evaluated for number of unique source-detector
%    distances rounded up to the nearest centimeter. The returned time
%    resolved curves DATA.tpsf has size of n by NUM_DT where n is the
%    number of all source-detector combinations as specified in MESH.link
%    and NUM_DT is the requested numnber of time bins. The DATA.time is a
%    matrix of the same size as DATA.tpsf and represent the variable sample
%    time. N_MEAN_TIME specfies for how many mean time of flight of photons
%    each curve should be calculated.
% 
%   [DATA] = FEMDATA_STND_TR(MESH, [], [], NUM_DT, N_MEAN_TIME, SOLVER)
% 
%   [DATA] = FEMDATA_STND_TR(MESH, [], [], NUM_DT, N_MEAN_TIME, SOLVER, OPTIONS)
%  
% See also GET_SOLVER, SOLVER_OPTIONS, FEMDATA_SPEC_TR, GET_FIELD_TR_CPU,
%          GET_FIELD_TR_CUDA, GEN_MASS_MATRIX_TR_CPU,
%          GEN_MASS_MATRIX_TR_CUDA, GEN_SOURCES.
% 
% Read: Arridge S.R. et al. A finite element approach for modeling photon
%       transport in tissue. Med. Phys., 20, 1993, p. 299–309
%       Arridge S.R. and M. Schweiger Photon-measurement density
%       functions. Part 2: Finite-element-method calculations. Applied
%       Optics, Vol. 34, No. 34, 1995, p. 8026-8037
%       Liebert A. et al. Evaluation of optical properties of highly
%       scattering media by moments of distributions of times of flight of
%       photons. Applied Optics, Vol. 42, No. 28, 2003, p. 5785-5792
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

%% BODY

% if the unstable first part of curves should be removed
isMakeItNice = 1;

% If not a workspace variable, load mesh
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

% evenly spaced time bins
if N_mean_time == 0
    
    % make time vector
    data.time = 0:dt:t-dt;
    % size of time dimension
    num_dt = length(data.time);
    
    % make FEM matrices
    if isCUDA
        if isfield(OPTIONS,'GPU')
            [i_index, j_index, value_1, value_2] = gen_mass_matrix_TR_CUDA(mesh,dt,OPTIONS.GPU);
        else
            [i_index, j_index, value_1, value_2] = gen_mass_matrix_TR_CUDA(mesh,dt);
        end
    else
        [i_index, j_index, value_1, value_2] = gen_mass_matrix_TR_CPU(mesh,dt);
    end
    
    % if GPU in use
    if strcmp(solver,solver_name_GPU)
        % get photon fluence rate
        [phi, ~] = get_field_TR_CUDA(i_index,j_index,value_1,value_2,qvec,t,dt,OPTIONS);
    % if no GPU in use, use CPU
    elseif strcmp(solver,solver_name_CPU)
        % get photon fluence rate
        [phi, ~] = get_field_TR_CPU(i_index,j_index,value_1,value_2,qvec,t,dt,OPTIONS);
    elseif strcmp(solver,solver_name_matlab_iterative)
        % get photon fluence rate
        [phi, ~] = get_field_TR_bicgstab_matlab(i_index, j_index, value_1, value_2, qvec, t, dt, OPTIONS);

    else
        % use MATLAB backslash
        % +1 as the i_index and j_index are zero-based and MATLAB uses 1-based indexing
        phi = zeros(size(qvec,1),size(qvec,2),num_dt); % nodal data per source per time point
        % loop through time bins
        for ind_time = 1:num_dt
            if ind_time == 1
                % A1 * PHI(t=0) =  QVEC
                phi(:,:,ind_time) = full(sparse(i_index+1, j_index+1, value_1)\qvec);
            else
                % A1 * PHI(t+DT) = -A2 * PHI(t)
                phi(:,:,ind_time) = full(sparse(i_index+1, j_index+1, value_1)\(-sparse(i_index+1, j_index+1, value_2)*squeeze(phi(:,:,ind_time-1))));
            end
        end
    end

    % add raw data if needed
    if isField
        data.phi = phi;
    end

    % get boundary time-resolved curves from spatial photon fluence rate
    data.tpsf = zeros(num_meas,num_dt);
    for ind = 1:num_dt
        data.tpsf(:,ind)= get_boundary_data(mesh,phi(:,:,ind));
    end
    % if you want the data without NaN for disabled pairs
    % data.tpsf = data.tpsf(~isnan(data.tpsf(:,1)),:);
    
    % set source-detector link
    data.link = mesh.link;

    % make it nice!
    if isMakeItNice
        % all source-detector distances 
        rsd = sqrt(sum((mesh.source.coord(mesh.link(:,1),:) - ...
              mesh.meas.coord(mesh.link(:,2),:)).^2,2));
        % what is the minimum time of flight from source to detector
        t_cut = rsd./mean(mesh.c);

        % go through s-d pairs
        for ind_pair = 1:size(data.tpsf,1)
            % do nothing if the source-detector pair is turned off in the MESH.link
            if ~isnan(data.tpsf(ind_pair,1))
                % set zeroes at time of flight shorter than minumum possible time of flight between source and detector 
                data.tpsf(ind_pair,data.time <= t_cut(ind_pair)) = 0;

                % Also, set zeroes for the unstable part of solution.
                % Take the maximum of curve and look when it started to be
                % positive befor it reached the maximu. All values before that
                % are unstable and in fact should be zero.
                begi = 1;
                [~, maxl] = max(data.tpsf(ind_pair,:));
                for ind_sample = maxl:-1:1
                    if data.tpsf(ind_pair,ind_sample) <= 0 
                        begi = ind_sample;
                        break;
                    end
                end
                data.tpsf(ind_pair,1:begi) = 0;
            end
        end
    end
%variable spaced time bins
else
    % mean speed of light within the model
    cm = mean(mesh.c);
    % mean absorption within the model
    mua = mean(mesh.mua);
    % mean diffusion coeffcient within the model
    kappa = mean(mesh.kappa);

    % prealocate the time-resolved curves
    data.time = zeros(num_meas,num_dt);
    data.tpsf = zeros(num_meas,num_dt);
        
    % get all source-detector distances
    rsd = sqrt(sum((mesh.source.coord(mesh.link(:,1),:) - mesh.meas.coord(mesh.link(:,2),:)).^2,2));
    % round to the nearst 10 millimeters
    rsd_round = round(rsd,-1);
    % make sure we do not have zero distances by setting them to 10 mm
    rsd_round(rsd_round == 0) = 10;
    % take unique values of rounded to nearest centimeter
    rsd_unique = unique(rsd_round);
      
    % for all unique source-detector distances rounded to the nearest centimeter
    for ind_rsd = 1:length(rsd_unique)

        % get the mean time of flight of photons for mean optical
        % properties within the model and the current source-detector
        % distance as shown in Liebert2003
        mean_time = (rsd_unique(ind_rsd)^2)./(2*cm*sqrt(kappa).*(rsd_unique(ind_rsd)*sqrt(mua) + sqrt(kappa)));
        % the current max calculation time
        T_ = N_mean_time * mean_time;
        % the current time bine width
        dt_ = T_/num_dt;
        
        % make FEM matrices
        if isCUDA
            if isfield(OPTIONS,'GPU')
                [i_index, j_index, value_1, value_2] = gen_mass_matrix_TR_CUDA(mesh,dt_,OPTIONS.GPU);
            else
                [i_index, j_index, value_1, value_2] = gen_mass_matrix_TR_CUDA(mesh,dt_);
            end
        else
            [i_index, j_index, value_1, value_2] = gen_mass_matrix_TR_CPU(mesh,dt_);
        end

        % if GPU in use
        if strcmp(solver,solver_name_GPU)
            % get photon fluence rate
            [phi, ~] = get_field_TR_CUDA(i_index,j_index,value_1,value_2,qvec,T_,dt_,OPTIONS);
        % if no GPU in use, use CPU
        elseif strcmp(solver,solver_name_CPU)
            % get photon fluence rate
            [phi, ~] = get_field_TR_CPU(i_index,j_index,value_1,value_2,qvec,T_,dt_,OPTIONS);
        elseif strcmp(solver,solver_name_matlab_iterative)
            % get photon fluence rate
            [phi, ~] = get_field_TR_bicgstab_matlab(i_index, j_index, value_1, value_2, qvec, T_, dt_, OPTIONS);
        else
            % use MATLAB backslash
            % +1 as the i_index and j_index are zero-based and MATLAB uses 1-based indexing
            phi = zeros(size(qvec,1),size(qvec,2),num_dt); % nodal data per source per time point
            % loop through time bins
            for ind_time = 1:num_dt
                if ind_time == 1
                    % A1 * PHI(t=0) =  QVEC
                    phi(:,:,ind_time) = full(sparse(i_index+1, j_index+1, value_1)\qvec);
                else
                    % A1 * PHI(t+DT) = -A2 * PHI(t)
                    phi(:,:,ind_time) = full(sparse(i_index+1, j_index+1, value_1)\(-sparse(i_index+1, j_index+1, value_2)*squeeze(phi(:,:,ind_time-1))));
                end
            end
        end
        
        % get the time-resolved curves
        link_copy = mesh.link;
        mask_current_sd_pairs = rsd_round==rsd_unique(ind_rsd);
        mesh.link = mesh.link(mask_current_sd_pairs,:);
        for ind = 1:num_dt
            data.tpsf(mask_current_sd_pairs,ind) = get_boundary_data(mesh,phi(:,:,ind));
        end
        mesh.link = link_copy;

        % make the current time vector
        data.time(mask_current_sd_pairs,:) = ones(sum(mask_current_sd_pairs),1) * (0:dt_:T_-dt_);
        
    end
    
    % make it nice!
    if isMakeItNice
        % what is the minimum time of flight from source to detector
        t_cut = rsd./cm;

        % go through s-d pairs
        for ind_pair = 1:size(data.tpsf,1)
            % do nothing if the source-detector pair is turned off in the MESH.link
            if ~isnan(data.tpsf(ind_pair,1))
                % set zeroes at time of flight shorter than minumum possible time of flight between source and detector 
                data.tpsf(ind_pair,data.time(ind_pair,:) <= t_cut(ind_pair)) = 0;

                % Also, set zeroes for the unstable part of solution.
                % Take the maximum of curve and look when it started to be
                % positive befor it reached the maximu. All values before that
                % are unstable and in fact should be zero.
                begi = 1;
                [~, maxl] = max(data.tpsf(ind_pair,:));
                for ind_sample = maxl:-1:1
                    if data.tpsf(ind_pair,ind_sample) <= 0 
                        begi = ind_sample;
                        break;
                    end
                end
                data.tpsf(ind_pair,1:begi) = 0;
            end
        end
    end
    
    % copy the mesh link
    data.link = mesh.link;
end


end

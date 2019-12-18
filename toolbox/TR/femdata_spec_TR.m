function [data] = femdata_spec_TR(mesh,t,dt,varargin)
% FEMDATA_SPEC_TR Returns time-resolved photon fluence rate for the
%   'spectral' FEM mesh in NIRFAST format. This function generates sparse
%   matrices representing the FEM problem (the time-resolved 'mass'
%   mastrices and sparse sources vectors). Then, it calls the time-resolved
%   solver at all wavelengths specified in the NIRFAST mesh structure.     
% 
% SYNTAX: 
%   [DATA] = FEMDATA_SPEC_TR(MESH, T, DT)
%   [DATA] = FEMDATA_SPEC_TR(MESH, T, DT, 'field')
%   [DATA] = FEMDATA_SPEC_TR(MESH, T, DT, 'field', SOLVER)
%   [DATA] = FEMDATA_SPEC_TR(MESH, T, DT, 'field', SOLVER, OPTIONS)
%   [DATA] = FEMDATA_SPEC_TR(MESH, [], [], NUM_DT, N_MEAN_TIME)
%   [DATA] = FEMDATA_SPEC_TR(MESH, [], [], NUM_DT, N_MEAN_TIME, SOLVER)
%   [DATA] = FEMDATA_SPEC_TR(MESH, [], [], NUM_DT, N_MEAN_TIME, SOLVER, OPTIONS)
% 
%   [DATA] = FEMDATA_SPEC_TR(MESH, T, DT)
%    Returns the time-resolved DATA for sources and detectors as specified
%    in the MESH structure (MESH.link) within the time frame 0:DT:T-DT at
%    all wavelengths as specified in the mesh structure MESH.wv. 
%    MESH - FEM mesh in NIRFAST format. Spectral mesh (multi wavelengths)
%    T - in s, obserwation time, e.g. 10e-9 seconds
%    DT - in s, time bin width, the T discretization step, e.g. 10e-9/128
%    seconds 
%    The output stucture DATA has following entries:
%       - DATA.phi - [optional] is 4D matrix (n by m by l by k) of the raw
%         photon fluence rate where n - number of mesh nodes, m - number of
%         sources, l - number of wavelength, k - number of the time bins.
%         This value is optional and is enabled by the optional input
%         parameter 'field'. 
%       - DATA.time - is the time vector in seconds calculated as 0:DT:T-DT
%       - DATA.tpsf - is the temporal point spread function (the
%         time-resolved curve) of size (p by l by k) where p represents
%         source-detector combinations as specified in the MESH.link, l is
%         wavelength as specified in MESH.wv and k is the number of the
%         time bins.
%       - DATA.link - is copy of the MESH.link showing available
%         source-detector combinations
%       - DATA.wv - is copy of the MESH.wv showing the wavelengths
% 
% WARNING: If the MESH.link specifies that a given source-detecor pair is
%          OFF on a sfecified wavelength, the DATA.tpsf values at that pair
%          and wavelength will be NaN. 
% 
%   [DATA] = FEMDATA_SPEC_TR(MESH, T, DT, 'field')
%    Also returns the raw photon fluence rate per mesh node at all sources
%    and time bins. Values other than 'field' are ignored.
% 
%   [DATA] = FEMDATA_SPEC_TR(MESH, T, DT, 'field', SOLVER) User specified
%    SOLVER is used. See help of the 'get_solvers' function for how to use
%    solvers.
% 
%   [DATA] = FEMDATA_SPEC_TR(MESH, T, DT, 'field', SOLVER, OPTIONS) OPTIONS
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
%    Entries within the structure are optional. The default OPTIONS
%    structure is returned by 'solver_options' function.
% 
%   [DATA] = FEMDATA_SPEC_TR(MESH, [], [], NUM_DT, N_MEAN_TIME)
%    Introduces variable width of the time bins. Usefull if the mesh has
%    sources and detectors at short and long separations at the same time.
%    However, it will be slower as the calculation routine is called
%    multiple times. The routine is evaluated for number of unique
%    source-detector distances rounded up to the nearest centimeter. The
%    returned time-resolved curves DATA.tpsf has size of p by l by NUM_DT
%    where p is the number of all source-detector combinations as specified
%    in MESH.link, l is wavelength as specified in MESH.wv and NUM_DT is
%    the requested numnber of time bins. The DATA.time is a matrix of the
%    same size as DATA.tpsf and represent the variable sample time.
%    N_MEAN_TIME specfies for how many mean time of flight of photons each
%    curve should be calculated.
% 
%   [DATA] = FEMDATA_SPEC_TR(MESH, [], [], NUM_DT, N_MEAN_TIME, SOLVER)
% 
%   [DATA] = FEMDATA_SPEC_TR(MESH, [], [], NUM_DT, N_MEAN_TIME, SOLVER, OPTIONS)
% 
% See also GET_SOLVER, SOLVER_OPTIONS, FEMDATA_STND_TR, GET_FIELD_TR_CPU,
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
        else
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

% check for spectral mesh and wavelengths field
if ~strcmp(mesh.type,'spec')
    warning('The mesh type is not ''spec''. This function was intended for ''spec'' mesh type. Make sure that the rest mesh structure fields are correct.');
end
if ~isfield(mesh,'wv')
    error('The mesh is missing ''wv'' field specifying wavelengths.');
end

%% Now calculate sources vector for those active in mesh.link

% real values of sources in time domain
qvec = real(gen_sources(mesh));

%% Calculate forward data
% per mesh node for all sources at all wavelengths and all time bins

% number of source-detector pairs
num_meas = size(mesh.link,1);
% number of wavelengths
nlambda = length(mesh.wv);

% evenly spaced time bins
if N_mean_time == 0

    % make time vector
    data.time = 0:dt:t-dt;
    % size of time dimension
    num_dt = length(data.time);
    % declare time-resolved curve
    data.tpsf = zeros(num_meas,nlambda,num_dt);
    % declare raw data if needed
    if isField
        data.phi = zeros(size(mesh.nodes,1),size(qvec,2),nlambda,num_dt);
    end
    
    % make it nice!
    if isMakeItNice
        % all source-detector distances 
        rsd = sqrt(sum((mesh.source.coord(mesh.link(:,1),:) - ...
              mesh.meas.coord(mesh.link(:,2),:)).^2,2));
        % what is the minimum time of flight from source to detector
        t_cut = rsd./mean(mesh.c);
    end
    
    % for all wavelengths
    for ind_wv = 1:nlambda
        
        % get optical properties from spectra
        [mesh.mua, mesh.mus, mesh.kappa] = calc_mua_mus(mesh, mesh.wv(ind_wv));
    
        % make FEM matrices
        if isCUDA
            [i_index, j_index, value_1, value_2] = gen_mass_matrix_TR_CUDA(mesh,dt);
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
        
        % assign raw data if needed
        if isField
            data.phi(:,:,ind_wv,:) = phi;
        end
    
        % get boundary time-resolved curves from spatial photon fluence rate
        link_copy = mesh.link;
        % leave in the link flags for the current wavelength only
        % just to calculate the boundary data as it is designed to work
        % on a single wavelength (size(mesh.link,2) equals 3)
        mesh.link = mesh.link(:,[1:2 ind_wv + 2]);
        for ind = 1:num_dt
            data.tpsf(:,ind_wv,ind)= get_boundary_data(mesh,phi(:,:,ind));
        end
        mesh.link = link_copy;

        % make it nice!
        if isMakeItNice
            % go through s-d pairs
            for ind_pair = 1:size(data.tpsf,1)
                % do nothing if the source-detector pair is turned off in the MESH.link
                if ~isnan(data.tpsf(ind_pair,ind_wv,1))
                    % set zeroes at time of flight shorter than minumum possible time of flight between source and detector 
                    data.tpsf(ind_pair,ind_wv,data.time <= t_cut(ind_pair)) = 0;

                    % Also, set zeroes for the unstable part of solution.
                    % Take the maximum of curve and look when it started to be
                    % positive befor it reached the maximu. All values before that
                    % are unstable and in fact should be zero.
                    begi = 1;
                    [~, maxl] = max(data.tpsf(ind_pair,ind_wv,:));
                    for ind_sample = maxl:-1:1
                        if data.tpsf(ind_pair,ind_wv,ind_sample) <= 0 
                            begi = ind_sample;
                            break;
                        end
                    end
                    data.tpsf(ind_pair,ind_wv,1:begi) = 0;
                end
            end
        end
    end
    
    % set source-detector link
    data.link = mesh.link;
    % copy eavelengths to data
    data.wv = mesh.wv;
%variable spaced time bins    
else
    
    % allocate time-resolved data
    data.time = zeros(num_meas,nlambda,num_dt);
    data.tpsf = zeros(num_meas,nlambda,num_dt);
    
    % mean speed of light within the model
    cm = mean(mesh.c);
    
    % get all source-detector distances
    rsd = sqrt(sum((mesh.source.coord(mesh.link(:,1),:) - mesh.meas.coord(mesh.link(:,2),:)).^2,2));
    % round to the nearst 10 millimeters
    rsd_round = round(rsd,-1);
    % make sure we do not have zero distances by setting them to 10 mm
    rsd_round(rsd_round == 0) = 10;
    % take unique values of rounded to nearest centimeter
    rsd_unique = unique(rsd_round);

    % for all wavelengths
    for ind_wv = 1:nlambda
        
        % get optical properties from spectra
        [mesh.mua, mesh.mus, mesh.kappa] = calc_mua_mus(mesh, mesh.wv(ind_wv));

        % mean absorption within the model
        mua = mean(mesh.mua);
        % mean diffusion coeffcient within the model
        kappa = mean(mesh.kappa);
        
        % for all unique source-detector distances rounded to the nearest centimeter
        for ind_rsd = 1:length(rsd_unique)
            % get the mean time of flight of photons for mean optical
            % properties within the model and the current source-detector
            % distance as shown in Liebert2003
            mean_time = (rsd_unique(ind_rsd)^2)./(2*cm*sqrt(kappa).*(rsd_unique(ind_rsd)*sqrt(mua) + sqrt(kappa)));
            % the current max calculation time
            T_ = N_mean_time*mean_time;
            % the current time bine width
            dt_ = T_/num_dt;

                    
            % make FEM matrices
            if isCUDA
                [i_index, j_index, value_1, value_2] = gen_mass_matrix_TR_CUDA(mesh,dt_);
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
            % get mask of current pairs
            mask_current_sd_pairs = rsd_round==rsd_unique(ind_rsd);
            % set link to current pairs only
            mesh.link = mesh.link(mask_current_sd_pairs,:);
            % leave in the link flags for the current wavelength only
            % just to calculate the boundary data as it is designed to work
            % on a single wavelength (size(mesh.link,2) equals 3)
            mesh.link = mesh.link(:,[1:2 ind_wv + 2]);
            for ind = 1:num_dt
                data.tpsf(mask_current_sd_pairs,ind_wv,ind) = get_boundary_data(mesh,phi(:,:,ind));
            end
            mesh.link = link_copy;
            
            % make the current time vector
            data.time(mask_current_sd_pairs,ind_wv,:) = ones(sum(mask_current_sd_pairs),1) * (0:dt_:T_-dt_);

        end
        
        % make it nice!
        if isMakeItNice
            % what is the minimum time of flight from source to detector
            t_cut = rsd./cm;

            % go through s-d pairs
            for ind_pair = 1:size(data.tpsf,1)
                % do nothing if the source-detector pair is turned off in the MESH.link
                if ~isnan(data.tpsf(ind_pair,ind_wv,1))
                    % set zeroes at time of flight shorter than minumum possible time of flight between source and detector 
                    data.tpsf(ind_pair,ind_wv,data.time(ind_pair,ind_wv,:) <= t_cut(ind_pair)) = 0;

                    % Also, set zeroes for the unstable part of solution.
                    % Take the maximum of curve and look when it started to be
                    % positive befor it reached the maximu. All values before that
                    % are unstable and in fact should be zero.
                    begi = 1;
                    [~, maxl] = max(data.tpsf(ind_pair,ind_wv,:));
                    for ind_sample = maxl:-1:1
                        if data.tpsf(ind_pair,ind_wv,ind_sample) <= 0 
                            begi = ind_sample;
                            break;
                        end
                    end
                    data.tpsf(ind_pair,ind_wv,1:begi) = 0;
                end
            end
        end
    end
    
    % copy source-detectors link to data
    data.link = mesh.link;
    % copy wavelengths to data
    data.wv = mesh.wv;
   
end
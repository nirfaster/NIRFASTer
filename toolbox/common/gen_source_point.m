function qvec = gen_source_point(mesh,source,varargin)
% GEN_SOURCE_POINT Calculates FEM sources vector for those active in mesh.link.
% 
% [QVEC] = GEN_SOURCE_POINT(MESH, SOURCE) MESH is NIRFAST mesh
%   structure. SOURCES is a cartesian coordinate (2D or 3D in mm) of a
%   source or the source logical mask (M length, where M - number of MESH
%   sources). When SOURCE is logical, it is used for logical indexing of
%   the source within the MESH.source.coord coordinates matrix. QVEC is a
%   sparse, complex vector of size N by 1 of initial photon fluence rate at
%   N nodes. Since the point source is located within a MESH element, the
%   QVEC will have 3 or 4 nonzero entries depending on the MESH 2D or 3D
%   dimmension. Spatially integrated photon fluence rate of the source is
%   eual to 1 + 1j*eps.
% 
% See also MYTSEARCHN, GEN_SOURCES, GEN_SOURCE, FEMDATA_STND_FD. 
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz and H. Dehghani 2018

%% check in/out

narginchk(2,2);
nargoutchk(1,1);


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


%% check the source variable type

% Allocate memory
qvec = spalloc(size(mesh.nodes,1),1,mesh.dimension+1);

% cecki if logical indexing
if islogical(source)
    % check if only one source
    if sum(source) > 1
        error('Please specify only one source. Check the input SOURCE argument.');
    end
else
    % checki if coordinates from mesh.source.coord and try to make logical indexing
    mask_source = sum(mesh.source.coord - source,2) == 0;
    % if at least one match
    if sum(mask_source) == 1
        % make the logical index
        source = mask_source;
    elseif sum(mask_source) > 1
        % check if only one source
        error('Please specify only one source. Check the input SOURCE argument.');
    end
end

%% calculate or load the integration functions
% at this point, if it's not logical it's not present in the mesh.source.coord
if ~islogical(source)
    % find elements where the source point belongs and calculate interpolation function
    [ind,int_func] = mytsearchn(mesh,source(:,1:mesh.dimension));
else
    % check if we already have the source integration functions
    if isfield(mesh.source,'int_func')
        ind = mesh.source.int_func(source,1);
        int_func = mesh.source.int_func(source,2:end);
    else
        % find elements where the source point belongs and calculate interpolation function
        % use the logical indexing of mesh.source.coord
        [ind,int_func] = mytsearchn(mesh,mesh.source.coord(source,1:mesh.dimension));
    end
end

%% assign the result
% assign the sparse FEM source (initial conditions)
qvec(mesh.elements(ind,:)) = int_func .* (1 + 1j*eps);



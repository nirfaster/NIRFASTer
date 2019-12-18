function fn = add_extension(infn, ext)
% ADD_EXTENSION Adds an extension to a file if needed (checks if the
%   requested extension already exists). 
% 
% [FN] = ADD_EXTENSION(INFN, EXT) Adds the extension EXT to the filename
%   INFN if needed (checks if the requested extension EXT already exists at
%   the end of the INFN). The EXT may or may not have the leading extension
%   dot. It returns the full file name FN. E.g. INFN = '/dir1/dir2/file'
%   and EXT = 'xxx' or EXT = '.xxx' gives FN = '/dir1/dir2/file.xxx'; INFN
%   = 'file' and EXT = 'xxx' or EXT = '.xxx' gives FN = 'file.xxx'.
% 
% See also REMOVE_EXTENSION. 
% 
%   Part of NIRFAST package.

%% check in/out

narginchk(2,2);
nargoutchk(1,1);

%% BODY

% check if we need the extension dot
if ~isempty(ext)
    if ~strcmp(ext(1),'.')
        ext = strcat('.',ext);
    end
end

[~, fooext] = remove_extension(infn);

if ~strcmp(fooext,ext)
    fn = strcat(infn,ext);
else
    fn = infn;
end

end
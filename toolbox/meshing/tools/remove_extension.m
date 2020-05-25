function [fn, extension] = remove_extension(fname)
% REMOVE_EXTENSION Removes a file extension from the string/char array.
% 
% [FN,EXTENSION] = REMOVE_EXTENSION(FNAME) Removes a file extension from
%   the filename FNAME. It returns the file root FN and the extension
%   itself. E.g. FNAME = '/dir1/dir2/file.xxx' gives FN = '/dir1/dir2/file'
%   and EXTENSION = '.xxx'; FNAME = 'file.xxx' gives FN = 'file' and
%   EXTENSION = '.xxx'; FNAME = 'file' gives FN = 'file' and EXTENSION is
%   empty.
% 
% See also ADD_EXTENSION. 
% 
%   Part of NIRFAST package.

%% check in/out

narginchk(1,1);
nargoutchk(1,2);

%% BODY

[fdir, fname, extension] = fileparts(fname);
fn = fullfile(fdir, fname);

end

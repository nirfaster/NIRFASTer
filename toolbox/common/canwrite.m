function flag = canwrite(loc)
% CANWRITE Determines if there is write access to the directory.
% 
% SYNTAX:
%  CANWRITE(LOC)
%  [FLAG] = CANWRITE(LOC)
% 
% [FLAG] = CANWRITE(LOC) Determines if there is write access to the
%   directory LOC. If LOC is a full file path e.g. '/dir1/dir2/file.xxx'
%   the file path '/dir1/dir2' is checked. If LOC is empty, the current
%   directory is checked (as retrned by pwd). FLAG returns true (1) or
%   false (0). False is also returned if the location does not exist. If
%   the location exist, an empty file is created and immediately deleted
%   for testing.
%  
%   Part of NIRFAST package.

%% check in/out

narginchk(0,1);
nargoutchk(0,1);

%% BODY

% if check for current folder
if (nargin == 0) || isempty(loc)
    loc = pwd;
end


if isdir(loc)
    % if within MATLAB search path
    fooloc = loc;
else
    % if not in MATLAB search path
    % check if only a filename
    if isempty(fileparts(loc))
        % no path specified in the 'loc'
        fooloc = pwd;
    else
        % there is a path in the 'loc'
        fooloc = fileparts(loc);
    end
end

% get an unique name of a temporary file
[~,testFile,~] = fileparts(tempname);
loc = fooloc;


if exist(loc, 'dir') == 7
    % if directory exist, try to create the temporary text file there
    [fid, message] = fopen(fullfile(loc,testFile),'wt');
    % if can open the file and no system error message
    if (fid ~= -1) && isempty(message)
        % close and delete the file
        fclose(fid);
        delete(fullfile(loc,testFile));
        % we can write in 'loc'
        flag = true;
    else
        % we can't write in 'loc'
        flag = false;
    end
else
    % loc doesn't exist
    flag = false;
end

end
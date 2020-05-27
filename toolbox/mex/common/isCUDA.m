% ISCUDA Determines if there is a NVidia CUDA enabled GPU in the system.
%   Also, checks if the compute capability is at least 3.0. 
% 
% SYNTAX:
%   ISCUDA
%   [IS] = ISCUDA
%   [IS, INFO] = ISCUDA
%   [IS, INFO, MESSAGE] = ISCUDA
% 
%   ISCUDA Returns logical true or false:
%       0 - no available NVidia GPUs with the compute capability at least 3.0
%       1 - at least one NVidia GPU found and the compute capability is at
%       least 3.0 
% 
%   [IS, INFO] = ISCUDA Also returns the INFO structure. The structure size
%   is n by 1, where n is the number of NVidia CUDA enable GPUs found,
%   regardles of the compute capability.
%       - INFO.DeviceName - the device name, e.g. 'Quadro M2000M'
%       - INFO.ComputeCapability - the compute capability, e.g. 5
% 
%   [IS, INFO, MESSAGE] = ISCUDA Also returns optional MESSAGE string. The
%   MESSAGE can be empty if all went OK. The MESSAGE shows errors/messages
%   generated while identyfying GPU hardware.
%
% WARNING: By default, this function always returns 'false' on MAC systems
%          as the NIRFAST CUDA procedures are not supported on MAC. 
% 
%   MEX File function.
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

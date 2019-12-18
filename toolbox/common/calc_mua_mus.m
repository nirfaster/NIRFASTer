function [mua, mus, kappa, E] = calc_mua_mus(mesh,wv_array)
% CALC_MUA_MUS Given a spectral mesh, and specific wavelengths, calculates
%  optical properties at these wavelengths.
% 
% SYNTAX:
%  [MUA, MUS, KAPPA, E] = CALC_MUA_MUS(MESH)
%  [MUA, MUS, KAPPA, E] = CALC_MUA_MUS(MESH, WV_ARRAY)
%  [MUA, MUS, KAPPA] = CALC_MUA_MUS(___)
%  [MUA, MUS] = CALC_MUA_MUS(___)
%  [MUA] = CALC_MUA_MUS(___)
%  CALC_MUA_MUS(___)
% 
% [MUA, MUS, KAPPA, E] = CALC_MUA_MUS(MESH) Calcualte optical propertues at
%   all wavlelength as specified in 'MESH.wv' field.
%   Returns:
%    - MUA - in mm^-1, absorption coeficient at mesh nodes
%    - MUS - in mm^-1, reduced scattering coeficient at mesh nodes
%    - KAPPA - in mm, diffusion coeficient at mesh nodes (1/(3(MUA+MUS)))
%   and optional
%    - E - spectral 'extincton' coefficients. Copy of 'MESH.excoef' at
%          WV_ARRAY wavelengths only. Expressed in ln(10)*mm^-1mM^-1 (mM -
%          milli molar (molar - moles/liter)) for haemoglobins and in mm^-1
%          for water.
% 
% [MUA, MUS, KAPPA, E] = CALC_MUA_MUS(MESH, WV_ARRAY) WV_ARRAY specifies
%   requested wavelengths for optical properties calculation. WV_ARRAY
%   entries must exist in the 'MESH.wv' field.
% 
%   Part of NIRFAST package.

%% check in/out

narginchk(1,2);
nargoutchk(0,4);

% use mesh only
if nargin == 1
    wv_array = mesh.wv;
end
% check wv_array_size
if size(wv_array,1) > 1
    wv_array = wv_array';
end

%% BODY

% mask of 'wv_array' present in 'mesh.wv'
mask_wv = sum(mesh.wv == wv_array,2) > 0;
% check if all from 'wv_array' present in 'mesh.wv'
if sum(mask_wv) ~= numel(wv_array)
    error('Not all wavelength present in ''mesh.wv'' field')
end

% get extinction coefficients
E = mesh.excoef(mask_wv,:);
% calculate absorption coefficients in (nodes,wavelengths) size
mua = (E*mesh.conc')';

% for mus, wavelength must be in micrometers. Thus /1000.
mus = mesh.sa.*(wv_array/1000).^(-mesh.sp);
% diffusion
kappa = 1./(3*(mus+mua));

end
function [phi] = semi_infinite_FD(mua, musp, n, frequency, rsd, varargin)
% SEMI_INFINITE_FD Calculates photon fluence rate at any given point of
%   semin-infinite medium assuming global absorption, scattering and
%   refractive index.
% 
% SYNTAX:
%   [PHI] = SEMI_INFINITE_FD(MUA, MUSP, N, FREQUENCY, RSD)
%   [PHI] = SEMI_INFINITE_FD(MUA, MUSP, N, FREQUENCY, RSD, BOUNDARY)
%   [PHI] = SEMI_INFINITE_FD(MUA, MUSP, N, FREQUENCY, RSD, BOUNDARY, N_AIR)
%   [PHI] = SEMI_INFINITE_FD(MUA, MUSP, N, FREQUENCY, X, Z)
%   [PHI] = SEMI_INFINITE_FD(MUA, MUSP, N, FREQUENCY, X, Z, BOUNDARY)
%   [PHI] = SEMI_INFINITE_FD(MUA, MUSP, N, FREQUENCY, X, Z, BOUNDARY, N_AIR)
% 
% [PHI] = SEMI_INFINITE_FD(MUA, MUSP, N, FREQUENCY, RSD) Calculates the
%   photon fluence rate PHI at the boundary of semi-infinite medium and
%   air: 
%    MUA - in mm^-1, global absorption coeffcient
%    MUSP - in mm^-1, global reduced scattering coefficient
%    N - unitless, medium refractive indes
%    FREQUENCY - in Hz, a vector of m source modulation frequencies at which
%     the field will be calculated (FREQUENCY>=0)
%    RSD - in mm, a vector of n distances from the source at which the
%     boundary fiels will be calculated 
%    PHI - in mm^-2s^-1, n by m matrix of photon fluence rate at medium-air
%     boundry. n - length of RSD, m - length of FREQUENCY.
%   By default, the extrapolated boundary condition is used with the exact
%   calculation of boundary attenuation as numerical integrals of polarised
%   reflectances (see implementation for details).
%   Air refractive index is assumed to be 1.
%  
% [PHI] = SEMI_INFINITE_FD(MUA, MUSP, N, FREQUENCY, RSD, BOUNDARY)
%   Calculates the photon fluence rate PHI (mm^-2s^-1) at the surfe of
%   semi-infinite medium using extrapolated boundary condition and the
%   boundary attenuation calculated based on the BOUNDARY parameter as:
%    - 'exact' - (default) exact calculation of the internal
%        reflectance as numerical integrals of polarised reflectances (see
%        implementation for details). Marginally slover but most accurate.
%        Used by default. 
%    - 'approx' - Groenhuis polynamial approximation of the exact
%        solution (see implementation for details).
%    - 'robin' - robin boundary condition, internal reflectance
%        derived from Fresnel's law (see implementation for details).
%   Air refractive index is assumed to be 1.
% 
% [PHI] = SEMI_INFINITE_FD(MUA, MUSP, N, FREQUENCY, RSD, BOUNDARY, N_AIR)
%   N_AIR is the refractive index of the space outside the semi-infinite
%   medium. By default it is assumed as air with N_AIR=1.
% 
% [PHI] = SEMI_INFINITE_FD(MUA, MUSP, N, FREQUENCY, X, Z)
%   Calculates the photon fluence rate PHI (mm^-2s^-1) at any given point
%   of the semi-infinite medium. I.e. at the (X,Z) point. The semi-infinite
%   space is symetrical around perpendicular axis (Z-axis) going throug the
%   source. Thus, all unique PHI values can be calculated at coordinates
%   X>0 and Z>0. The assumed convention is that the medium semi-infinite
%   space has positive Z coordinates. X and Z are vectors of the same size
%   n, where the PHI (n by m) is calculated at n requestd points and m
%   frequencies. 
%   Returns boundary data if Z=0.
% 
% [PHI] = SEMI_INFINITE_FD(MUA, MUSP, N, FREQUENCY, X, Z, BOUNDARY)
% [PHI] = SEMI_INFINITE_FD(MUA, MUSP, N, FREQUENCY, X, Z, BOUNDARY, N_AIR)
%  
% See also BOUNDARY_ATTENUATION, SEMI_INFINITE_TR
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

%% error check and serve variable input

narginchk(5,8);
nargoutchk(0,1);

% default values
% default boundary
boundary = 'exact';
% default depth (the Z coordinate)
depth = 0;
% default air reefractive index
n_air = 1;


if ~isempty(varargin)
    if length(varargin) == 1
        if ischar(varargin{1}) || isstring(varargin{1})
            boundary = varargin{1};
        else
            depth = varargin{1};
        end
    elseif length(varargin) == 2
        if ischar(varargin{1}) || isstring(varargin{1})
            boundary = varargin{1};
            n_air = varargin{2};
        elseif ischar(varargin{2}) || isstring(varargin{2})
            depth = varargin{1};
            boundary = varargin{2};
        else
            error('Bad 5th/6th argument values. Please see help for details on how to use this function.')
        end
    elseif length(varargin) == 3
        if ischar(varargin{2}) || isstring(varargin{2})
            depth = varargin{1};
            boundary = varargin{2};
            n_air = varargin{3};
        else
            error('Bad 5th/6th/7th argument values. Please see help for details on how to use this function.')
        end
    else
        error('Bad arguments. Please see the help for details on how to use this function.')
    end
        
end

%% BODY

omega = 2*pi*frequency;

z = depth;

z0 = 1/(musp); % /mm  - mean free path
c0 = 299792458000; %/mm*s-1 - speed of light in air/vaccum
cm = c0/n; %/mm*s-1 - speed of light in medium
D = z0/3; % /mm - diffusion coeffcient

% result
phi = zeros(length(rsd),length(omega));

% get boundary attenuation
A = boundary_attenuation(n,n_air,boundary);


r1 = sqrt((z - z0).^2 + rsd.^2);
zb = 2*z0*A/3;
rb = sqrt((z + 2*zb + z0).^2 + (rsd).^2);

% Using analytical expression for semi-infinite model as shown in Durduran2010
for ind_w = 1:length(omega)
    k = sqrt((mua - 1j*omega(ind_w)/cm)./D);
    phi(:,ind_w) = (4*pi*D)^(-1) .* (exp(-k.*r1)./r1 - exp(-k.*rb)./rb);
end

end

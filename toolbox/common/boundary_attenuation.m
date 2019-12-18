function [A] = boundary_attenuation(n_incidence, varargin)
% BOUNDARY_ATTENUATION returns attenuation of photon fluence rate on a
%   boundary as caused by the mismatch in refractive indexes.
% 
% SYNTAX:
%  [A] = BOUNDARY_ATTENUATION(N_INCIDENCE)
%  [A] = BOUNDARY_ATTENUATION(N_INCIDENCE,METHOD)
%  [A] = BOUNDARY_ATTENUATION(N_INCIDENCE,N_TRANSMISSION)
%  [A] = BOUNDARY_ATTENUATION(N_INCIDENCE,N_TRANSMISSION,METHOD)
% 
%   [A] = BOUNDARY_ATTENUATION(N_INCIDENCE) Returns the boundary
%    attenuation A of diffusely scattered light that penetrates medium with
%    refractive index N_INCIDENCE. The medium borders with air described by
%    refractive index of 1.
%    By default, a Robin boundary condition is used where the attenuation A
%    is derived from Fresnel's law.
%       N_INCIDENCE - refractive indexes of the medium. It can be scalar,
%       vector or matrix of arbitrary size.
%   [A] = BOUNDARY_ATTENUATION(N_INCIDENCE,METHOD) The medium borders with
%    air described by refractive index of 1. However, the MATHOD is used to
%    change method of derivation the attenuation A. MATHOD can be specified
%    as follows:
%       'robin'  - (default) internal reflectance derived from Fresnel's law
%       'approx' - Groenhuis internal reflectance approximation (1.440n^-2
%                  + 0.710n^-1 + 0.668 + 0.00636n)
%       'exact'  - exact internal reflectance (integrals of polarised
%                  reflectances, etc. - see implementation for details) 
%   [A] = BOUNDARY_ATTENUATION(N_INCIDENCE,N_TRANSMISSION) The medium
%    borders with a medium described by refractive index of N_TRANSMISSION.
%    The default 'robin' method is used. N_TRANSMISSION should be a scalar.
%   [A] = BOUNDARY_ATTENUATION(N_INCIDENCE,N_TRANSMISSION,METHOD) The medium
%    borders with a medium described by refractive index of N_TRANSMISSION.
%    Te attenuaton derivation method as specified in METHOD is used.
% 
% Read: Durduran T. et al. Diffuse optics for tissue monitoring and
%       tomography. Rep. Prog. Phys., 73, 2010
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

%% error check

narginchk(1,3);
nargoutchk(1,1);

%% serve the variable input 

% used for variable sampling related to source-detector distance, mean
% optical properties, etc.
n_transmission = 1; % default refractive index of surrounding medium (air)
method = 'robin'; % default boundary attenuation derivation method

if ~isempty(varargin)
    if length(varargin) == 1
        if ischar(varargin{1})
            method = varargin{1};
        else
            n_transmission = varargin{1};
        end
    elseif length(varargin) == 2
        if ischar(varargin{1})
            error(['Bad 2nd argument value: ''' varargin{1} '''. A number expected. '...
                'Please see the help for details on how to use this function.'])
        end
        n_transmission = varargin{1};
        if ischar(varargin{2})
            method = varargin{2};
        else
            error(['Bad 3nd argument value: ''' varargin{1} '''. A character array expected. '...
                'Please see the help for details on how to use this function.'])
        end
    else
        error('Bad arguments. Please see the help for details on how to use this function.')
    end
end

%% BODY

if strcmp(method,'robin')
    A = robin(n_incidence,n_transmission);
elseif strcmp(method,'approx')
    A = approx(n_incidence,n_transmission);
elseif strcmp(method,'exact')
    A = exact(n_incidence,n_transmission);
else
    error(['Boundary option ''' num2str(method) ''' not recognized. '...
        'Please see the help for details on how to use this function.'])
end

end

function A = robin(n_incidence, n_transmission)
    n = n_incidence./n_transmission;
    R0 = ((n - 1).^2)./((n + 1).^2);
    % critical angle (total internal reflection (from tissue to air))
    theta_incidence = asin(1./n);
    A = (2./(1-R0) - 1 + abs(cos(theta_incidence)).^3)./(1 - abs(cos(theta_incidence)).^2);
end

function A = approx(n_incidence, n_transmission)
    n = n_incidence./n_transmission;
    Reff = -1.44*n.^-2 + 0.71*n.^-1 + 0.668 + 0.0636*n;
    A = (1 + Reff)./(1 - Reff);
end

function A = exact(n_incidence, n_transmission)

    % reflectance fluence rate at tissue-something boundary (from n_incidence to n_transmission)
    fun_fi = @(x,ni,nt) sin(2*x).*fresnel(x,ni,nt);
    % reflectance current density at tissue-something boundary (from n_incidence to n_transmission)
    fun_j = @(x,ni,nt) 3*sin(x).*(cos(x).^2).*fresnel(x,ni,nt);

    % declare output
    A = ones(size(n_incidence));
    % get unique n_incidences to avoid blind calculations
    n_incidence_unique = unique(n_incidence);
    % loop through unique refractive indexes
    for ind_n = 1:length(n_incidence_unique)
        % make a mask to know where to put the result
        mask_unique = n_incidence == n_incidence_unique(ind_n);
        
        % reflectance fluence rate at tissue-something boundary (from n_incidence to n_transmission)
        Rfi = integral(@(x)fun_fi(x,n_incidence_unique(ind_n),n_transmission),0,pi/2);
        % reflectance current density at tissue-something boundary (from n_incidence to n_transmission)
        Rj = integral(@(x)fun_j(x,n_incidence_unique(ind_n),n_transmission),0,pi/2);
        % reflectance at tissue-something boundary (from n_incidence to n_transmission)
        Reff = (Rfi + Rj)/(2 - Rfi + Rj);
        % extrapolation boundary (fluence rate = 0) (if times 2 * D)
        A(mask_unique) = (1 + Reff)/(1 - Reff);
    end
end

function Rf = fresnel(angle,n_incid,n_transmit)
    % s-polarised reflectance (perpendicular)
    Rs = abs((n_incid*cos(angle) - n_transmit*sqrt(1-(n_incid/n_transmit*sin(angle)).^2))./(n_incid*cos(angle) + n_transmit*sqrt(1-(n_incid/n_transmit*sin(angle)).^2))).^2;
%     Rf = Rs;
    % p-polarised reflectance (parallel)
    Rp = abs((n_incid*sqrt(1-(n_incid/n_transmit*sin(angle)).^2) - n_transmit*cos(angle))./(n_incid*sqrt(1-(n_incid/n_transmit*sin(angle)).^2) + n_transmit*cos(angle))).^2;
    % unpolarised reflectance
    Rf = (Rs + Rp)/2;
end
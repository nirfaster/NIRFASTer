function [phi] = semi_infinite_TR(mua, musp, n, rsd, T, dt, boundary)
% SEMI_INFINITE_TR Returns time-resolved boundary reflectance as collected
%  on semi infinite surface.
% 
% SYNTAX:
%   [PHI] = SEMI_INFINITE_TR(MUA, MUSP, N, RSD, T, DT)
%   [PHI] = SEMI_INFINITE_TR(MUA, MUSP, N, RSD, T, DT, BOUNDARY)
% 
% [PHI] = SEMI_INFINITE_TR(MUA, MUSP, N, RSD, T, DT, BOUNDARY) Calculates the time
%  resolved reflectance (see warning) at boundary of semi-infinite medium.
%  PHI - is a n by 2 matrix, where n is number of the time bins
%        (floor(T/DT) or length of T i T is a vector). First column is
%        time, second is the reflectance. 
%  MUA - absorption in mm^-1
%  MUS - reduced scatternig in mm^-1
%  N   - medium refractive index, unitless
%  RSD - source-detector distance on the surface in mm
%  T   - the observation time in seconds
%  DT  - time bin width in seconds
%  BOUNDARY is optional and could be:
%   PCB - partial current boundary condition
%    boundary = 'PCB-exact' - exact internal reflectance (integrals of polarised reflectances, etc. - see implementation)
%    boundary = 'PCB-approx' - Groenhuis internal reflectance approximation (?1.440n^-2 + 0.710n^-1 + 0.668 + 0.00636n)
%    boundary = 'PCB-Robin' - internal reflectance derived from Fresnel's law
%   EBC - extrapolated boundary condition
%    boundary = 'EBC-exact' - as for PCB
%    boundary = 'EBC-approx' - as for PCB
%    [default] boundary = 'EBC-Robin' - as for PCB
%   ZBC - zero boundary condition
%    boundary = 'ZBC'
% 
% WARNING:
%  REFLECTANCE it returns!
%  Do not compare directly with FEM time-resolved as NIRFAST FEM returne
%  photon fluence rate, not the reflectance!
% 
% See also BOUNDARY_ATTENUATION, SEMI_INFINITE_FD
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz 2018

%% check in/out

narginchk(6,7);
nargoutchk(0,1);


% default boundary condition
if nargin == 6
    boundary = 'EBC-Robin';
end

%% BODY

% z0 = 1/(mua + musp); % /mm  - mean free path
z0 = 1/(musp); % /mm  - mean free path
D = z0/3; % /mm - diffusion coeffcient
c0 = 299792458000; %/mm*s-1 - speed of light in air/vaccum
cm = c0/n; %/mm*s-1 - speed of light in medium

% TPSF time steps
if length(T) == 1
    steps = floor(T/dt);
    % TPSF calculation time
%     t = (1:steps-1)*dt; % /s
    % values for the centre of time bins
    t = (1:steps)*dt - dt/2; % /s
else
    steps = length(T);
    % TPSF calculation time
    t = T; % /s
    dt = t(2) - t(1);
end

phi = zeros(steps,2);
% phi(2:end,1) = t;
phi(:,1) = t;

if strcmp(boundary,'ZBC')
    phi(:,2) = (4*pi*D*cm)^(-3/2) * t.^(-5/2) * z0 .* exp(-mua*cm*t) .* exp(-(z0^2 + rsd^2)./(4*D*cm*t));
%     phi(:,2) = cm*(4*pi*D*cm*t).^(-3/2) .* exp(-mua*cm*t) .* (exp(-(4*z0^2 + rsd^2)./(4*D*cm*t)) - exp(-(rsd^2)./(4*D*cm*t)));
elseif strcmp(boundary,'EBC-exact')
    A = exact(n);
    zp = z0 + 4*A*D;
    phi(:,2) = 0.5 * (4*pi*D*cm)^(-3/2) * t.^(-5/2) .* exp(-mua*cm*t) .* (z0 * exp(-(z0^2 + rsd^2)./(4*D*cm*t)) + zp * exp(-(zp^2 + rsd^2)./(4*D*cm*t)));
elseif strcmp(boundary,'EBC-approx')
    A = approx(n);
    zp = z0 + 4*A*D;
    phi(:,2) = 0.5 * (4*pi*D*cm)^(-3/2) * t.^(-5/2) .* exp(-mua*cm*t) .* (z0 * exp(-(z0^2 + rsd^2)./(4*D*cm*t)) + zp * exp(-(zp^2 + rsd^2)./(4*D*cm*t)));
elseif strcmp(boundary,'EBC-Robin')
    A = robin(n);
    zp = z0 + 4*A*D;
    phi(:,2) = 0.5 * (4*pi*D*cm)^(-3/2) * t.^(-5/2) .* exp(-mua*cm*t) .* (z0 * exp(-(z0^2 + rsd^2)./(4*D*cm*t)) + zp * exp(-(zp^2 + rsd^2)./(4*D*cm*t)));
%     phi(:,2) = cm * (4*pi*D*cm*t).^(-3/2) .* exp(-mua*cm*t) .* (exp(-(z0^2 + rsd^2)./(4*D*cm*t)) - exp(-(zp^2 + rsd^2)./(4*D*cm*t)));
elseif strcmp(boundary,'PCB-exact')
    A = exact(n);
    a = z0 * A./(cm * t);
    b = 4/3 * A;
    Tpcb = 1./a .* (1 - sqrt(pi./(a*b)) .* exp(((1 + a).^2)./(a*b)) .* erfc((1+a)./(sqrt(a*b))));
    phi(:,2) = (4*pi*D*cm)^(-3/2) * t.^(-5/2) * z0 .* exp(-mua*cm*t) .* exp(-(z0^2 + rsd^2)./(4*D*cm*t)) .* Tpcb;
elseif strcmp(boundary,'PCB-approx')
    A = approx(n);
    a = z0 * A./(cm * t);
    b = 4/3 * A;
    Tpcb = 1./a .* (1 - sqrt(pi./(a*b)) .* exp(((1 + a).^2)./(a*b)) .* erfc((1+a)./(sqrt(a*b))));
    phi(:,2) = (4*pi*D*cm)^(-3/2) * t.^(-5/2) * z0 .* exp(-mua*cm*t) .* exp(-(z0^2 + rsd^2)./(4*D*cm*t)) .* Tpcb;
%     phi(1:end,2) = (4*pi*D*cm)^(-3/2) * t.^(-5/2) * z0 .* exp(-mua*cm*t) .* exp(-(z0^2 + rsd^2)./(4*D*cm*t)) .* Tpcb;
elseif strcmp(boundary,'PCB-Robin')
    A = robin(n);
    a = z0 * A./(cm * t);
    b = 4/3 * A;
    Tpcb = 1./a .* (1 - sqrt(pi./(a*b)) .* exp(((1 + a).^2)./(a*b)) .* erfc((1+a)./(sqrt(a*b))));
    phi(:,2) = (4*pi*D*cm)^(-3/2) * t.^(-5/2) * z0 .* exp(-mua*cm*t) .* exp(-(z0^2 + rsd^2)./(4*D*cm*t)) .* Tpcb;
else
    error(['Boundary option ''' num2str(boundary) ''' not recognized. See the function or type ''help diff_eq_semiinf'' for available options.'])
end

% scale to mm^-2s-1
phi(:,2) = phi(:,2) * dt;

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

function A = exact(n)

    % reflectance fluence rate at tissue-air boundary (from tissue to air)
    fun = @(x,ni,nt) sin(2*x).*fresnel(x,ni,nt);
    Rfi = integral(@(x)fun(x,n,1),0,pi/2);
    % reflectance current density at tissue-air boundary (from tissue to air)
    fun = @(x,ni,nt) 3*sin(x).*(cos(x).^2).*fresnel(x,ni,nt);
    Rj = integral(@(x)fun(x,n,1),0,pi/2);
    % reflectance at tissue-air boundary (from tissue to air)
    Reff = (Rfi + Rj)/(2 - Rfi + Rj);
    % extrapolation boundary (fluence rate = 0) (if times 2 * D)
    A = (1 + Reff)/(1 - Reff);
end

function A = approx(n)
    % (n_air = 1)
    Reff = -1.44*n^-2 + 0.71*n^-1 + 0.668 + 0.0636*n;
    A = (1 + Reff)/(1 - Reff);
end

function A = robin(n)
    % (n_air = 1)
    R0 = (n - 1)^2/(n + 1)^2;
    % critical angle (total internal reflection (from tissue to air))
    theta_incidence = asin(1/n);
    A = (2/(1-R0) - 1 + abs(cos(theta_incidence))^3)/(1 - abs(cos(theta_incidence))^2);
end
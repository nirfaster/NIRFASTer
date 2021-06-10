function [J] = build_jacobian_TR_element(mesh,data,varargin)
% BUILD_JACOBIAN_FD Calculates spatial distributions of sensitivity of
%   field registerd on the mesh boundary to changes of optical properties
%   per mesh node. It deliveres up to 4 distributions:
%   If data are real (CW - continuous wave):
%    - J_{mua}^{A}     - expressesd in mm, sensitivity of boundary
%                      attenuation (A) to absorption (mua [mm^-1]). This is
%                      the mean partial pathlength photons travel through
%                      the volume unit. The sum of this Jacobian gives mean
%                      photons pathlength in mm. 
%    - J_{kappa}^{A}   - (OPTIONAL) expressed in mm^-1, sensitivity of
%                      boundary attenuation (A) to diffusion coefficient
%                      'kappa' ((3*(mua+mus))^-1 [mm]). Returned if the
%                      optional parameter is set to 'all'. 
%   If data are complex (FD - frequency domain):
%    - J_{mua}^{A}     - af for CW.
%    - J_{kappa}^{A}   - af for CW but always returned.
%    - J_{mua}^{phi}   -  expressed in mm, sensitivity of boundary phase
%                      shift (phi), as related to the source phase, to
%                      absorption (mua [mm]).
%    - J_{kappa}^{phi} - expressed in mm^-1, sensitivity of boundary phase
%                      shift (phi) to diffusion coefficient 'kappa'
%                      ((3*(mua+mus))^-1 [mm]).
% 
% SYNTAX:
%  [J] = BUILD_JACOBIAN_FD(MESH, DATA)
%  [J] = BUILD_JACOBIAN_FD(MESH, DATA, 'all')
% 
% [J] = BUILD_JACOBIAN_FD(MESH, DATA) MESH is NIRFAST mesh
%   structure, DATA is a structure with following required entries: 
%    - DATA.PHI - (direct field) N by M matrix of photon fluence rate at N
%       nodes for M sources as used for source-detector pairs enabled in
%       'MESH.link'. Expressed in mm^-2s^-1 (3D) and mm^-1s^-1 (2D).
%    - DATA.APHI - (adjoint field) N by K matrix of photon fluence rate at
%       N nodes for K detectors used as sources. Detectors as used for
%       source-detector pairs enabled in 'MESH.link' field. Expressed in
%       mm^-2s^-1 (3D) and mm^-1s^-1 (2D).
%    - DATA.COMPLEX - (boundary data) P elements vector of photon fluence
%       rate for P source-detectors pairs as enabled in 'MESH.link'.
%       Expressed in mm^-2s^-1 (3D) and mm^-1s^-1 (2D).
%   Returned sensitivity sistribution has folowing fields depending on
%   sources type: CW - continuous wave (real numbers) or FD - frequency
%   domain (complex numbers):
%   If DATA.PHI and DATA.APHI are real numbers:
%    - J.complex is P (pairs) by N (nodes) matrix of spatial convolutions
%       of fields (DATA.PHI and DATA.APHI) expressed in photons (photon
%       packets to be exact) per second per millimeter (phot*s^-1*mm^-1).
%    - J.complete is P (pairs) by N (nodes) matrix of the J_{mua}^{A}
%       sensitivity expressed in millimeters.
%   If DATA.PHI and DATA.APHI are complex numbers:
%    - J.complex is P (pairs) by 2*N (2*nodes) matrix composed of two P by
%       N matrices [J_mus J_mua]. The J_mus is spatial convolution of
%       gradients of fields (DATA.PHI and DATA.APHI) expressed in photon
%       packets per second per millimeter cubed (phot*s^-1*mm^-3). J_mua is
%       the spatial convolution of fields (DATA.PHI and DATA.APHI) as in the
%       real numbers case but carried out using complex numbers arithmetics.
%    - J.complete is 2*P (pairs) by 2*N (2*nodes) matrix of: J_{mua}^{A}
%       (in mm, P by N), J_{kappa}^{A} (in mm^-1, P by N), J_{mua}^{phi}
%       (in mm, P by N) and J_{kappa}^{phi} (in mm^-1, P by N). J.complete
%       is organized as follows:
%       J.complete(1:2:end,:) = [J_{kappa}^{A} J_{mua}^{A}];
%       J.complete(2:2:end,:) = [J_{kappa}^{phi} J_{mua}^{phi}];
% 
% [J] = BUILD_JACOBIAN_FD(MESH, DATA, 'all') Has effect only if data are real
%   (CW - continuous wave). If specified for the CW data, the attenuation
%   sensitivity to scattering (the 'kappa') is returned as well. 
%    - J.complex is P (pairs) by 2*N (2*nodes) matrix composed of two P by
%      N matrices [J_mus J_mua]. The J_mus and J_mua are the same as in the
%      complex numbers case. However, calculations are carried out using
%      real numbers arithmetic.
%    - J.complete is P (pairs) by 2*N (2*nodes) matrix of: J_{kappa}^{A}
%      (in mm, P by N) and J_{mua}^{A} (in mm^-1, P by N). J.complete is
%      organized as follows: J.complete = [J_{kappa}^{A} J_{mua}^{A}]; 
%
% EQUATIONS:
%   Jacobians delivered using perturbation theory and Rytov approach.
%      - J_{mua}^{A} = dA/dmua = Phi_0(rd,rs)^-1 * int(G(r,rd)*Phi_0(r,rs),dV)
%        units 3D: [mm] = [(phot*s^-1*mm^-2)^-1] * ([mm^-1] * [(phot*s^-1*mm^-3)] * [mm^3])
%      - J_{kappa}^{A} = dA/dkappa = Phi_0(rd,rs)^-1 * int(grad(G(r,rd))*grad(Phi_0(r,rs)),dV)
%        units 3D: [mm^-1] = [(phot*s^-1*mm^-2)^-1] * ([mm^-1] * [mm^-1] * [mm^-1] * [(phot*s^-1*mm^-3)] * [mm^3])
%      - J_{mua}^{phi} = dphi/dmua = Im(J_{mua}^{A})
%        units 3D: [rad*mm]
%      - J_{kappa}^{phi} = dphi/dkapa = Im(J_{kappa}^{A})
%        units 3D: [rad*mm^-1]
%    Where:
%     - dA - change in light attenuation (no unit) measured on boundary
%     - dphi - change in phase (rad) measured on boundary
%     - dmua - change in absorption per node
%     - dkappa - change in 'absorption' (diffusion coefficient) per node
%     - rs - source position (Euclidean coordinate, in mm)
%     - rd - detector position (Euclidean coordinate, in mm)
%     - r - coordinate in Euclidean space (in mm)
%     - dV - elementary volume (or surface in 2D case)
%     - G(r,rd) - photon fluence rate (phot*mm^-2s^-1) for a source located
%       in the rd (detector position) - calculated with 'femdata_FD'
%     - Phi_0(r,rs) - photon fluence rate (phot*mm^-2s^-1) for a source
%       located in rs (source position) - calculated with 'femdata_FD'
%     - Phi_0(rd,rs) - photon fluence rate (phot*mm^-2s^-1) for a source
%       located in rs (source position) and detector in rd (detector
%       positiom) - calculated with 'get_boundary_data' as used in
%       'femdata_FD'.
%
% HINT
%   If you want to use the perturbation theory with the Born approach,
%   please use the J.complex and divide node-wise by the 'mesh.kappa'.
% 
% See also IntFG_tet4_CPU, IntgradFgradG_tet4_CPU, JACOBIAN_FD.
% 
%   Part of NIRFAST package.
%   S. Wojtkiewicz and H. Dehghani 2018

%% check in/out

narginchk(2,3);
nargoutchk(0,1);

%% get optional inputs, handle variable input

% default
isAll = false;

if ~isempty(varargin)
    if length(varargin) == 1
        if ischar(varargin{1}) || isstring(varargin{1})
            if strcmp(varargin{1},'all')
                isAll = true;
            end
        else
            error('Bad 3rd argument value. Text expected. Please see help for details on how to use this function.')
        end
    else
        error('Bad arguments. Please see the help for details on how to use this function.')
    end
end

%% BODY

% indexes of active pairs
index_pairs = find(logical(mesh.link(:,3)));
% [labels of active sources, index of source per active source-detector pair, ~]
[index_source_num, index_sources, ~] = sources_active_index(mesh);
% [labels of active detectors, index of detector per active source-detector pair, ~]
[index_detector_num, index_detectors, ~] = detectors_active_index(mesh);

% required as phi and aphi columns exist only for active and unique sources
% and detectors! Forward and adjoint fields are calculated for sources and
% detectors as calculated by 'sources_active_index' and
% 'detectors_active_index'
% indexes of columns in phi (sources permutatin matrix)
source_index = zeros(length(index_pairs),1);
% indexes of columns in aphi (detectors permutatin matrix)
detector_index = zeros(length(index_pairs),1);

% loop through active source-detectors pairs
for ind_active_pair = 1:length(index_pairs)
    % required as phi columns exist only for active and unique sources! 
    mask_source = index_source_num == index_sources(index_pairs(ind_active_pair));
    mask_detector = index_detector_num == index_detectors(index_pairs(ind_active_pair));
    source_index(ind_active_pair) = find(mask_source);
    detector_index(ind_active_pair) = find(mask_detector);
end


if isreal(data.phi) && isreal(data.aphi)
    % if the CW sensitivity to scatttering should be returned
    if isAll
        % spatial convolution for all enabled source-dtector pairs
        J_mua = -IntFG_tet4_CPU(mesh,data.phi,data.aphi,source_index,detector_index)';
        % spatial convolution of gradients of fields for all enabled source-dtector pairs
        J_mus = -IntgradFgradG_tet4_CPU(mesh,data.phi,data.aphi,source_index,detector_index)';

        % the full convolution
        J.complex = [J_mus J_mua];
    else
        % spatial convolution for all enabled source-dtector pairs
        J.complex = -IntFG_tet4_CPU(mesh,data.phi,data.aphi,source_index,detector_index)';
    end
    
    % Scale to the Rytow approximation (natural log of boundary data change)
    % J.complete = J.complex./data.complex;
    % We do not scale at this point. Conversely to the Cw or FD for which
    % the full Jacobian would be know at this point, we only know the
    % Jacobian when all its elements are computed by calling this function
    % several times.
else
    error('Both direct and adjoint fields should be real or complex at the same time.');
end

end
function [theta_0] = SkylineFunction( x, z, Lh, alpha, phi , flag)
% SkylineFunction.m is a function for calculating the skyline shielding
%   integration limit as a function of azimuth
% 
% Inputs:
%   x:          horizontal distance from ridgeline (m)
%   z:          vertical distance below surface (m)
%   Lh:         hillslope length in mapview (m)
%   alpha:      surface slope (radians) 
%   phi:        vector of incident radiation azimuthal directions (radians)
%   flag:       flag to determine "interior" (flag = 1) or "exterior" (flag = 0) catchment,
%                   as shown in Fig. 2 of DiBiase (2018)
%
% Outputs:
%   theta_0:    vector of skyline elevations (radians)
%
% Note: only valid for the range 0 ? phi ? pi
%
% Reference:
%    DiBiase, R.A., 2018. Short Communication: Increasing vertical attenuation length
%        of cosmogenic nuclide production on steep slopes negates topographic shielding
%        corrections for catchment erosion rates. Earth Surf. Dynam. Discuss., in review,
%        https://doi.org/10.5194/esurf-2018-48
%
% Written by Roman DiBiase rad22@psu.edu 8/15/2018
%

% Check to see that phi is everywhere between 0 and pi
if sum(phi<0|phi>pi)>0
    error('invalid range - phi must be in the interval [0:pi]')
end


if x==0
    theta_0 = zeros(size(phi)); % Assume no skyline shielding on ridge
    
else
    % Initialize skyline elevation angle array
    theta_0 = zeros(size(phi));

    % Calculate skyline elevation angle using Equation 15 of DiBiase (2018)
    
    if flag == 1
        theta_0(phi<(pi/2)) = abs(atan((x.*tan(alpha)+z).*cos(phi(phi<(pi/2)))./(2.*Lh-x)));
    end        
    theta_0(phi>=(pi/2)) = abs(atan(-tan(alpha).*cos(phi(phi>=(pi/2)))-(z./x).*cos(phi(phi>=(pi/2)))));
end

end
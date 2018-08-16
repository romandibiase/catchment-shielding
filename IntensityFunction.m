function [ normalized_flux ] = IntensityFunction(z, density, lambda_uni, alpha, phi, theta, m)
% IntensityFunction.m is a function for calculating the normalized radiation flux as a
%     function of slope position and depth below the surface
%
% Inputs
%   z:                  vertical distance below the surface (m)
%   density:            rock density (g cm^-2)
%   lambda_uni:         mass attenuation length for unidirectional incoming radiation (g cm^-2)
%   alpha:              surface slope (radians)
%   phi:                incident radiation azimuthal direction (radians)
%   theta:              incident radiation elevation angle (radians)
%   gamma:              apparent dip of surface in direction of phi (radians)
%   m:                  exponent of intensity function (typically assumed to be 2.3)
%
% Outputs
%   normalized flux:    total incident cosmic radiation flux normalized to horizontal surface with no shielding
%
% Reference:
%    DiBiase, R.A., 2018. Short Communication: Increasing vertical attenuation length
%        of cosmogenic nuclide production on steep slopes negates topographic shielding
%        corrections for catchment erosion rates. Earth Surf. Dynam. Discuss., in review,
%        https://doi.org/10.5194/esurf-2018-48
%
% Written by Roman DiBiase rad22@psu.edu 8/15/2018
%

% Calculate apparent dip, gamma
gamma = -atan(cos(phi).*tan(alpha));            % Equation 7 of DiBiase (2018)

% Calculate the mass distance traveled through rock
d = (density.*(z.*100).*cos(gamma))./sin(theta-gamma);  % Equation 8 of DiBiase (2018)

% Calculate normalized flux (integrand in Equation 3 of DiBiase (2018)
normalized_flux = ((m+1)./(2.*pi)).*((sin(theta).^m).*exp(-d./lambda_uni).*cos(theta));

end
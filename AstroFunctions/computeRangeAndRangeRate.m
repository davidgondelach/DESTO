function [rangeAndRangeRate] = computeRangeAndRangeRate(rSatECI,vSatECI,rSiteECI,vSiteECI)
% Compute range and ragne rate from satellite and observer ECI positions 
% and velocities
%
% Based on: rv2tradc.m by Vallado
% Reference: Vallado, 2013, Fundamentals of Astrodynamics, p.260, Algorithm 26

% Compute range vector from site to satellite
rhoECI  = rSatECI - rSiteECI;
drhoECI = vSatECI - vSiteECI;
% Compute range from site to satellite
rho     = sqrt( sum( rhoECI.^2, 1 ));
% Compute range rate from site to satellite
drho= dot(rhoECI,drhoECI)./rho;

rangeAndRangeRate = [rho;drho];

end


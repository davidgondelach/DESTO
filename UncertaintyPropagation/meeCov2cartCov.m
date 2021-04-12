function [posVel,posVelCov] = meeCov2cartCov(mee,meeCov,GM_kms,varargin)
% meeCov2cartCov - Convert covariance from modified equinoctial elements to
% Cartesian coordinates using unscented transform
%
%  Copyright (C) 2021 by David Gondelach
%
% This code is licensed under the GNU General Public License version 3.
%
% Author: David Gondelach
% Massachusetts Institute of Technology, Dept. of Aeronautics and Astronautics
% email: davidgondelach@gmail.com
% Aug 2020; Last revision: 31-Aug-2020

%------------- BEGIN CODE --------------

if nargin > 3
    Qmee2cart = varargin{1};
else
    Qmee2cart = zeros(6,6);
end

transformationFnc = @(xx) ep2cart(xx,GM_kms);
[posVel,posVelCov] = unscentedUncertaintyTransform(mee,transformationFnc,meeCov,Qmee2cart);

end

%------------- END OF CODE --------------
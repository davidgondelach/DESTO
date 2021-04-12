function [mee,meeCov] = cartCov2meeCov(cart,cartCov,GM_kms,varargin)
% cartCov2meeCov - Convert covariance from Cartesian coordinates to
% modified equinoctial elements using unscented transform
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
    Qcart2mee = varargin{1};
else
    Qcart2mee = zeros(6,6);
end

transformationFnc = @(xx) cart2ep_wrapL(xx,GM_kms);
[mee,meeCov] = unscentedUncertaintyTransform(cart,transformationFnc,cartCov,Qcart2mee);

end

%------------- END OF CODE --------------
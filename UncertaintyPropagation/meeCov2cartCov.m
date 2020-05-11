function [posVel,posVelCov] = meeCov2cartCov(mee,meeCov,GM_kms,varargin)

if nargin > 3
    Qmee2cart = varargin{1};
else
    Qmee2cart = zeros(6,6);
end

transformationFnc = @(xx) ep2cart(xx,GM_kms);
[posVel,posVelCov] = unscentedUncertaintyTransform(mee,transformationFnc,meeCov,Qmee2cart);

end


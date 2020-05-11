function [mee,meeCov] = cartCov2meeCov(cart,cartCov,GM_kms,varargin)

if nargin > 3
    Qcart2mee = varargin{1};
else
    Qcart2mee = zeros(6,6);
end

transformationFnc = @(xx) cart2ep_wrapL(xx,GM_kms);
[mee,meeCov] = unscentedUncertaintyTransform(cart,transformationFnc,cartCov,Qcart2mee);

end


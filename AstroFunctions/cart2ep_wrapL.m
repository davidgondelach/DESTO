function [EP] = cart2ep_wrapL(cart, mu)

EP = cart2ep(cart, mu);

% Make sure the difference in true longitude L of the sigma points wrt
% the nominal true longitude L0 is minimal (i.e. L-L0 <= pi)
% If the nominal true longitude L0 is close to pi (i.e. pi/2<L0 or L0<-pi/2) then wrap
% all L to [0,2pi] domain, so all difference in L <=pi (by default L is
% on [-pi,pi] domain).
if EP(6,1) > pi/2 || EP(6,1) < -pi/2
    EP(6,:) = wrapTo2Pi(EP(6,:));
end

end


% This function converts equinoctial parameters to eci position
% and velocity vectors 

function [cart] = ep2cart(EP, mu)

% input

%  mu     = gravitational constant (km**3/sec**2)
%  EP(1) = semilatus rectum of orbit (kilometers)
%  EP(2) = f equinoctial element
%  EP(3) = g equinoctial element
%  EP(4) = h equinoctial element
%  EP(5) = k equinoctial element
%  EP(6) = true longitude (radians)

% output

%  rr = eci position vector (kilometers)
%  vv = eci velocity vector (kilometers/second)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% unload equinoctial orbital elements

[m,n] = size(EP);

if m~=6
    error('Number of rows of MEE input must be 6');
end

cart = zeros(m,n);

for i = 1:n
    
clearvars -except i m n EP mu cart

p = EP(1,i);
f = EP(2,i);
g = EP(3,i);
h = EP(4,i);
k = EP(5,i);
l = EP(6,i);

sqrtmup = sqrt(mu / p);

cosl = cos(l);

sinl = sin(l);

q = 1 + f * cosl + g * sinl;

r = p / q;

alphasqrd = h^2 - k^2;

ssqrd = 1 + h^2 + k^2;

% compute eci position vector

rr	= r / ssqrd * ...
            [   cosl + alphasqrd * cosl + 2 * h * k * sinl;
                sinl - alphasqrd * sinl + 2 * h * k * cosl;
                2 * (h * sinl - k * cosl);
            ];
        
% compute eci velocity vector

vv	= sqrtmup  / ssqrd * ...
            [   - sinl - alphasqrd * sinl + 2 * h * k * cosl - g ...
                	+ 2 * f * h * k - alphasqrd * g;
                cosl - alphasqrd * cosl - 2 * h * k * sinl + f ...
                    - 2 * g * h * k - alphasqrd * f;
                2 * (h * cosl + k * sinl + f * h ...
                    + g * k)
            ];

cart(:,i) = [rr ; vv];

end
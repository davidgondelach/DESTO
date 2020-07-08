function [singlePosition] = extractSinglePosition(Xp,meas,svs,mu)
% - Return only position of single object without BC or ROM state

% Get object number from measurement
objectNumber = meas(1);

[cart] = ep2cart(Xp(svs*(objectNumber-1)+1:svs*(objectNumber-1)+6,:), mu);
% Get position
singlePosition = cart(1:3,:);

end
function [singleState] = extractSingleState(Xp,meas,svs)
% - Return only state of single object without BC or ROM state

% Get object number from measurement
objectNumber = meas(1);
% Get state
singleState = Xp(svs*(objectNumber-1)+1:svs*(objectNumber-1)+6,:);

end
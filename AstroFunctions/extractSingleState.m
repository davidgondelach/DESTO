function [singleState] = extractSingleState(Xp,objectNumber,svs)
% - Return only state of single object without BC or ROM state

singleState = Xp(svs*(objectNumber-1)+1:svs*(objectNumber-1)+6,:);

end
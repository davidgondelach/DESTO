function [yres] = getResidualMEE(ym,meas)
% ym = predicted measurement
% meas = measurement including object number

yres = meas(2:end)-ym;
yres(6:6:end) = wrapToPi(yres(6:6:end)); % Wrap difference in true longitude to [-pi,pi]

end


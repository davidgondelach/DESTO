function [yres] = getResidualDefault(ym,meas)
% ym = predicted measurement
% meas = measurement including object number

yres = meas(2:end)-ym;

end


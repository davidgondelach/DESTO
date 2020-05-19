function [yres] = getResidualRangeRangeRate(ym,meas)
% ym = predicted measurement: range and range rate
% meas = measurement: 1) object number, 2) range, 3) range rate, 4-9) Observer location

yres = meas(2:3)-ym;

end


function [reci, veci] = convertECEFtoJ2000(recef, vecef, jdate)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global EOPMat;
[ xp, yp, dut1, lod, ddpsi, ddeps, dat ] = computeEOP_Celestrak( EOPMat, jdate );

date = jed2date(jdate);
year = date(1);
mon = date(2);
day = date(3);
hr = date(4);
min = date(5);
sec = date(6);

timezone = 0;
[~, ~, jdut1, ~, ~, ~, ttt, ~, ~, ~, ~, ~, ~, ~, ~ ] ...
         = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );

eqeterms = 2; % use the kinematic eqe terms after 1997
[reci,veci,~] = ecef2eci  ( recef,vecef,zeros(3,1),ttt,jdut1,lod,xp,yp,eqeterms,ddpsi,ddeps );

end


function [plotTime,rho_sat,rho_rom,rho_jb2,rho_msise,SWdata,varargout] = getDensitiesRealRomJBMSISE(data,etROM,romStates,r,F_U,M_U,SOLdata,DTCdata,eopdata,SWmatDaily,SWmatDailyTIEGCM,varargin)

computeCov = nargin > 11;
if computeCov
    romCovs = varargin{1};
end

n = length(data.longitudes);
rho_rom = zeros(n,1);
rho_sat = zeros(n,1);
rho_jb2 = zeros(n,1);
rho_msise = zeros(n,1);
plotTime = zeros(n,1);
SWdata = zeros(n,4);
if computeCov
    rho_rom_std = zeros(n,1);
end

datetype = -1;
if isfield(data,'seconds') && isfield(data,'doys')
    datetype = 1;
elseif isfield(data,'dates')
    datetype = 2;
else
    error('Data does not contain valid date format');
end

for i=1:n
    
    lon = data.longitudes(i);
    slt = data.localtimes(i);
    lat = data.latitudes(i);
    alt = data.altitudes(i);
    
    if datetype == 1
        secondsGPS = data.seconds(i);
        secondsTDT = secondsGPS + 19 + 32.184; % TDT = TAI + 32.184, TAI = GPS time + 19
        [hrTDT,minTDT,secTDT] = sec2hms( secondsTDT );
        yr = data.years(i) + 2000;
        doy_ = data.doys(i);
        if hrTDT >= 24
            hrTDT = hrTDT - 24;
            doy_ = doy_+1;
        end
        dv = datevec(datenum(yr, 1, doy_));
        mth = dv(2);
        dy = dv(3);
        et  = cspice_str2et(strcat([num2str([yr mth dy hrTDT minTDT secTDT],'%d %d %d %d %d %.10f') 'TDT']));
    elseif datetype == 2
        dv = datevec(data.dates(i));
        et  = cspice_str2et(strcat([num2str([dv(1) dv(2) dv(3) dv(4) dv(5) dv(6)],'%d %d %d %d %d %.10f') 'UTC']));
    end
        
    
    % ROM density
    romState = interp1(etROM, romStates', et, 'linear','extrap');
    slt(slt>24) = slt(slt>24)-24; slt(slt<0) = slt(slt<0)+24;
    UhI = zeros(1,r);
    for j = 1:r
        UhI(:,j) = F_U{j}(slt,lat,alt);
    end
    MI = M_U(slt,lat,alt);
    rho_rom(i) = 10.^(UhI*romState'+MI');
    % ROM covariance
    if computeCov
        romCov = interp1(etROM, romCovs', et);
        Pyy = diag(UhI * diag(romCov) * UhI');
%         rho_rom_std(i) = 100*Pyy.^0.5*log(10); % percentage std
        rho_rom_std(i) = Pyy.^0.5*log(10)*rho_rom(i); % true std
    end
    
    % Satellite real density
    rho_sat(i) = data.densities(i);
    
    % JB2008 density
    jdatestrUTC    = cspice_et2utc( et, 'J', 12 );
    jdateUTC       = str2double(jdatestrUTC(4:end)); % Cut trailing 'JD ' off from string
    [yyUTC, mmUTC, ddUTC, hhUTC, mnmnUTC, ssUTC] = datevec(jdateUTC-1721058.5);
    doyUTC = day(datetime(yyUTC, mmUTC, ddUTC),'dayofyear');
    [MJD,GWRAS,SUN,F10,F10B,S10,S10B,XM10,XM10B,Y10,Y10B,DSTDTC] = computeJB2000SWinputs(yyUTC,doyUTC,hhUTC,mnmnUTC,ssUTC,SOLdata,DTCdata,eopdata);
    
    XLON = deg2rad(lon); % Lon
    SAT(1) = mod(GWRAS + XLON, 2*pi);
    SAT(2) = deg2rad(lat); % Lat
    SAT(3) = alt;
    [~,rho_jb2(i)] = JB2008(MJD,SUN,SAT,F10,F10B,S10,S10B,XM10,XM10B,Y10,Y10B,DSTDTC);
    
    % NRLMSISE-00 density
    [ f107A, f107, ap ] = computeSWnrlmsise( SWmatDaily, [], jdateUTC );
    UTsec = hhUTC*3600+mnmnUTC*60+ssUTC;
    [output] = nrlmsise(alt,lat,lon,yyUTC,doyUTC,UTsec,slt,f107A,f107,ap);
    rho_msise(i) = output.d(6);

    %     plotTime(i) = et;
    if datetype == 1
        plotTime(i) = data.doys(i) + data.seconds(i)/86400;
    else
        yearStartDateNum = datenum(yyUTC,1,1);
        plotTime(i) = data.dates(i) - yearStartDateNum + 1;
    end
    
    % Space weather data
    [ ~, f107, Kp ] = computeSWnrlmsise( SWmatDailyTIEGCM, [], jdateUTC, true );
    SWdata(i,1) = F10;
    SWdata(i,2) = DSTDTC;
    SWdata(i,3) = f107;
    SWdata(i,4) = Kp(2);
end

if computeCov
    varargout{1} = rho_rom_std;
end

end
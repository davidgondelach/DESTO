function [gpsMeas] = getGPSpositionMeasurements(gpsDataPath,selectedObjects,et0,etf)

% Planet satellites with GPS data
planetSkysats = struct;
planetSkysats(1).planetID = 's1';     planetSkysats(1).noradID = 39418; % SKYSAT 1 [MISLABEL]
planetSkysats(2).planetID = 's2';     planetSkysats(2).noradID = 40072; % SKYSAT 2 [MISLABEL]
planetSkysats(3).planetID = 's3';     planetSkysats(3).noradID = 41601; % SKYSAT C1
planetSkysats(4).planetID = 's4';     planetSkysats(4).noradID = 41773; % SKYSAT C2
planetSkysats(5).planetID = 's103';   planetSkysats(5).noradID = 41774; % SKYSAT C3
planetSkysats(6).planetID = 's104';   planetSkysats(6).noradID = 41771; % SKYSAT C4
planetSkysats(7).planetID = 's105';   planetSkysats(7).noradID = 41772; % SKYSAT C5
planetSkysats(8).planetID = 's106';   planetSkysats(8).noradID = 42992; % SKYSAT C6
planetSkysats(9).planetID = 's107';   planetSkysats(9).noradID = 42991; % SKYSAT C7
planetSkysats(10).planetID = 's108';  planetSkysats(10).noradID = 42990; % SKYSAT C8
planetSkysats(11).planetID = 's109';  planetSkysats(11).noradID = 42989; % SKYSAT C9
planetSkysats(12).planetID = 's110';  planetSkysats(12).noradID = 42988; % SKYSAT C10
planetSkysats(13).planetID = 's111';  planetSkysats(13).noradID = 42987; % SKYSAT C11
planetSkysats(14).planetID = 's112';  planetSkysats(14).noradID = 43797; % SKYSAT C12
planetSkysats(15).planetID = 's113';  planetSkysats(15).noradID = 43802; % SKYSAT C13

% Get GPS measurements
gpsMeas = [];
for i=1:length(selectedObjects)
    index = find([planetSkysats.noradID]==selectedObjects(i));
    
    % Get GPS data
    [gpsData] = getGPSdataFromJsonFileET(gpsDataPath,planetSkysats(index).planetID,planetSkysats(index).noradID,et0,etf);
    
    % Remove half of measurements closer than 5 seconds from each other
    smallerThan5sec = find(diff([gpsData(:).tET]) <= 5);
    gpsData(smallerThan5sec(1:2:end)) = [];
    % Remove half of measurements closer than 10 seconds from each other
    smallerThan10sec = find(diff([gpsData(:).tET]) <= 10);
    gpsData(smallerThan10sec(1:2:end)) = [];
    % Remove half of measurements closer than 20 seconds from each other
    smallerThan20sec = find(diff([gpsData(:).tET]) <= 20);
    gpsData(smallerThan20sec(1:2:end)) = [];
    % Remove half of measurements closer than 40 seconds from each other
    smallerThan40sec = find(diff([gpsData(:).tET]) <= 40);
    gpsData(smallerThan40sec(1:2:end)) = [];
    % Remove measurements closer than 40 seconds from each other
    smallerThan40sec = find(diff([gpsData(:).tET]) <= 40);
    gpsData(smallerThan40sec) = [];
    % Remove half of measurements closer than 80 seconds from each other
    smallerThan80sec = find(diff([gpsData(:).tET]) <= 80);
    gpsData(smallerThan80sec(1:2:end)) = [];
    
    gps_states = [gpsData(:).xx_J2000];
    
    % GPS measurements: [Epoch, Obj nr, position]
    gpsMeasObject = [   [gpsData(:).tET];   % et epoch
                        i*ones(1,length(gpsData)); % obj no.
                        gps_states(1:3,:)]; % position
    gpsMeas = [gpsMeas, gpsMeasObject];
end

% Sort measurements on epoch
[~, order] = sort(gpsMeas(1,:));
gpsMeas = gpsMeas(:,order);

end


function [fighandle] = plotOrbitAverage(latitudes,time,varargin)

[~,orbitIndices] = findpeaks(-latitudes);
nofFullOrbits = length(orbitIndices)-1;

validOrbits = ~isoutlier(diff(orbitIndices),'mean');

time_mean = zeros(1,nofFullOrbits);
for i=1:nofFullOrbits
    time_mean(i) = mean(time(orbitIndices(i):orbitIndices(i+1)-1));
end
time_mean = time_mean(validOrbits);

fighandle = figure;
for k = 1:nargin-2
    data = varargin{k};
    data_mean = zeros(1,nofFullOrbits);
    for i=1:nofFullOrbits
        data_mean(i)  = mean(data(orbitIndices(i):orbitIndices(i+1)-1));
    end
    data_mean = data_mean(validOrbits);
    plot(time_mean,data_mean); hold on;
end

end


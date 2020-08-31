function [fighandle,orbitAvgErr,time_mean] = plotOrbitAveragedDensityError(latitudes,time,trueDens,varargin)

[~,orbitIndices] = findpeaks(-latitudes);
nofFullOrbits = length(orbitIndices)-1;

validOrbits = ~isoutlier(diff(orbitIndices),'mean');

time_mean = zeros(1,nofFullOrbits);
for i=1:nofFullOrbits
    time_mean(i) = mean(time(orbitIndices(i):orbitIndices(i+1)-1));
end
time_mean = time_mean(validOrbits);

trueDens_mean = zeros(1,nofFullOrbits);
for i=1:nofFullOrbits
    trueDens_mean(i) = mean(trueDens(orbitIndices(i):orbitIndices(i+1)-1));
end
trueDens_mean = trueDens_mean(validOrbits);

fighandle = figure;
for k = 1:nargin-3
    dens = varargin{k};
    dens_mean = zeros(1,nofFullOrbits);
    for i=1:nofFullOrbits
        dens_mean(i)  = mean(dens(orbitIndices(i):orbitIndices(i+1)-1));
    end
    dens_mean = dens_mean(validOrbits);
    plot(time_mean,abs(dens_mean-trueDens_mean)./trueDens_mean*100); hold on;
    orbitAvgErr(k,:) = ((dens_mean-trueDens_mean)./trueDens_mean*100);
end

end


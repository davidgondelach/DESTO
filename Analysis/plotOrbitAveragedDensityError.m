function [fighandle,orbitAvgErr,time_mean] = plotOrbitAveragedDensityError(latitudes,time,trueDens,varargin)

[~,orbitIndices] = findpeaks(-latitudes);
nofFullOrbits = length(orbitIndices)-1;

time_mean = zeros(1,nofFullOrbits);
for i=1:nofFullOrbits
    time_mean(i) = mean(time(orbitIndices(i):orbitIndices(i+1)-1));
end

trueDens_mean = zeros(1,nofFullOrbits);
for i=1:nofFullOrbits
    trueDens_mean(i) = mean(trueDens(orbitIndices(i):orbitIndices(i+1)-1));
end

fighandle = figure;
for k = 1:nargin-3
    dens = varargin{k};
    dens_mean = zeros(1,nofFullOrbits);
    for i=1:nofFullOrbits
        dens_mean(i)  = mean(dens(orbitIndices(i):orbitIndices(i+1)-1));
    end
    plot(time_mean,abs(dens_mean-trueDens_mean)./trueDens_mean*100); hold on;
    orbitAvgErr(k,:) = ((dens_mean-trueDens_mean)./trueDens_mean*100);
end

end


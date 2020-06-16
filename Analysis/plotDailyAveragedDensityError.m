function [fighandle,dailyAvgerr] = plotDailyAveragedDensityError(time,trueDens,varargin)

timeDoy = floor(time);
doyStart = timeDoy(1);
doyEnd = timeDoy(end);
nofDoy = doyEnd-doyStart+1;

i=0;
trueDens_mean = zeros(1,nofDoy);
for doy=doyStart:doyEnd
    i=i+1;
    trueDens_mean(i) = mean(trueDens(timeDoy==doy));
end

fighandle = figure;
cmap=colormap('lines');
for k = 1:nargin-2
    dens = varargin{k};
    dens_mean = zeros(1,nofDoy);
    i=0;
    for doy=doyStart:doyEnd
        i=i+1;
        dens_mean(i)  = mean(dens(timeDoy==doy));
    end
    plot(doyStart:doyEnd,abs(dens_mean-trueDens_mean)./trueDens_mean*100,'Color',cmap(k,:)); hold on;
    dailyAvgerr(k,:) = (dens_mean-trueDens_mean)./trueDens_mean*100;
%     sum((dens_mean-trueDens_mean)./trueDens_mean*100 < 10)
%     sum((dens_mean-trueDens_mean)./trueDens_mean*100 < 5)
end

end


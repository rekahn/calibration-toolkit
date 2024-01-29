%% Filter out centered pings that are too low...
% use max mean TS

function [centerpings] = filtercenter(data)

fid1 = find(data.f_data >= data.FreqStart,1);
fid2 = find(data.f_data >= data.FreqEnd,1);

f = data.f_data(fid1:fid2);
TSdata = data.TS_data(fid1:fid2,:);

% query model to frequency vector
TSmodel = interp1(data.f_model,data.TS_model,f);

% calculate MSE between pings and theory
meanTS = 20*log10(mean(10.^(TSdata./20),1));

% take the best 15% of pings (minimum 10% of rmse)
k = round(length(meanTS)*0.15);
[~,pingids] = maxk(meanTS,k);

% rewrite usable pings
centerpings.TS_data = data.TS_data(:,pingids);
centerpings.TS_avg = 20*log10(mean(10.^(centerpings.TS_data./20),2));
centerpings.f_data = data.f_data;
centerpings.TS_model = data.TS_model;
centerpings.f_model = data.f_model;


figure; hold on; grid on;
plot(centerpings.f_model/1000,centerpings.TS_model,'r','LineWidth',2)
plot(data.f_data/1000,data.TS_data,'b')
plot(centerpings.f_data/1000,centerpings.TS_data,'g')
plot(centerpings.f_data/1000,centerpings.TS_avg,'k','LineWidth',2)
xlim([data.FreqStart data.FreqEnd]/1000)

end




























































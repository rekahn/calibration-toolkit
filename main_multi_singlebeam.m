% main script to run target localization, beam pattern calculation, and calibration
% SINGLE BEAM TRANSDUCER
% READS MULTIPLE FILES AT ONCE

% June 2020


% To change channels you need to change
    % - chan
    % - strings in data loading loop (lines 27-28)
    % - aT
    % - calibration sphere parameters
    % - range
    % - might have to change threshold in pingfilt function 
            % (right now this line is commented out so n/a)

%% LOAD DATA AND INITIALIZE

% cal.chan = 2; % which channel?

[fname,fpath] = uigetfile('MultiSelect','on');

for f = 1:length(fname)

if contains(fname{f},'ES') || contains(fname{f},'Airmar')
   continue
else
   load([fpath fname{f}])
   if exist([fpath fname{f}(1:end-4) '_ES38.mat'],'file')
       load([fpath fname{f}(1:end-4) '_ES38.mat'])
   else
       continue
   end
end


for ducer = 1:5
    if contains(data.echodata(ducer,1).channelID,'ES38')
        cal.chan1 = ducer;
    else
        continue
    end
end
for ducer = 1:2
    if contains(data.config.transceivers(4).channels(ducer).ChannelID,'ES38')
        cal.chan = ducer;
    else
        continue
    end
end




if f == 1
    cal.CompressedVoltage = CompressedVoltage; 
    cal.Range = Range;
    cal.FreqStart = FreqStart;
    cal.FreqEnd = FreqEnd;
    cal.echodata = data.echodata(cal.chan1,:);
    cal.config = data.config.transceivers(4);
    cal.param = data.param(cal.chan1,:);
    cal.Z = data.parameters.Ztrd;
    cal.parameters = data.parameters;
    cal.filters = data.filters(cal.chan1,:);
else
    cal.CompressedVoltage = [cal.CompressedVoltage CompressedVoltage]; 
    cal.echodata = [cal.echodata data.echodata(cal.chan1,:)];
    cal.param = [cal.param data.param(cal.chan1,:)];
end
clear CompressedVoltage
end


%% define parameters
depth = 0; % depth of transducer in m
D = 21.2; % diameter of the target [mm]; 38.1 mm or 21.2 mm WC sphere
rmin = 7.5; rmax = rmin+0.5; % range encompassing the target
    % for Mobile Bay 2020 ES70 cal: 8-8.5 m
    % for 2019 Deep-See in situ: 21-22 m
S=31.5; T=19.7; pH=8; % salinity, temperature, pH
offset = 47; % Mobile Bay: 51; DeepSee 2019: 61

% put parameters into a structure
cal.fnom = data.config.transceivers(4).channels(cal.chan).transducer.Frequency;
cal.fc = (cal.FreqStart+cal.FreqEnd)/2;
cal.depth = depth;
cal.D = D;
cal.rmin = rmin; cal.rmax = rmax;
cal.S = S; cal.T = T; cal.pH = pH;
cal.offset = offset;


%% PLOT AN ECHOGRAM

echogram(cal, 1);
%ylim([rmax+1 rmin-1])
clim([-40 -10])

%% Calculate the sphere model

% Calibration sphere info:
    % Tungsten carbide (WC) sphere either 38.1 mm or 21.2 mm diameter
    % Aluminum (Al) sphere 200 mm diameter

% Tungsten carbide sphere
[cal.TS_model,cal.f_model] = WC_TS_inputparams(cal.D, cal.S, cal.T, cal.depth, 39, 0);

% Aluminum sphere
%[cal.TS_model,cal.f_model] = Al_TS_inputparams(cal.D, cal.S, cal.T, cal.depth, 39, 0);   


%% CALCULATE TS OF top 10% OF PINGS
% remove top %1 and bottom 90% of pings
cal = calcTS_singlebeam(cal, cal.CompressedVoltage, cal.fnom/1000, 1,1); 

%% CALCULATE CALIBRATION CURVE FROM TOP PINGS

cal = calcG_singlebeam(cal);
ylim([16 26])

figure; hold on; grid on;
plot(cal.f_model/1000, cal.TS_model, 'r');
plot(cal.f_data/1000,cal.TS_avg-cal.offset,'g','LineWidth',2);
xlim([FreqStart/1000 FreqEnd/1000]); 
ylabel('TS (dB)'); xlabel('Frequency (kHz)');
legend('Theory','Calculated from top pings','Location','southeast');


%% MAKE SMALLER DATA STRUCTURE FOR SAVING AND RELOADING QUICKLY

gain.G = cal.G;
gain.Gsmooth = cal.Gsmooth;
gain.freq = cal.freq;
gain.FreqStart = cal.FreqStart;
gain.FreqEnd = cal.FreqEnd;
gain.fnom = cal.fnom;
gain.depth = cal.depth;
gain.S = cal.S; gain.T = cal.T; gain.pH = cal.pH;
gain.offset = cal.offset;
gain.D = cal.D;
gain.TS = cal.TS_avg;
gain.TS_model = cal.TS_model;
gain.f_model = cal.f_model;
gain.rmin = cal.rmin; gain.rmax = cal.rmax;


%% Check that this curve works



figure; hold on; grid on;
plot(cal.f_model/1000, cal.TS_model, 'r');
plot(cal.f_data/1000,cal.TS_avg-2*gain.Gsmooth,'g','LineWidth',2);
xlim([FreqStart/1000 FreqEnd/1000]); 
ylabel('TS (dB)'); xlabel('Frequency (kHz)');
legend('Theory','Calculated from top pings','Location','southeast');






% main script to run target localization, beam pattern calculation, and calibration

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

cal.chan = 3; % which channel?

[fname,fpath] = uigetfile('MultiSelect','on');

for f = 1:length(fname)

if contains(fname{f},'ES') || contains(fname{f},'Airmar')
   continue
else
   load([fpath fname{f}])
   if exist([fpath fname{f}(1:end-4) '_ES70.mat'],'file')
       load([fpath fname{f}(1:end-4) '_ES70.mat'])
   else
       continue
   end
end



if f == 1
    cal.CompressedVoltage = CompressedVoltage; 
    cal.Range = Range;
    cal.FreqStart = FreqStart;
    cal.FreqEnd = FreqEnd;
    cal.echodata = data.echodata(cal.chan,:);
    cal.config = data.config;
    cal.param = data.param(cal.chan,:);
    cal.Z = data.parameters.Ztrd;
    cal.parameters = data.parameters;
    cal.filters = data.filters;
else
    cal.CompressedVoltage = [cal.CompressedVoltage CompressedVoltage]; 
    cal.echodata = [cal.echodata data.echodata(cal.chan,:)];
    cal.param = [cal.param data.param(cal.chan,:)];
end
clear CompressedVoltage
end


%% define parameters
depth = 0; % depth of transducer in m
D = 38.1; % diameter of the target [mm]; 38.1 mm or 21.2 mm WC sphere
rmin = 7.7; rmax = rmin+0.5; % range encompassing the target
    % for Mobile Bay 2020 ES70 cal: 8-8.5 m
    % for 2019 Deep-See in situ: 21-22 m
S=32; T=19.7; pH=8; % salinity, temperature, pH
aT = 0.063/2; % transducer face radius - actual 0.1803/2
    % ES70: 0.181, ES120: 0.105, ES200: 0.063, ES333: 0.038
offset = 56; % Mobile Bay: 51; DeepSee 2019: 61

% put parameters into a structure
cal.fnom = str2num(data.config.transceivers(cal.chan).channels.transducer.Frequency);
cal.fc = (cal.FreqStart+cal.FreqEnd)/2;
cal.depth = depth;
cal.D = D;
cal.rmin = rmin; cal.rmax = rmax;
cal.S = S; cal.T = T; cal.pH = pH;
cal.aT = aT;
cal.offset = offset;


%% CHECK THE DATA TO SEE IF A QUADRANT IS DAMAGED

checkdata(cal, rmin, rmax);

%% PLOT AN ECHOGRAM

echogram(cal, 1);
%ylim([rmax+1 rmin-1])

%% Calculate the sphere model

% Calibration sphere info:
    % Tungsten carbide (WC) sphere either 38.1 mm or 21.2 mm diameter
    % Aluminum (Al) sphere 200 mm diameter

% Tungsten carbide sphere
[cal.TS_model,cal.f_model] = WC_TS_inputparams(cal.D, cal.S, cal.T, cal.depth, 39, 0);

% Aluminum sphere
%[cal.TS_model,cal.f_model] = Al_TS_inputparams(cal.D, cal.S, cal.T, cal.depth, 39, 0);   


%% FIND TARGET LOCATIONS AND OFF-AXIS ANGLES

cal = localize_multi(cal, 0);

% CALCULATE UNCOMPENSATED TS OF ALL THE DATA
%   and filter out bad pings
% take only pings with target within 3.5 degrees off-axis
cal = calcTS_multi(cal, cal.CompressedVoltage, 3.5, cal.fnom/1000, 0, 0);

%% PLOT BEAM PATTERN

crosshair(cal);
ylim([-2 2]); xlim([-2 2]);
%caxis([-56 -36]);
plotbeampattern(cal, offset);
%ylim([-56 -36]);


%% CALCULATE CALIBRATION CURVE FROM CENTERED PINGS

% calculate measured TS from centered pings (within 0.5 degrees off-axis)
cal.centerpings = calcTS(cal, cal.CompressedVoltage, 0.5, cal.fnom/1000, 1,1);

if isempty(cal.centerpings.TS_data)
    warning('No centered pings found');
else   
    % calculate calibration curve
    cal = calcG(cal, 'centerpings');
end

%% CALCULATE CALIBRATION CURVE FROM ALL PINGS, CORRECT FOR BEAM PATTERN

cal = beamcorrect(cal, 1);

if isempty(cal.centerpings.TS_data)
    cal = calcG(cal,'beamcorrected');
end


%% Compare TS calculated from just centered pings to beam-corrected pings
figure; hold on; grid on;
plot(cal.f_model/1000, cal.TS_model, 'r');
plot(cal.freq/1000,cal.beamcorrected.TS_avg-cal.offset,'g','LineWidth',2);
if ~isempty(cal.centerpings.TS_data)
    plot(cal.freq/1000,cal.centerpings.TS_avg-cal.offset,'k','LineWidth',2); 
end
xlim([FreqStart/1000 FreqEnd/1000]); 
ylabel('G (dB)'); xlabel('Frequency (kHz)');
legend('Theory','Beam-corrected','Calculated from centered pings','Location','southeast');


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
gain.aT = cal.aT; gain.D = cal.D;
gain.TS = cal.TS_avg;
gain.TS_model = cal.TS_model;
gain.f_model = cal.f_model;
gain.rmin = cal.rmin; gain.rmax = cal.rmax;


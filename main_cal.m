% main script to run target localization, beam pattern calculation, and calibration
% written specifically for EK80-ES70 transducer, but might work for others

% READS A SINGLE FILE

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

[fname,fpath] = uigetfile;
load([fpath fname])
load([fpath fname(1:end-4) '_ES200.mat'])

%%

cal.chan = 3; % CHECK THE CHANNEL NUMBER!
% define parameters
D = 21.2; % radius of the target [mm]; 38.1 mm or 21.2 mm
rmin = 7.7; rmax = rmin+0.5; % range encompassing the target
S=31.5; T=19.7; pH=8; % salinity, temperature, pH
aT = 0.063/2; % transducer face radius - actual 0.181/2
    % ES70: 0.181, ES120: 0.105, ES200: 0.063, ES333: 0.038
depth = 0;
offset = 57.5;
cal.fnom = 200000;


% put parameters into a structure
cal.D = D;
cal.depth = depth;
cal.rmin = rmin; cal.rmax = rmax;
cal.S = S; cal.T = T; cal.pH = pH;
cal.aT = aT;
cal.CompressedVoltage = CompressedVoltage;
cal.Range = Range;
cal.FreqStart = FreqStart;
cal.FreqEnd = FreqEnd;
%cal.fnom = str2num(data.config.transceivers(cal.chan).channels.transducer.Frequency);
cal.offset = offset;
cal.echodata = data.echodata;
cal.config = data.config;
cal.param = data.param;
cal.Z = data.parameters.Ztrd;
cal.parameters = data.parameters;
cal.filters = data.filters;

%% CHECK THE DATA TO MAKE SURE IT'S NOT FAULTY

checkdata(data, rmin, rmax);

%% PLOT AN ECHOGRAM

echogram(cal, 1);


%% Calculate the sphere model

% Calibration sphere info:
    % Tungsten carbide (WC) sphere either 38.1 mm or 21.2 mm diameter
    % Aluminum (Al) sphere 200 mm diameter

% Tungsten carbide sphere
[cal.TS_model,cal.f_model] = WC_TS_inputparams(cal.D, cal.S, cal.T, cal.depth, 39, 0);

% Aluminum sphere
%[cal.TS_model,cal.f_model] = Al_TS_inputparams(cal.D, cal.S, cal.T, cal.depth, 39, 0);   



%% FIND TARGET LOCATIONS AND OFF-AXIS ANGLES

cal = localize(cal, 1);

% CALCULATE TS OF ALL THE DATA
% and filter out any bad pings

cal = calcTS(cal, cal.CompressedVoltage, 3.5, cal.fnom/1000, 1,1); % use <=3.5

% plot crosshair colormap
crosshair(cal);

% PLOT BEAM PATTERN

plotbeampattern(cal, offset);

%% CALCULATE CALIBRATION CURVE FROM CENTERED PINGS

% calculate measured TS from centered pings (within 0.5 degrees off-axis)
cal.centerpings = calcTS(cal, cal.CompressedVoltage, 0.5, cal.fnom/1000, 1, 1);

if isempty(cal.centerpings.TS_data)
    warning('No centered pings found');
else   
    % calculate calibration curve
    cal = calcG(cal, 'centerpings');
end

%% CALCULATE CALIBRATION CURVE FROM ALL PINGS, CORRECT FOR BEAM PATTERN

cal = beamcorrect(cal, 0);

%% Compare TS calculated from just centered pings to beam-corrected pings
figure; hold on; grid on;
plot(cal.f_model/1000, cal.TS_model, 'r');
plot(cal.freq/1000,cal.beamcorrected.TS_avg-cal.offset,'g','LineWidth',2);
plot(cal.freq/1000,cal.centerpings.TS_avg-cal.offset,'k','LineWidth',2); 
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



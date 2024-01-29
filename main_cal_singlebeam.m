% main script to run target localization, beam pattern calculation, and calibration
% written specifically for EK80-ES70 transducer, but might work for others
% SINGLE BEAM TRANSDUCER
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
load([fpath fname(1:end-4) '_ES120.mat'])

%%

cal.chan = 1; % CHECK THE CHANNEL NUMBER!
% define parameters
D = 21.2; % radius of the target [mm]; 38.1 mm or 21.2 mm
rmin = 8; rmax = rmin+0.3; % range encompassing the target
S=31.5; T=19.7; pH=8; % salinity, temperature, pH
depth = 0;
offset = 63;
cal.fnom = 333000;


% put parameters into a structure
cal.D = D;
cal.depth = depth;
cal.rmin = rmin; cal.rmax = rmax;
cal.S = S; cal.T = T; cal.pH = pH;
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


%% PLOT AN ECHOGRAM

echogram(cal, 5);


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



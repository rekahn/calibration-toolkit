%% Calibration with automatic target detection
    % works for detecting just one target within given range

% create options structure
options.rmin = 9.5; % set range window to look for targets in
options.rmax = 10.5;
options.pldl = 6; % pulse length determination level (dB)
options.pl_max = 1; % min and max pulse length (m)
options.pl_min = 0.01;
options.max_stdalong = 0.1; % max standard deviation in angles
options.max_stdathw = 0.1;
options.TSthresh = -70; % minimum TS of a target (dB)
options.max_phi = 0.5; % max off-axis angle (degrees)
options.maxd = 0.5; % max distance between pulses (m)
% size of window for FFT calculation (m)
options.rabove = 0.1;
options.rbelow = 0.15;


%% Load in files


autocal.chan = 1; % which channel?

[fname,fpath] = uigetfile('MultiSelect','on'); 

for f = 1:length(fname)

% read in just the full data files
if contains(fname{f},'ES') || contains(fname{f},'Airmar')
   continue
else
    load([fpath fname{f}])
end

if f == 1
    autocal.echodata = data.echodata(autocal.chan,:);
    autocal.config = data.config;
    autocal.param = data.param(autocal.chan,:);
    autocal.Z = data.parameters.Ztrd;
    autocal.parameters = data.parameters;
else
    autocal.echodata = [autocal.echodata data.echodata(autocal.chan,:)];
    autocal.param = [autocal.param data.param(autocal.chan,:)];
end
end

% define parameters
depth = 0; % depth of transducer in m
D = 21.2; % diameter of the target [mm]; 38.1 mm or 21.2 mm WC sphere

S=34; T=20; pH=8; % salinity, temperature, pH
aT = 0.063/2; % transducer face radius - actual 0.1803/2
    % ES70: 0.181, ES120: 0.105, ES200: 0.063, ES333: 0.038
offset = 56; % Mobile Bay: 51; DeepSee 2019: 61

% put parameters into a structure
options.fnom = autocal.config.transceivers(autocal.chan).channels.transducer.Frequency;
options.FreqStart = autocal.param(1).FrequencyStart;
options.FreqEnd = autocal.param(1).FrequencyEnd;
options.fc = (options.FreqStart + options.FreqEnd)/2;
options.depth = depth;
options.D = D;
options.S = S; 
options.T = T; 
options.pH = pH;
options.aT = aT;
options.offset = offset;


% save options into autocal structure
autocal.options = options;


%% CHECK THE DATA TO SEE IF A QUADRANT IS DAMAGED

checkdata(autocal, options.rmin, options.rmax);


%% Detect targets in all pings
npings = size(autocal.echodata,2);
for n = 1:npings
    n
    autocal.targetdata(n).singletarget = findsingletargets(autocal, autocal.chan, n, options, 0);
end
fprintf('Done finding targets. \n')


%% Calculate calibration curve

% Which pings have on-axis targets?
targetpings = [];
for n = 1:length(autocal.targetdata)
    if ~isempty(autocal.targetdata(n).singletarget.targets)
        targetpings = [targetpings; n];
    end    
end


% Calculate theory curve
% Tungsten carbide sphere
[autocal.TS_model,autocal.f_model] = WC_TS_inputparams(options.D, options.S, options.T, options.depth, 39, 0);

% Aluminum sphere
%[cal.TS_model,cal.f_model] = Al_TS_inputparams(cal.D, cal.S, cal.T, cal.depth, 39, 0);

rawpings = [];
for n = 1:length(targetpings)
    ploc = targetpings(n);
    rawpings = [rawpings autocal.targetdata(ploc).singletarget.ping];
end
pingavg = mean(rawpings,2);


autocal = autoG(autocal, options, 1);
figure;hold on;
plot(autocal.freq./1000,autocal.G);
plot(autocal.freq./1000,autocal.Gsmooth,'k','LineWidth',2);
ylabel('G (dB)')
xlabel('freq (kHz)')




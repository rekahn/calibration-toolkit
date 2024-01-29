%% Trying to figure out how to make my cal match the EK80 cal tool


%% load the EK80 calibration curve

dir = 'F:\jellyfish\ek80\calibration\ES70_38mm\depth_1.5m\';
load([dir 'f_ek80.mat'])
load([dir 'G_ek80.mat'])

% load old gain structure and interpolate to match EK80
load ([dir 'gain70_1.5.mat'])
myoldG = interp1(gain70_15.freq,gain70_15.Gsmooth,f_ek80);

%plot
figure; hold on; grid on;
plot(f_ek80/1000,G_ek80)
plot(f_ek80/1000,myoldG)
legend('EK80 cal','my old cal')


%% redo my calibration calculation using just centered pings

load([dir 'cal_ES70_38mm-D20210817-T154538.mat'])
load([dir 'cal_ES70_38mm-D20210817-T154538_ES70.mat'])



%% create structure of parameters

cal.chan = 1; % CHECK THE CHANNEL
% define parameters
D = 38.1; % radius of the target [mm]; 38.1 mm or 21.2 mm
rmin = 1.6; rmax = rmin+0.5; % range encompassing the target
S=31.4; T=24.8; pH=8; % salinity, temperature, pH
aT = 0.181/2; % transducer face radius - actual 0.181/2
    % ES70: 0.181, ES120: 0.105, ES200: 0.063, ES333: 0.038
depth = 0;
offset = 24;
cal.fnom = 70000;
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
cal.offset = offset;


cal.echodata = data.echodata;
cal.config = data.config;
cal.param = data.param;


cal.Z = data.parameters.Ztrd;
cal.parameters = data.parameters;
cal.filters = data.filters;

% PLOT AN ECHOGRAM
echogram(cal, 1);

% Tungsten carbide sphere
[cal.TS_model,cal.f_model] = WC_TS_inputparams(cal.D, cal.S, cal.T, cal.depth, 39, 1);
xlim([40 95])



%% FIND TARGET LOCATIONS AND OFF-AXIS ANGLES

cal = localize(cal, 1);
cal = calcTS(cal, cal.CompressedVoltage, 3.5, cal.fnom/1000, 1,1);
plotbeampattern(cal,0);



%% CALCULATE CALIBRATION CURVE FROM CENTERED PINGS

% calculate measured TS from centered pings (within 0.5 degrees off-axis)
cal.centerpings = calcTS(cal, cal.CompressedVoltage, 0.5, cal.fnom/1000, 1, 1);

if isempty(cal.centerpings.TS_data)
    warning('No centered pings found');
else   
    % calculate calibration curve
    cal = calcG(cal, 'centerpings');
    scatter(f_ek80/1000,G_ek80); grid on;
end




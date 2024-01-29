%% Calculate volume scattering spectrum from calibration data
% to check against TS 
% run the main calibration script first

%function [Sv,f] = Sv_check(cal,Nsectors)
Nsectors = 4;

subrange = cal.Range(find(cal.Range>=cal.rmin,1):find(cal.Range>=cal.rmax,1));
CompressedVoltage = cal.CompressedVoltage(find(cal.Range>=cal.rmin,1):find(cal.Range>=cal.rmax,1),:);
S = cal.S;
T = cal.T;
ducerdepth = 0; %calibration


%CompressedVoltage = analysis.roi_data(chan).compressed;
range = subrange; %MFdata(chan).Range(rangeid1:rangeid2);
depth = range + ducerdepth;
fstart = cal.FreqStart;
fend = cal.FreqEnd;
fc = 0.5*(fstart+fend);
pH = cal.pH;
%r = range(end)-range(1);

% calculate the FFT and frequency vector
fsdec = 1/cal.param(1).SampleInterval;

% % calculate autocorrelated transmit signal spectrum
ytx = EK80calcSentSignal_1channel(cal, 1); 
ytxx = xcorr(ytx)/norm(ytx)^2;

% range compensation pre-FFT
CompressedVoltage = CompressedVoltage;

nfft = 2^13;
win = tukeywin(size(CompressedVoltage,2),0.05)' .* ones(size(CompressedVoltage));
CompressedVoltagewin = CompressedVoltage .* win;
CompressedVoltageFFT = (fft(CompressedVoltagewin,nfft));
ytxfft = fft(ytxx,nfft);
% create frequency vector for plotting
f_data = fsdec .* linspace(0,1-1/nfft,nfft);% + fsdec; % basebanded within fc 
% unwrap the signal and shift frequency vector to proper location
    % could take this into account by ffthsift and fsdec/2 to fvec in previous line
    if fc > fsdec/2
        CompressedVoltageFFT = [CompressedVoltageFFT; CompressedVoltageFFT; CompressedVoltageFFT];
        ytxfft = [ytxfft; ytxfft; ytxfft];
        f_data = [f_data fsdec+f_data 2*fsdec+f_data];
        idxmin = round((fc-fsdec/2)/fsdec*nfft);
        CompressedVoltageFFT = CompressedVoltageFFT(idxmin:idxmin+nfft-1,:);
        ytxfft = ytxfft(idxmin:idxmin+nfft-1,:);
        f_data = f_data(idxmin:idxmin+nfft-1);
    end
f = f_data';

% normalize by autocorrelated transmit signal spectrum
CompressedVoltageFFT = CompressedVoltageFFT./ytxfft;

% calculate the correction terms
alpha = alpha_sea(mean(depth), S, T, pH, f/1000);

%G = interp1(MFdata(chan).gain.freq, MFdata(chan).gain.Gsmooth, f);






% rimin = find(analysis.rmin <= volscat(chan).Range,1);
% rimax = find(analysis.rmax <= volscat(chan).Range,1);
zer = str2double(cal.config.Impedance);
zet = cal.Z;
   % Nu = [1 4 1 4 4];
   Nu = Nsectors;
Per = Nu*(abs(CompressedVoltageFFT)/(2*sqrt(2))).^2 *((zer+zet)/zer)^2.*1/zet;
%cw=sw_svel(analysis.S,analysis.T,r+analysis.dsdepth);
cw = sw_svel(S,T,ducerdepth);
Pet = cal.param(1).TransmitPower;
%psinom = cal.config.transducer.EquivalentBeamAngle;
psinom = str2double(cal.config.channels.transducer.EquivalentBeamAngle);
%psi = psinom + 20*log10(cal.config.transducer.Frequency./f); %cal.config.channels.transducer.Frequency
psi = psinom + 20*log10(str2double(cal.config.channels.transducer.Frequency)./f);

% take the top pings at center frequency
fcix = find(f >= fc,1);
[~,maxpings] = maxk(Per(fcix,:),round(size(Per,2)/10));

% Sv_all = 10*log10(Per) + 2*alpha*r - 10*log10(Pet*cw^3./2./(4*pi*f).^2)...
%     - 2*G - 10*log10((rangeid2-rangeid1)/fsdec) - psi + 20*log10(r);
% Sv_avg = 10*log10(mean(Per,2,'omitnan')) + 2*alpha*r - 10*log10(Pet*cw^3./2./(4*pi*f).^2)...
%     - 2*G - 10*log10((rangeid2-rangeid1)/fsdec) - psi + 20*log10(r);

% Sv_all = 10*log10(Per) + 2*alpha*r - 10*log10(Pet*cw^3./2./(4*pi*f).^2)...
%     - 2*G - 10*log10(length(range)/fsdec) - psi ;
% Sv_avg = 10*log10(mean(Per,2,'omitnan')) + 2*alpha*r - 10*log10(Pet*cw^3./2./(4*pi*f).^2)...
%     - 2*G - 10*log10(length(range)/fsdec) - psi ;

% uncalibrated Sv
Sv_all = 10*log10(Per(:,maxpings)) + 2*alpha*mean(range) - 10*log10(Pet*cw^3./2./(4*pi*f).^2)...
    - 10*log10(length(range)/fsdec) - psi + 20*log10(mean(range));
Sv_avg = 10*log10(mean(Per(:,maxpings),2,'omitnan')) + 2*alpha*mean(range) - 10*log10(Pet*cw^3./2./(4*pi*f).^2)...
    - 10*log10(length(range)/fsdec) - psi + 20*log10(mean(range)) ;

Sv_avg_novolume = Sv_avg + 10*log10(cw*length(range)/fsdec/2) + psi + 20*log10(mean(range));% -cal.offset;
TS_avg = 10*log10(mean(Per(:,maxpings),2,'omitnan')) + 40*log10(mean(range)) + 2*alpha*mean(range)...
    - 10*log10(Pet*cw^2./((4*pi*f).^2));% - cal.offset;

%Sv.Sv_all = Sv_all;
%Sv.Sv_avg = Sv_avg;



%% truncate Sv to be just the frequency range of the transducer

% fid1 = find(f >= gainstruct.usef1,1);
% fid2 = find(f >= gainstruct.usef2, 1);
% %fid1 = 1; fid2 = length(gainstruct.Gsmooth);
% 
% Sv(chan).freq = f_data(fid1:fid2)';
% Sv(chan).Sv_avg = analysis.roi_data(chan).Sv_avg(fid1:fid2);
% Sv(chan).Sv_all = analysis.roi_data(chan).Sv_all(fid1:fid2,:);



%end

%% more plotting w/ gain

figure;
subplot(2,1,2); hold on; grid on;
plot(gain.f_model/1000,gain.TS_model,'k')
plot(f/1000,TS_avg - 2*gain.Gsmooth,'LineWidth',2)
plot(f/1000,Sv_avg_novolume - 2*gain.Gsmooth,'LineWidth',2,'LineStyle','--')
legend('Theory','TS','S_v w/out volume')
xlim([gain.FreqStart gain.FreqEnd]/1000)
ylabel('dB'); xlabel('Frequency (kHz)')

subplot(2,1,1)
hold on; grid on;
plot(f/1000,TS_avg - 2*gain.Gsmooth,'LineWidth',2)
plot(f/1000,Sv_avg - 2*gain.Gsmooth,'LineWidth',2)
xlim([gain.FreqStart gain.FreqEnd]/1000)
ylabel('dB'); 
legend('TS','S_v')


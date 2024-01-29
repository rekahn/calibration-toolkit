% Calculate the TS of calibration target and plot with theoretical model
%% CORRECTED 4/2023 TO MATCH EK80 TOOL
% FOR A 4-SECTOR TRANSDUCER

function cal = calcTS(cal, compressed, ang_thresh, filter_freq, plot_flag,thresh_flag)  
% cal - data structure containing pertinent data, parameters, off-axis angles, etc.
% compressed - matrix of pulse-compressed data values
% ang_thresh - max off-axis angle of pings to keep
% filter_freq - frequency (kHz) to use to filter out bad pings
    % set to 0 to keep all pings     
% plot_flag - set 0 or 1 to plot the TS

rmin = cal.rmin; rmax = cal.rmax;
S = cal.S; T = cal.T; pH = cal.pH;
offset = cal.offset;


cal.phis = cal.localize.phis_polar;
cal.phis_along = cal.localize.phis_along;
cal.phis_athw = cal.localize.phis_athw;


% D =  diameter of the target [mm]; 38.1 mm or 21.2 mm
% rmin, rmax = range of depth containing the target
% offset = arbitrary - adjust as necessary so data match theory
% S, T, pH: salinity, temperature, pH


% take only data in range of target
r1 = rmin; r2 = rmax;
r = (r1+r2)/2;
Range = cal.echodata(cal.chan,1).range;
ind1 = find(Range > r1, 1);
ind2 = find(Range > r2, 1);
targetpings = compressed(ind1:ind2,:);

fsdec = 1/cal.param(cal.chan,1).SampleInterval;

FreqStart = cal.FreqStart;
FreqEnd = cal.FreqEnd;

% % calculate autocorrelated transmit signal spectrum
ytx = EK80calcSentSignal(cal, cal.chan, 1);
ytxx = xcorr(ytx)/norm(ytx)^2;


nfft = 2^13;
win = tukeywin(size(targetpings,2),0.05)' .* ones(size(targetpings));
targetpingswin = targetpings .* win;
pingsFFT = (fft(targetpingswin,nfft));
ytxfft = fft(ytxx,nfft);

% create frequency vector for plotting
f_data = fsdec .* linspace(0,1-1/nfft,nfft);% + fsdec; % basebanded within fc 


% unwrap the signal and shift frequency vector to proper location
    % could take this into account by ffthsift and fsdec/2 to fvec in previous line
    fc = (FreqStart+FreqEnd)/2;
    if fc > fsdec/2
        pingsFFT = [pingsFFT; pingsFFT; pingsFFT];
        ytxfft = [ytxfft; ytxfft; ytxfft];
        f_data = [f_data fsdec+f_data 2*fsdec+f_data];
        idxmin = round((fc-fsdec/2)/fsdec*nfft);
        pingsFFT = pingsFFT(idxmin:idxmin+nfft-1,:);
        ytxfft = ytxfft(idxmin:idxmin+nfft-1,:);
        f_data = f_data(idxmin:idxmin+nfft-1);
    end
cal.f_data = f_data';
 
% normalize by autocorrelated transmit signal spectrum
pingsFFT = pingsFFT./ytxfft;


% option to filter out targets outside specified off-axis angle
if ang_thresh ~=0
    [pingsFFT, ping_inds] = angfilt(cal, pingsFFT, ang_thresh);
    cal.phis(ping_inds) = [];
    cal.phis_along(ping_inds) = [];
    cal.phis_athw(ping_inds) = [];
end

% option to filter out the bad pings
if filter_freq ~= 0
    [pingsFFT, rmpings] = pingfilt(cal,pingsFFT,filter_freq); 
    cal.phis(find(rmpings == 1)) = [];
    cal.phis_along(find(rmpings == 1)) = [];
    cal.phis_athw(find(rmpings == 1)) = [];
end  

if thresh_flag == 1
% remove top 1% and bottom 10% of pings
    fidx = find(f_data >= fc,1);
    pvec = 20*log10(abs(pingsFFT(fidx,:)));
    percent1 = ceil(0.01 * length(pvec));
    percent10 = ceil(0.1 * length(pvec));
    [~,sortidx] = sort(pvec);
    top1 = sortidx(end-percent1:end);
    bot10 = sortidx(1:percent10);
    pingsFFT(:,[top1 bot10]) = [];
    cal.phis([top1 bot10]) = [];
    cal.phis_along([top1 bot10]) = [];
    cal.phis_athw([top1 bot10]) = [];
end
    
    
% calculate theoretical/model TS for tungsten carbide
%[TS_model,f_model] = WC_TS_inputparams(D, S, T, rmin, 39, 0);
% target diameter, salinity, temperature, depth, latitude, plot_flag
TS_model = cal.TS_model;
f_model = cal.f_model;

% calculate attenuation coefficient
alpha = alpha_sea(r, S, T, pH, f_data/1000)';

% %correct for power
% Kt = 2*cal.Z *(1500/cal.fnom)^2 /(16*pi^2);
% Kt2 = 2*cal.Z .*(1500./f_data').^2 ./(16*pi^2); % in line with Demer formulation

% %calculate empirical TS
% TS_data = 10*log10(abs(pingsFFT).^2)+ 2*alpha*r +10*log10(r^4)...
%     - 10*log10(Kt2*cal.param(1).TransmitPower.*ones(size(pingsFFT))); 
% TS_avg = 10*log10(mean(abs(pingsFFT).^2,2)) + 2*alpha*r +10*log10(r^4)...
%    - 10*log10(Kt2*cal.param(1).TransmitPower);


% % DEMER FORMULATION
% normalize FFT of MF output (pingsFFT) by autocorrelated transmit signal
% spectrum
zer = str2double(cal.config.transceivers(cal.chan).Impedance); % should be 108000 for WBT, 5400 for shipbased
%zer = cal.parameters.Rwbtrx; % what is this??
zet = cal.Z;
Per = (abs(pingsFFT)/sqrt(2)).^2 *((zer+zet)/zer)^2.*1/zet;
cw=sw_svel(cal.S,cal.T,cal.depth);
Pet = cal.param(cal.chan).TransmitPower;
TS_data = 10*log10(Per) + 40*log10(r) + 2*alpha*r...
    - 10*log10(Pet*cw^2./((4*pi*f_data').^2));
TS_avg = 10*log10(mean(Per,2)) + 40*log10(r) + 2*alpha*r...
    - 10*log10(Pet*cw^2./((4*pi*f_data').^2));


if plot_flag == 1 && ~isempty(TS_data)
    figure; hold on;
    plot(f_model/1000,TS_model,'r');
    plot(f_data/1000,TS_data-offset,'b');    
    plot(f_data/1000, TS_avg-offset, 'k','LineWidth',2);
    ylabel('TS (dB)'); xlabel('Frequency (kHz)');
    legend('Theory',['Measured, offset = ' num2str(offset) ' dB'],'Location','south');
    xlim([FreqStart/1000 FreqEnd/1000]); 
    grid on 
end

% pack output into data structure 
cal.TS_avg = TS_avg;
cal.TS_data = TS_data;

end
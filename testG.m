% apply my gain curve to the data to see how well it aligns with the
% theoretical TS of the sphere

%  testG(cal, gain70, 0.5, 70)


function testG(cal, gain, ang_thresh, filter_freq)


rmin = cal.rmin; rmax = cal.rmax;
S = cal.S; T = cal.T; pH = cal.pH;
offset = cal.offset;
compressed = cal.CompressedVoltage;

cal.phis = cal.localize.phis_polar;
cal.phis_along = cal.localize.phis_along;
cal.phis_athw = cal.localize.phis_athw;

% take only data in desired range
r1 = rmin; r2 = rmax;
r = (r1+r2)/2;
Range = cal.echodata(1).range;
ind1 = find(Range > r1, 1);
ind2 = find(Range > r2, 1);
targetpings = compressed(ind1:ind2,:);

fsdec = 1/cal.param(1).SampleInterval;

FreqStart = cal.FreqStart;
FreqEnd = cal.FreqEnd;

nfft = 2^13;
pingsFFT = (fft(targetpings,nfft));

% create frequency vector for plotting
f_data = fsdec .* linspace(0,1-1/nfft,nfft);% + fsdec; % basebanded within fc 

% unwrap the signal and shift frequency vector to proper location
    fc = (FreqStart+FreqEnd)/2;
    if fc > fsdec/2
        pingsFFT = [pingsFFT; pingsFFT; pingsFFT];
        f_data = [f_data fsdec+f_data 2*fsdec+f_data];
        idxmin = round((fc-fsdec/2)/fsdec*nfft);
        pingsFFT = pingsFFT(idxmin:idxmin+nfft-1,:);
        f_data = f_data(idxmin:idxmin+nfft-1);
    end
cal.f_data = f_data';
    
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
    
    
    
% calculate theoretical/model TS for tungsten carbide
[TS_model,f_model] = WC_TS_inputparams(cal.D, cal.S, cal.T, r, 39, 0);
% target diameter, salinity, temperature, depth, latitude, plot_flag

% calculate attenuation coefficient
alpha = alpha_sea(r, S, T, pH, f_data/1000)';

% correct for power
% Kt = 2*cal.Z *(1500/cal.fnom)^2 /(16*pi^2);
% Kt2 = 2*cal.Z .*(1500./f_data').^2 ./(16*pi^2); % in line with Demer formulation

% %calculate empirical TS
% TS_data = 10*log10(abs(pingsFFT).^2)+ 2*alpha*r +10*log10(r^4)...
%     - 10*log10(Kt2*cal.param(1).TransmitPower.*ones(size(pingsFFT))); 
% TS_avg = 10*log10(mean(abs(pingsFFT).^2,2)) + 2*alpha*r +10*log10(r^4)...
%    - 10*log10(Kt2*cal.param(1).TransmitPower);


% % DEMER FORMULATION for calibrated TS
zer = str2double(cal.config.transceivers(1).Impedance);
zet = cal.Z;
Per = (abs(pingsFFT)/sqrt(2)).^2 *((zer+zet)/zer)/zet;
cw=sw_svel(16,12,0);
Pet = cal.param(1).TransmitPower;
% TS_data = 10*log10(Per) + 40*log10(r) + 2*alpha*r...
%     - 10*log10(Pet*cw^2./(4*pi*f_data').^2);
TS_avg = 10*log10(mean(Per,2)) + 40*log10(r) + 2*alpha*r...
    - 10*log10(Pet*cw^2./(4*pi*f_data').^2) - 2*gain.Gsmooth;


% calculate calibrated TS using EK80 calibration
% ek80gain = interp1(ek80cal(:,1),ek80cal(:,2),f_data);
% TS_avg_ek80cal = 10*log10(mean(Per,2)) + 40*log10(r) + 2*alpha*r...
%     - 10*log10(Pet*cw^2./(4*pi*f_data').^2) - 2*ek80gain';


    figure; hold on;
    plot(f_model/1000,TS_model,'k');
    plot(f_data/1000,TS_avg,'LineWidth',2);    
    %plot(f_data/1000, TS_avg_ek80cal,'LineWidth',2);
    ylabel('TS (dB)'); xlabel('Frequency (kHz)');
    legend('Theory',"Rachel's cal",'EK80 cal','Location','south');
    xlim([FreqStart/1000 FreqEnd/1000]); 
    grid on 



end
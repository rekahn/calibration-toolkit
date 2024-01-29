% Calculate the calibration curve

function autocal = autoG(autocal, options, plot_flag)

% Which pings have on-axis targets?
targetpings = [];
for n = 1:length(autocal.targetdata)
    if ~isempty(autocal.targetdata(n).singletarget.targets)
        targetpings = [targetpings; n];
    end    
end
autocal.targetpings = targetpings;


rawpings = []; rs = [];
for n = 1:length(targetpings)
    ploc = targetpings(n);
    rawpings = [rawpings autocal.targetdata(ploc).singletarget.ping];
    rs = [rs autocal.targetdata(ploc).singletarget.targets(1).r];
end
rtarget = mean(rs);


%% calculate the average TS
fsdec = 1/autocal.param(1).SampleInterval;
FreqStart = options.FreqStart;
FreqEnd = options.FreqEnd;
nfft = 2^13;
pingsFFT = (fft(rawpings,nfft));

% create frequency vector for plotting
f_data = fsdec .* linspace(0,1-1/nfft,nfft)';

% unwrap the signal and shift frequency vector to proper location
    % could take this into account by ffthsift and fsdec/2 to fvec in previous line
    fc = options.fc;
    if fc > fsdec/2
        pingsFFT = [pingsFFT; pingsFFT; pingsFFT];
        f_data = [f_data fsdec+f_data 2*fsdec+f_data];
        idxmin = round((fc-fsdec/2)/fsdec*nfft);
        pingsFFT = pingsFFT(idxmin:idxmin+nfft-1,:);
        f_data = f_data(idxmin:idxmin+nfft-1);
    end
autocal.f_data = f_data;

TS_model = autocal.TS_model;
f_model = autocal.f_model;

% calculate attenuation coefficient
alpha = alpha_sea(rtarget, options.S, options.T, options.pH, f_data/1000);

% % correct for power
% Kt = 2*autocal.Z *(1500/options.fnom)^2 /(16*pi^2);
% 
% %calculate empirical TS
% TS_data = 10*log10(abs(pingsFFT).^2)+ 2*alpha'*rtarget +10*log10(rtarget^4)...
%     - 10*log10(Kt*autocal.param(1).TransmitPower); 
% TS_avg = 10*log10((mean(abs(pingsFFT),2).^2)) + 2*alpha'*rtarget +10*log10(rtarget^4)...
%     - 10*log10(Kt*autocal.param(1).TransmitPower);

% DEMER FORMULATION
zer = autocal.parameters.Rwbtrx;
zet = autocal.Z;
Per = (abs(pingsFFT)/sqrt(2)).^2 *((zer+zet)/zer)/zet;
cw=sw_svel(34,20,0);
Pet = autocal.param(1).TransmitPower;
TS_data = 10*log10(Per) + 40*log10(rtarget) + 2*alpha'*rtarget...
    - 10*log10(Pet*cw^2./(4*pi*f_data').^2);
TS_avg = 10*log10(mean(Per,2)) + 40*log10(rtarget) + 2*alpha'*rtarget...
    - 10*log10(Pet*cw^2./(4*pi*f_data').^2);


% pack output into data structure 
autocal.TS_avg = TS_avg;
autocal.TS_data = TS_data;

%% calculate calibration curve

TS_theory = interp1(f_model,TS_model,f_data);
G = 0.5*(TS_avg-TS_theory');
Gsmooth = smoothG(f_data, G, options);

% pack data into structure
autocal.G = G;
autocal.Gsmooth = Gsmooth;
autocal.freq = f_data';

if plot_flag == 1 && ~isempty(TS_data)
    figure; hold on;
    plot(f_model/1000,TS_model,'r');
    plot(f_data/1000,TS_data-options.offset,'b');    
    plot(f_data/1000, TS_avg-options.offset, 'k','LineWidth',2);
    ylabel('TS (dB)'); xlabel('Frequency (kHz)');
    legend('Theory',['Measured, offset = ' num2str(options.offset) ' dB'],'Location','south');
    xlim([FreqStart/1000 FreqEnd/1000]); 
    grid on 
    
    
    figure;
    plot(f_data/1000,G,'blue'); hold on;
    plot(f_data/1000,Gsmooth,'black','LineWidth',1.5);
    xlim([FreqStart/1000 FreqEnd/1000]);
    ylabel('G (dB)'); xlabel('Frequency (kHz)');
    ylim([10 33]);   
end

end
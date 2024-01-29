% Calculate the beam pattern corrected TS

function cal = beamcorrect(cal, plot_flag)

TS_avg = cal.TS_avg;
TS_data = cal.TS_data;
f_data = cal.f_data;
TS_model = cal.TS_model;
f_model = cal.f_model;
phis_polar = cal.phis;
aT = cal.aT;
offset = cal.offset;



% calculate beam pattern
D_T = zeros(length(f_data),length(phis_polar));
k_data = 2*pi*f_data/1500;
for n=1:length(phis_polar)
    
    theta = phis_polar(n) * pi/180;
    D_T(:,n) = (2*besselj(1,k_data*aT*sin(theta))./(k_data*aT*sin(theta))).^2;
    
end

% get rid of NaN values
offax_correct = 10*log10(D_T.^2);
offax_correct(isnan(offax_correct)) = 0;

% sum to get beam-corrected TS
TS_beamcorrect = TS_data - offax_correct;


% Do it for the average TS
phiavg = mean(phis_polar);
theta_avg = phiavg * pi/180;
D_Tavg = (2*besselj(1,k_data*aT*sin(theta_avg))./(k_data*aT*sin(theta_avg))).^2;
TS_avgbeamcorr = TS_avg - 20*log10(D_Tavg);


FreqStart = cal.FreqStart;
FreqEnd = cal.FreqEnd;


if plot_flag == 1;
    figure; hold on;
    plot(f_model/1000,TS_model,'r');
    plot(f_data/1000,TS_beamcorrect-offset,'b'); 
    plot(f_data/1000,TS_avgbeamcorr-offset,'k','LineWidth',2);
    ylabel('TS (dB)'); xlabel('Frequency (kHz)');
    legend('Theory',['Beam-corrected data, offset = ' num2str(offset) ' dB'],'Location','south');
    xlim([FreqStart/1000 FreqEnd/1000]); 
    grid on 
end

% pack files into data structure
cal.beamcorrected.TS_avg = TS_avgbeamcorr;
cal.beamcorrected.TS = TS_beamcorrect;
cal.beamcorrected.freq = f_data;


end
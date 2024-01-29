%% calculate beam corrected TS for a given angle, single beam

function cal = singlebeam_beamcorrect(cal,offangle, aT)

D_T = zeros(length(cal.f_data),1);
k_data = 2*pi*cal.f_data/1500;
    
    theta = offangle * pi/180;
    D_T = (2*besselj(1,k_data*aT*sin(theta))./(k_data*aT*sin(theta))).^2;
  

% get rid of NaN values
offax_correct = 10*log10(D_T.^2);
offax_correct(isnan(offax_correct)) = 0;

% sum to get beam-corrected TS
TS_beamcorrect = cal.TS_avg + offax_correct;

cal.TS_beamcorrect = TS_beamcorrect;


end






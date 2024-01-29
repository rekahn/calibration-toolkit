% calculate the measured beam pattern and compare with theory

function plotbeampattern(cal, offset)

aT = cal.aT;
phis_polar = cal.phis;
TS_data = cal.TS_data;
f_data = cal.f_data;
TS_model = cal.TS_model;
f_model = cal.f_model;


% measured beam pattern
FreqStart = cal.FreqStart;
FreqEnd = cal.FreqEnd;
fc = (FreqStart+FreqEnd)/2;

ind_fc1 = find(f_data > fc, 1);
TS_beam = TS_data(ind_fc1-1, :);

%calculate theoretical beampattern
ind_fc2 = find(f_model > fc, 1);
TSc = TS_model(ind_fc2-1);

phivec = linspace(0,max(phis_polar),300);

k = 2*pi*fc/1500;
theta = phivec * pi/180;
D_T = (2*besselj(1,k*aT*sin(theta))./(k*aT*sin(theta))).^2; % intensity beam pattern 

% get rid of NaN values
offax_correct = 20*log10(D_T); % two-way beam pattern
offax_correct(isnan(offax_correct)) = 0;

% sum to get off-axis TS
TS_offaxis = TSc + offax_correct;


figure; hold on;
scatter(phis_polar, TS_beam - offset,'filled');
line(phivec, TS_offaxis, 'Color', 'r', 'LineWidth',2);
ylabel('TS (dB)'); xlabel('Off-axis angle (\circ)');
legend(['Measured, offset = ', num2str(offset), ' dB'],'Theory','Location','southwest');


end




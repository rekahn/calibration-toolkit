% calculate calibration curve

function cal = calcG_singlebeam(cal,beamcorrect)

    TS_avg = cal.TS_avg;
    f_data = cal.f_data;
    TS_model = cal.TS_model;
    f_model = cal.f_model;

    if beamcorrect ==1
        TS_avg = cal.TS_beamcorrect;
    end


TS_theory = interp1(f_model,TS_model,f_data);
G = 0.5*(TS_avg-TS_theory);
Gsmooth = smoothG(f_data, G, cal);

FreqStart = cal.FreqStart;
FreqEnd = cal.FreqEnd;

figure;
plot(f_data/1000,G,'blue'); hold on;
plot(f_data/1000,Gsmooth,'black','LineWidth',1.5);
xlim([FreqStart/1000 FreqEnd/1000]);
ylabel('G (dB)'); xlabel('Frequency (kHz)');
ylim([20 35]);

% pack data into structure
cal.G = G;
cal.Gsmooth = Gsmooth;
cal.freq = f_data;

end
% calculate calibration curve

function cal = calcG(cal, method)

switch method
    case 'centerpings'
    TS_avg = cal.centerpings.TS_avg;
    f_data = cal.centerpings.f_data;
    TS_model = cal.centerpings.TS_model;
    f_model = cal.centerpings.f_model;
    case 'beamcorrected'
    TS_avg = cal.beamcorrected.TS_avg;
    f_data = cal.beamcorrected.freq;
    TS_model = cal.TS_model;
    f_model = cal.f_model;
    otherwise
    error(['''' method '''' ' is not an option. Try ''centerpings'' or ''beamcorrected''.'])
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
ylim([10 33]);

% pack data into structure
cal.G = G;
cal.Gsmooth = Gsmooth;
cal.freq = f_data;

end
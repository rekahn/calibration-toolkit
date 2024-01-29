% plot an echogram, indicate the specified range


function echogram(datastruct, downsamp)
% downsamp = factor to downsample by to facilitate plotting

Range = datastruct.Range;
fc = (datastruct.FreqStart + datastruct.FreqEnd)/2;
pingnos = 1:size(datastruct.CompressedVoltage,2);

alpha = alpha_sea(datastruct.rmin, datastruct.S, datastruct.T, datastruct.pH, fc/1000);

% create range matrix
rangeVec = 10.*log10(datastruct.Range.^2);
rangeMat = ones(size(datastruct.CompressedVoltage)) .* rangeVec;

Vc_plot = 10*log10(abs(datastruct.CompressedVoltage).^2) + 10*log10(rangeMat.^2) + 2*alpha.*rangeMat;

figure; hold on;
if downsamp == 0
    pcolor(pingnos,datastruct.Range,Vc_plot);
else
   pcolor(pingnos,datastruct.Range(1:downsamp:end),Vc_plot(1:downsamp:end,:)); 
end
shading flat;
set(gca, 'YDir', 'reverse')
colormap('jet')
ylabel('Depth (m)'); xlabel('Ping number');
colorbar
caxis([-70 -10])
xlim([pingnos(1) pingnos(end)]); ylim([Range(end) Range(1)])

plot(pingnos,ones(size(pingnos))*datastruct.rmin, 'magenta', 'LineWidth',2);
plot(pingnos,ones(size(pingnos))*datastruct.rmax, 'magenta', 'LineWidth',2);
end
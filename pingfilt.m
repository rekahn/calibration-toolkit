% filters out bad pings based on TS at the nominal center frequency

function [pingsFFT, rmpings] = pingfilt(cal,pingsFFT,filter_freq)

% filter out bad pings by removing outliers at nominal center frequency
%fnom = cal.config.transceivers(1).channels.transducer.Frequency;
%[~,idx] = find(f_data > fnom, 1);

% filter out bad pings by removing outliers at specified frequency
% filter_freq
[idx,~] = find(cal.f_data > filter_freq*1000, 1);
FFT_filt = 20*log10(abs(pingsFFT(idx-1,:)));
%toss_ind = find(FFT_filt < -67);
%FFT_filt(toss_ind) = NaN;

[FFT_filt2, rmpings] = rmoutliers(FFT_filt);
%rmpings = zeros(size(FFT_filt));
%rmpings(toss_ind) = 1;

%phis_polar(find(rmpings == 1)) = [];
pingsFFT(:,find(rmpings~=0)) = [];

end
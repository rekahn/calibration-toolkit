% filter out pings outside of specified angle

function [pingsFFT, ping_inds] = angfilt(cal, pingsFFT, ang_thresh)

phis = cal.localize.phis_polar;

ping_inds = find(abs(phis) > ang_thresh);

pingsFFT(:,ping_inds) = [];

end



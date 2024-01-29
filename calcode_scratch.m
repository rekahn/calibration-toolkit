% calibration-related code
% will probably need to be organized into separate scripts...


% load some data
load('C:\Users\rkahn\Google Drive\acoustics_code\summer2020\Calibration70k\MB_70k_cal_38p1-D20200310-T191258.mat')
%load('C:\Users\rkahn\Google Drive\acoustics_code\summer2020\Calibration70k\MB_70k_cal_38p1-D20200310-T191258_ES70.mat')

%%

% Emma's code:
% yfore   = sum(yc(:,3:4),2)/2;
% yaft    = sum(yc(:,1:2),2)/2;
% ystar   = (yc(:,1) + yc(:,4))/2;
% yport   = sum(yc(:,2:3),2)/2;
%         
% phialong  = angle(yfore.*conj(yaft)) *180/pi / (sensalong * fc/fnom);
% phiathw  = angle(ystar.*conj(yport)) *180/pi / (sensathw  * fc/fnom);
%         
% phi = sqrt(phialong.^2 + phiathw.^2);


% create a 4xn matrix with the 4-channel compressed voltage for one ping
pingno = 100;
range = data.echodata(1,pingno).range;
pingdata = data.echodata(1,pingno).compressed;

% % sampling period
% T = data.param(1,pingno).SampleInterval;

% center frequencies
fc = (data.param(1,pingno).FrequencyStart + data.param(1,pingno).FrequencyEnd)/2;
fnom = 70e3;

% sensitivities?
sensalong = -190 ;
sensathw = -190 ;

% take only data in range of target
r1 = 7.8; r2 = r1+1; % m
ind1 = find(range > r1, 1);
ind2 = find(range > r2, 1);
range = range(ind1:ind2);
pingdata = pingdata(ind1:ind2,:);

yfore = (pingdata(:,1) + pingdata(:,2))/2;
yaft = (pingdata(:,3) + pingdata(:,4))/2;
yport = (pingdata(:,2) + pingdata(:,3))/2;
ystar = (pingdata(:,1) + pingdata(:,4))/2;

% raise an error if one of the channels is shot
if abs(mean(abs(yfore)) - mean(abs(yaft))) > 1.0 * mean(abs(yfore))  % CHANGE THIS FACTOR TO SOMETHING LIKE 0.15, 1.0 overrides the error
    error('Error: significant discrepancy between quadrants')
end

% calculate off-axis angles
 phialong  = angle(yfore.*conj(yaft)) *180/pi / (sensalong * fc/fnom);
 phiathw  = angle(ystar.*conj(yport)) *180/pi / (sensathw  * fc/fnom);
%         
phivec = sqrt(phialong.^2 + phiathw.^2);
phi = mean(phivec)




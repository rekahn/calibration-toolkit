% Find target locations

function cal = localize_multi(cal, plot_flag)

% unpack parameters
CompressedVoltage = cal.CompressedVoltage;
rmin = cal.rmin; rmax = cal.rmax;


npings = size(cal.echodata,2);
fnom = cal.fnom;
    if ischar(fnom)
        fnom = str2double(fnom);
    end


pings = 1:npings;
phis_along = 1:length(pings);
phis_athw = 1:length(pings);
phis_polar = 1:length(pings);

for n = 1:npings
    
    pingno = n;
    range = cal.echodata(pingno).range;
    pingdata = cal.echodata(pingno).compressed;

    if isempty(pingdata)
        continue
    end

    % center frequencies
    fc = (cal.param(pingno).FrequencyStart + cal.param(pingno).FrequencyEnd)/2;
    
    % sensitivities?
    %sensalong = cal.config.transceivers(cal.chan).channels.transducer.AngleSensitivityAlongship;
    %sensathw = cal.config.transceivers(cal.chan).channels.transducer.AngleSensitivityAthwartship;
    sensalong = cal.config.channels.transducer.AngleSensitivityAlongship;
    sensathw = cal.config.channels.transducer.AngleSensitivityAthwartship;

    if ischar(sensalong)
        sensalong = str2double(sensalong);
        sensathw = str2double(sensathw);
    end
    
    % take only data in range of target
    r1 = rmin; r2 = rmax;
    ind1 = find(range > r1, 1);
    ind2 = find(range > r2, 1);
    target_range = range(ind1:ind2);
    pingdata = pingdata(ind1:ind2,:);
    target_Vc = CompressedVoltage(ind1:ind2,:);

    yfore = (pingdata(:,1) + pingdata(:,2))/2;
    yaft = (pingdata(:,3) + pingdata(:,4))/2;
    yport = (pingdata(:,2) + pingdata(:,3))/2;
    ystar = (pingdata(:,1) + pingdata(:,4))/2;

    % calculate off-axis angles in degrees
    phialong  = angle(yfore.*conj(yaft)) *180/pi / (sensalong * fc/fnom);
    phiathw  = angle(ystar.*conj(yport)) *180/pi / (sensathw  * fc/fnom);         
    phivec = sqrt(phialong.^2 + phiathw.^2); 
    
    % take just the angle at the voltage peak
    yavg = (yfore + yaft + yport + ystar)/4;
    [~,peakind] = max(abs(yavg));
    phis_along(n) = phialong(peakind);
    phis_athw(n) = phiathw(peakind);
    phis_polar(n) = phivec(peakind);
       
end

% plot
if plot_flag == 1
    figure; hold on;
    yline(0,'Color','k'); xline(0,'Color','k');
    xcirc1 = 0.5*cos(0:0.01:2*pi); ycirc1 = 0.5*sin(0:0.01:2*pi);
    line(xcirc1,ycirc1,'Color','k','LineStyle','--');
    for a = 1:5
    xcirc2 = a*cos(0:0.01:2*pi); ycirc2 = a*sin(0:0.01:2*pi);
    line(xcirc2,ycirc2,'Color','k');
    end
    xlim([-5 5]); ylim([-5 5]);
    pbaspect([1 1 1]);
    
    plot1 = scatter(phis_athw,phis_along,'filled');
    plot1.MarkerFaceAlpha = 0.2;
    xlabel('Athwartship angle (\circ)');
    ylabel('Alongship angle (\circ)');

end

% pack output into data structure
cal.localize.range = target_range;
cal.localize.compressed = target_Vc;
cal.localize.phis_along = phis_along;
cal.localize.phis_athw = phis_athw;
cal.localize.phis_polar = phis_polar;

end

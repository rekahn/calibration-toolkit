% get a vector of along and athward angles from a single ping

function [phialong,phiathw,phivec] = angle_vector(channel, pingno, range, data)

fnom = data.config.transceivers(channel).channels.transducer.Frequency;

    % center frequencies
    fc = (data.param(channel).FrequencyStart + data.param(channel).FrequencyEnd)/2;
    
    % sensitivities?
    sensalong = data.config.transceivers(channel).channels.transducer.AngleSensitivityAlongship;
    sensathw = data.config.transceivers(channel).channels.transducer.AngleSensitivityAthwartship;

    % take only data in range of target
    ind1 = find(data.echodata(channel,pingno).range == range(1));
    ind2 = find(data.echodata(channel,pingno).range == range(end));
    pingdata = data.echodata(channel,pingno).compressed(ind1:ind2,:);

    yfore = (pingdata(:,1) + pingdata(:,2))/2;
    yaft = (pingdata(:,3) + pingdata(:,4))/2;
    yport = (pingdata(:,2) + pingdata(:,3))/2;
    ystar = (pingdata(:,1) + pingdata(:,4))/2;

    % calculate off-axis angles in degrees
    phialong  = angle(yfore.*conj(yaft)) *180/pi / (sensalong * fc/fnom);
    phiathw  = angle(ystar.*conj(yport)) *180/pi / (sensathw  * fc/fnom);         
    phivec = sqrt(phialong.^2 + phiathw.^2); 



end
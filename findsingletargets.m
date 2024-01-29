% function to find single targets in an individual ping and output
% a data structure containing the spectra of all detected targets

function singletarget = findsingletargets(datastruct, channel, pingno, options, plot_flag)

% Set plot_flag to 0 or 1 to display plots


singletarget = [];
% unpack options
rmin = options.rmin; % set range window to look for targets in
rmax = options.rmax;
pldl = options.pldl; % pulse length determination level (dB)
pl_max = options.pl_max; % min and max pulse length (m)
pl_min = options.pl_min;
max_stdalong = options.max_stdalong; % max standard deviation in angles
max_stdathw = options.max_stdathw;
TSthresh = options.TSthresh; % minimum TS of a target (dB)
max_phi = options.max_phi; % max off-axis angle (degrees)
maxd = options.maxd; % max distance between pulses (m)
rabove = options.rabove;
rbelow = options.rbelow;

% get compressed voltage for each quadrant and average of all quadrants
pingdata = datastruct.echodata(channel,pingno).compressed;
ping = mean(pingdata, 2); % mean compressed voltage
fullrange = datastruct.echodata(pingno).range;
FreqStart = datastruct.param(channel,pingno).FrequencyStart;
FreqEnd = datastruct.param(channel,pingno).FrequencyEnd;

% take just the data in the specified range window
rix1 = find(fullrange > rmin,1);
rix2 = find(fullrange > rmax,1);
singletarget.range = fullrange(rix1:rix2);
singletarget.ping = ping(rix1:rix2);
singletarget.pingdata = pingdata(rix1:rix2,:);
singletarget.P = 20*log10(abs(singletarget.ping)); % power

% create figure to display results
singletarget.pingplot = figure('visible','off'); hold on;
plot(fullrange,20*log10(abs(ping)),'Color',[0.5 0.5 0.5]);
plot(singletarget.range, singletarget.P);

% ----------------------------------------------------------------

% FIND PEAKS
    [pks,locs] = findpeaks(singletarget.P,'minpeakwidth',3);

    singletarget.pks = pks;
    singletarget.locs = locs;

    % structure the targets
    for n = 1:length(singletarget.pks)
       singletarget.targets(n).P = pks(n);
       singletarget.targets(n).loc = locs(n);
       singletarget.targets(n).r = singletarget.range(locs(n));
       singletarget.targets(n).TVG = 40*log10(singletarget.range(locs(n))) ...
           + 2*singletarget.range(locs(n))*alpha_sea(fullrange(locs(n)),24,20,8,70); % PUT THESE IN OPTIONS
       singletarget.targets(n).TS = singletarget.targets(n).P + singletarget.targets(n).TVG;
    end

    scatter(singletarget.range(locs),pks,'filled');

% Filter by length of pulse envelope
    rm = [];
    for n = 1:length(pks)
        plevel = pks(n) - pldl;
        loc = singletarget.targets(n).loc;

        idx1 = loc+1 - find(flip(singletarget.P(1:loc) <= plevel),1);
        idx2 = find(singletarget.P(loc:end) <= plevel , 1) + loc-1;

        if ~isempty(idx1) && ~isempty(idx2)
            singletarget.targets(n).pulselength = singletarget.range(idx2) - singletarget.range(idx1);
            singletarget.targets(n).envstart = idx1; 
            singletarget.targets(n).envend = idx2;
            singletarget.targets(n).P_envelope = singletarget.P(idx1:idx2);
            singletarget.targets(n).r_envelope = singletarget.range(idx1:idx2);
        else
            rm = [rm n];
        end
    end

    % keep only targets within the desired pulse range
    pl_idx = find([singletarget.targets.pulselength] <= pl_max & ...
        [singletarget.targets.pulselength] >= pl_min);
    pl_idx = setdiff(pl_idx, rm);

    singletarget.targets = singletarget.targets(pl_idx);
    scatter(singletarget.range([singletarget.targets.loc]),[singletarget.targets.P],'filled');

    if isempty(singletarget.targets)
        return
    end

% FILTER BY STANDARD DEVIATION OF ANGLE
    for n = 1:length(singletarget.targets)
        subrange = singletarget.targets(n).r_envelope;
        [phialong, phiathw] = angle_vector(channel, pingno, subrange, datastruct);
        singletarget.targets(n).phialong = phialong;
        singletarget.targets(n).phiathw = phiathw;
        singletarget.targets(n).stdev_along = std(phialong);
        singletarget.targets(n).stdev_athw = std(phiathw);
    end

singletarget.targets = singletarget.targets([singletarget.targets.stdev_along] <= max_stdalong ...
         & [singletarget.targets.stdev_athw] <= max_stdathw);
    scatter(singletarget.range([singletarget.targets.loc]),[singletarget.targets.P],'filled');

    if isempty(singletarget.targets)
        return
    end

% FILTER BY TS
    singletarget.targets = singletarget.targets([singletarget.targets.TS] >= TSthresh);
    scatter(singletarget.range([singletarget.targets.loc]),[singletarget.targets.P],'filled');

    if isempty(singletarget.targets)
        return
    end
    
% REJECT TARGETS OUTSIDE MAX OFF-AXIS ANGLE
    [~,~,phivec] = angle_vector(channel, pingno, singletarget.range, datastruct);
    for n = 1:length(singletarget.targets)
        singletarget.targets(n).phi = phivec(singletarget.targets(n).loc);
    end

    singletarget.targets = singletarget.targets([singletarget.targets.phi] <= max_phi);
    scatter(singletarget.range([singletarget.targets.loc]),[singletarget.targets.P],'filled');

    if isempty(singletarget.targets)
        return
    end
    
 
% REJECT PULSES THAT ARE TOO CLOSE TOGETHER
    % (keep the one with the higher TS)
    combos = combnk(1:length(singletarget.targets),2); % all combos of 2 targets
    rm = [];
    for c = 1:size(combos,1)   
        n1 = combos(c,1);
        n2 = combos(c,2);
        r1 = singletarget.targets(n1).r;
        r2 = singletarget.targets(n2).r;

        % find end of first pulse envelope and start of the second one
        if r1 < r2
            d1 = singletarget.range(singletarget.targets(n1).envend);
            d2 = singletarget.range(singletarget.targets(n2).envstart);
        else
            d1 = singletarget.range(singletarget.targets(n2).envend);
            d2 = singletarget.range(singletarget.targets(n1).envstart);
        end

        if d2-d1 <= maxd
            if singletarget.targets(n1).TS > singletarget.targets(n2).TS
                rm = [rm n2];
            else
                rm = [rm n1];
            end  
        end   
    end

    rm = unique(rm);
    singletarget.targets(rm) = [];
    
    if isempty(singletarget.targets)
        return
    end

% finish plotting
scatter(singletarget.range([singletarget.targets.loc]),[singletarget.targets.P],...
    140,'filled','p','MarkerFaceColor','y','MarkerEdgeColor','k');
title('Spatial (temporal) domain')
legend('Whole ping','Search window','All peaks','Filtered by pulse length','Filtered by StDev angle',...
    'Filtered by TS','Only on-axis','Detected targets'...
    ,'Location','northeast');
ylabel('Power (dB)'); xlabel('Range (m)');



%% CALCULATE SPECTRA OF TARGETS


singletarget.TSplot = figure('visible','off'); hold on;
    Legend = cell(1,length(singletarget.targets));

% DOES NOT INCLUDE BEAM COMPENSATION
for n = 1:length(singletarget.targets)
    
    % take section of ping in target window
    r = singletarget.targets(n).r;
    rid1 = find(singletarget.range >= r - rabove ,1);
    rid2 = find(singletarget.range >= r + rbelow ,1);
    subping = singletarget.ping(rid1:rid2);
    
    
    fsdec = 1/datastruct.param(channel,pingno).SampleInterval;
    %FreqStart = cal.FreqStart;
    %FreqEnd = cal.FreqEnd;
    nfft = 2^13;
    pingFFT = (fft(subping,nfft));
    % create frequency vector for plotting
    freq = fsdec .* linspace(0,1-1/nfft,nfft);

    fc = (FreqStart+FreqEnd)/2;
    if fc > fsdec/2
        pingFFT = [pingFFT; pingFFT; pingFFT];
        freq = [freq fsdec+freq 2*fsdec+freq];
        idxmin = round((fc-fsdec/2)/fsdec*nfft);
        pingFFT = pingFFT(idxmin:idxmin+nfft-1,:);
        freq = freq(idxmin:idxmin+nfft-1);
    end
    
    Kt2 = 2*datastruct.Z .*(1500./freq').^2 ./(16*pi^2);
    TSvec = 10*log10(abs(pingFFT).^2)...
        + singletarget.targets(n).TVG - 10*log10(Kt2*datastruct.param(1).TransmitPower);
    singletarget.targets(n).TSvec = TSvec;
    singletarget.targets(n).freq = freq;
    
    plot(freq/1000, TSvec,'LineWidth',2);
    Legend{n} = ['r = ' num2str(r) ' m'];     
end

title('Spectra of detected targets')
ylabel('TS (dB)');
xlabel('Frequency (kHz)');
xlim([FreqStart FreqEnd]/1000);
legend(Legend,'Location','southeast');

  

if plot_flag == 1
   figure(singletarget.pingplot);
   figure(singletarget.TSplot);
end

end
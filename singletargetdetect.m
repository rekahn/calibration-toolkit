%% SINGLE-TARGET DETECTION ALGORITHM

% Takes a ping and returns the spectra of all detected targets 

% input data structure

singletarget = [];

pingno = 3;
ping = CompressedVoltage(:,pingno);
pingdata = data.echodata(1,pingno).compressed;

% create options structure
singletarget.options.rmin = 10;
singletarget.options.rmax = 20;
singletarget.options.pldl = 6; % pulse length determination level (dB)
singletarget.options.pl_max = 1; % min and max pulse length (m)
singletarget.options.pl_min = 0.03;
singletarget.options.max_stdalong = 0.1; % max standard deviation in angles
singletarget.options.max_stdathw = 0.1;
singletarget.options.TSthresh = -85;
singletarget.options.max_phi = 0.5; % max off-axis angle
singletarget.options.maxd = 0.5; % max distance between pulses
% size of window for FFT calculation
singletarget.options.rabove = 0.1;
singletarget.options.rbelow = 0.15;


rix1 = find(Range > singletarget.options.rmin,1);
rix2 = find(Range > singletarget.options.rmax,1);
singletarget.range = Range(rix1:rix2);
singletarget.ping = ping(rix1:rix2);
singletarget.pingdata = pingdata(rix1:rix2,:);
singletarget.P = 10*log10(abs(singletarget.ping).^2);

% plot starting data
figure; hold on;
plot(Range,20*log10(abs(ping)),'Color',[0.5 0.5 0.5]);

plot(singletarget.range,singletarget.P);

%% Find peaks

[pks,locs] = findpeaks(singletarget.P,'minpeakwidth',3);

singletarget.pks = pks;
singletarget.locs = locs;

% structure the targets
for n = 1:length(singletarget.pks)
   singletarget.targets(n).P = pks(n);
   singletarget.targets(n).loc = locs(n);
   singletarget.targets(n).r = singletarget.range(locs(n));
   singletarget.targets(n).TVG = 40*log10(singletarget.range(locs(n))) ...
       + 2*singletarget.range(locs(n))*alpha_sea(range(locs(n)),24,20,8,70);
   singletarget.targets(n).TS = singletarget.targets(n).P + singletarget.targets(n).TVG;
end

scatter(singletarget.range(locs),pks,'filled');

%% Filter by length of pulse envelope

rm = [];
for n = 1:length(pks)
    plevel = pks(n) - singletarget.options.pldl;
    loc = singletarget.targets(n).loc;
    
    idx1 = loc - find(flip(singletarget.P(1:loc) <= plevel),1);
    idx2 = find(singletarget.P(loc:end) <= plevel , 1) + loc;
    
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
pl_idx = find([singletarget.targets.pulselength] <= singletarget.options.pl_max & ...
    [singletarget.targets.pulselength] >= singletarget.options.pl_min);
pl_idx = setdiff(pl_idx, rm);

singletarget.targets = singletarget.targets(pl_idx);

% plot
scatter(singletarget.range([singletarget.targets.loc]),[singletarget.targets.P],'filled');


%% Filter by standard deviation of angle

% calculate standard deviation of angles within the pulse envelopes

for n = 1:length(singletarget.targets)
   
    subping = singletarget.targets(n).P_envelope;
    subrange = singletarget.targets(n).r_envelope;
    [phialong, phiathw] = angle_vector(1, pingno, subrange, data);
    singletarget.targets(n).phialong = phialong;
    singletarget.targets(n).phiathw = phiathw;
    singletarget.targets(n).stdev_along = std(phialong);
    singletarget.targets(n).stdev_athw = std(phiathw);
    
end

std_idx = find([singletarget.targets.stdev_along] <= singletarget.options.max_stdalong ...
    & [singletarget.targets.stdev_athw] <= singletarget.options.max_stdathw);

singletarget.targets = singletarget.targets(std_idx);

% plot
scatter(singletarget.range([singletarget.targets.loc]),[singletarget.targets.P],'filled');


%% Filter by TS

TS_idx = find([singletarget.targets.TS] >= singletarget.options.TSthresh);

singletarget.targets = singletarget.targets(TS_idx);

% plot
scatter(singletarget.range([singletarget.targets.loc]),[singletarget.targets.P],'filled');


%% Reject targets outside max off-axis angle


max_phi = singletarget.options.max_phi;
[phi1,phi2] = angle_vector(1,pingno, singletarget.range, data);
phi = sqrt(phi1.^2 + phi2.^2);
for n = 1:length(singletarget.targets)
    singletarget.targets(n).phi = phi(singletarget.targets(n).loc);
end

phi_idx = find([singletarget.targets.phi] <= max_phi);

singletarget.targets = singletarget.targets(phi_idx);

% plot
scatter(singletarget.range([singletarget.targets.loc]),[singletarget.targets.P],'filled');


%% Reject pulses that are too close together
    % keep the one with the higher TS

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
    
    if d2-d1 <= singletarget.options.maxd
        if singletarget.targets(n1).TS > singletarget.targets(n2).TS
            rm = [rm n2];
        else
            rm = [rm n1];
        end  
    end
    
end

rm = unique(rm);
singletarget.targets(rm) = [];

% plot
scatter(singletarget.range([singletarget.targets.loc]),[singletarget.targets.P],...
    140,'filled','p','MarkerFaceColor','y','MarkerEdgeColor','k');

title('Time-series')
legend('Whole ping','Search window','All peaks','Filtered by pulse length','Filtered by StDev angle',...
    'Filtered by TS','Only on-axis','Detected targets'...
    ,'Location','southwest');
ylabel('Power (dB)'); xlabel('Range (m)');


%% Plot spectra of targets

figure; hold on;
title('Spectra of detected targets')


% DOES NOT INCLUDE BEAM COMPENSATION
for n = 1:length(singletarget.targets)
    
    % take section of ping in target window
    r = singletarget.targets(n).r;
    rid1 = find(singletarget.range >= r - singletarget.options.rabove ,1);
    rid2 = find(singletarget.range >= r + singletarget.options.rbelow ,1);
    subping = singletarget.ping(rid1:rid2);
    
    
    fsdec = 1/data.param(1,pingno).SampleInterval;
    %FreqStart = cal.FreqStart;
    %FreqEnd = cal.FreqEnd;
    nfft = 2^13;
    pingFFT = (fft(subping,nfft));
    % create frequency vector for plotting
    freq = fsdec .* linspace(0,1-1/nfft,nfft);% + fsdec; % basebanded within fc 

    fc = (FreqStart+FreqEnd)/2;
    if fc > fsdec/2
        pingFFT = [pingFFT; pingFFT; pingFFT];
        freq = [freq fsdec+freq 2*fsdec+freq];
        idxmin = round((fc-fsdec/2)/fsdec*nfft);
        pingFFT = pingFFT(idxmin:idxmin+nfft-1,:);
        freq = freq(idxmin:idxmin+nfft-1);
    end
    
    TSvec = 10*log10(abs(pingFFT).^2)...
        + singletarget.targets(n).TVG;
    singletarget.targets(n).TSvec = TSvec;
    singletarget.targets(n).freq = freq;
    
    plot(freq/1000, TSvec,'LineWidth',2);
    Legend{n} = ['r = ' num2str(r) ' m'];     
end


ylabel('TS (dB)');
xlabel('Frequency (kHz)');
xlim([FreqStart FreqEnd]/1000);
legend(Legend,'Location','southeast');



















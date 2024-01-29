% SINGLE-TARGET DETECTION ALGORITHM

% Takes a ping and returns the spectra of all detected targets 


% input: all data, ping number
% output: number of targets, spectra and range of targets
% option to specify range to look in


ping = CompressedVoltage(:,1);
rmin = 0.5; rmax = 12;
rix1 = find(Range > rmin,1);
rix2 = find(Range > rmax,1);
range = Range(rix:rix2);
ping = ping(rix:rix2);


% find peaks in pulse compressed power
P = 10*log10(abs(ping).^2);
[pks,locs] = findpeaks(P,'minpeakwidth',3); % must be 3-point peak




pldl = 25; %pulse length determination level

plength = []; pks2 = []; locs2=[];
env1 = []; env2 = [];
for n = 2:length(pks)
    plevel = pks(n) - pldl;
    i1 = locs(n) - find(flip(P(1:locs(n)) < plevel),1);
    i2 = find(P(locs(n):end) < plevel , 1) + locs(n);
    
    if ~isempty(i1) && ~isempty(i2)
        plength = [plength Range(i2) - Range(i1)];
        pks2 = [pks2 pks(n)];
        locs2 = [locs2 locs(n)];
        env1 = [env1 i1]; env2 = [env2 i2];
    end
end

% keep only targets within the right pulse range
pl_max = 1; pl_min = 0.1;
pks3 = pks2(find(plength < pl_max & plength > pl_min));
locs3 = locs2(find(plength < pl_max & plength > pl_min));

% calculate standard deviation of angles within the pulse envelopes
std_along = zeros(size(pks3)); std_athw = zeros(size(pks3));

for n = 1:length(pks3)  
    subping = ping(env1(n):env2(n)); % THESE ENVELOPES ARE WRONG
    subrange = range(env1(n):env2(n));
    [phialong,phiathw] = angle_vector(subping, subrange, data);
    std_along(n) = std(phialong);
    std_athw(n) = std(phiathw);
end

max_stdalong = 5; max_stdathw = 5;
pks4 = pks3(find(std_along <= max_stdalong & std_athw <= max_stdathw));
locs4 = locs3(find(std_along <= max_stdalong & std_athw <= max_stdathw));


% filter by TS
TS = zeros(size(pks4));
for n = 1:length(pks4)
    Pn = pks4(n);
    locn = locs(n);
    TS(n) = Pn + 40*log10(range(locn)) + ...
        2*range(locn)*alpha_sea(range(locn),24,20,8,70); 
end

TSthresh = -40;
pks5 = pks4(find(TS >= TSthresh));
locs5 = locs4(find(TS >= TSthresh));

% reject targets outside of the beamwitch of 7 degrees
phimax = 3.5;
[phi1,phi2] = angle_vector(ping, range, data);
phi = sqrt(phi1.^2 + phi2.^2);
phi5 = phi(locs5);
pks6 = pks5(find(phi5 <= phimax));
locs6 = locs5(find(phi5 <= phimax));

% rejects pulses that are too close together
%   keep the one with hightest TS
combos = combnk(1:length(pks6),2); % find all combinations of two targets

[~,ipks2,~] = intersect(pks2,pks6);
maxd = 0.5;
rm = [];
for c = 1:size(combos,1)
   
    i1 = combos(c,1); i2 = combos(c,2);
    
    p1 = range(locs6(i1));
    p2 = range(locs6(i2));
    
    if p1 < p2
        r1 = range(env2(ipks2(i1)));
        r2 = range(env1(ipks2(i2)));
    else
        r1 = range(env2(ipks2(i2)));
        r2 = range(env1(ipks2(i1)));
    end
    
    
    if abs(r2-r1) <= maxd
        TS1 = pks6(i1) + 40*log10(range(locs6(i1))) +...
            2*range(locs6(i1))*alpha_sea(range(locs6(i1)),24,20,8,70);
        TS2 = pks6(i2) + 40*log10(range(locs6(i2))) +...
            2*range(locs6(i2))*alpha_sea(range(locs6(i2)),24,20,8,70);

            if TS1 > TS2
               rm = [rm i2];
            else
               rm = [rm i1];                
            end
        
    end
    
    
end

rm = unique(rm);

locs6rm = locs6(rm);
pks6rm = pks6(rm);

locs7 = setdiff(locs6,locs6rm);
pks7 = setdiff(pks6,pks6rm);













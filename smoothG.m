% function to smooth G

function [Gsmooth] = smoothG(fvec, G, cal)

FreqStart = cal.FreqStart;
FreqEnd = cal.FreqEnd;

% find indices of end sections of the data
ind1 = find(fvec > FreqStart+7000 ,1);
ind2 = find(fvec > FreqEnd-7000, 1);

% split fvec and G into side and middle sections

fvec1 = fvec(1:ind1);
fvec2 = fvec(ind1+1:ind2);
fvec3 = fvec(ind2+1:end);

G1 = G(1:ind1);
G2 = G(ind1+1:ind2);
G3 = G(ind2+1:end);

% smooth G over each section

G1smooth = smoothdata(G1,'rloess',200);
G2smooth = smoothdata(G2,'rlowess',800); % try 'rloess' or 'movmean'
G3smooth = smoothdata(G3,'rloess',200);

% concatenate back together

Gsmooth = [G1smooth; G2smooth; G3smooth];
Gsmooth = smoothdata(Gsmooth,'rloess',300);
Gsmooth = smoothdata(Gsmooth,'movmean',300);

end

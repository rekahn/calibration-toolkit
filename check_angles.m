% check alongship and athwartship offsets


pks = zeros(size(cal.localize.phis_along));
for n = 1:size(cal.localize.compressed,2)
    
   pks(n) = max(abs(cal.localize.compressed(:,n))); 
    
end

figure; hold on;
scatter(cal.localize.phis_along,20*log10(pks))
xline(0); 
title('')
xlabel('degrees')
ylabel('max |MF|')
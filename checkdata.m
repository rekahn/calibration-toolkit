% check data to make sure it isn't faulty

% do this by making sure data from all four quadrants are of similar
% magnitude

function checkdata(cal, rmin, rmax)

npings = 9;
pings = 1:npings;
y1 = zeros(npings,1);
y2 = zeros(npings,1);
y3 = zeros(npings,1);
y4 = zeros(npings,1);

for n = 1:npings
    
    pingno = n;
    range = cal.echodata(pingno).range;
    pingdata = cal.echodata(pingno).compressed;

    % take only data in range of target
    r1 = rmin; r2 = rmax;
    ind1 = find(range > r1, 1);
    ind2 = find(range > r2, 1);
    range = range(ind1:ind2);
    pingdata = pingdata(ind1:ind2,:);

    y1(n) = mean(abs(pingdata(:,1)));
    y2(n) = mean(abs(pingdata(:,2)));
    y3(n) = mean(abs(pingdata(:,3)));
    y4(n) = mean(abs(pingdata(:,4)));
     
end

y11 = mean(y1);
y22 = mean(y2);
y33 = mean(y3);
y44 = mean(y4);
    
% raise an error if one of the channels is shot
er = 0.05; % percent error between channels allowed
if abs(y11+y22 - (y33+y44)) > 2*er *(y11+y22)
    out = 0;  
    
    % figure out which quadrant is shot
    [~,badquad] = min([y11 y22 y33 y44]);      

else
    out = 1;
end



if out == 0
    fprintf(['Looks like quadrant ' num2str(badquad) ' is shot. \n']);
else
    fprintf('All quadrants look good. \n');
end



end
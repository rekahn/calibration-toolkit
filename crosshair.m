% plot crosshairs as a colormap

function crosshair(cal)

    FreqStart = cal.FreqStart;
    FreqEnd = cal.FreqEnd;
    fc = (FreqStart+FreqEnd)/2;

    TS_data = cal.TS_data;
    f_data = cal.f_data;
    ind_fc1 = find(f_data > fc, 1);
    TS_beam = TS_data(ind_fc1-1, :) - cal.offset; 


    figure; hold on;
    yline(0,'Color','k'); xline(0,'Color','k');
    xcirc1 = 0.5*cos(0:0.01:2*pi); ycirc1 = 0.5*sin(0:0.01:2*pi);
    line(xcirc1,ycirc1,'Color','k','LineStyle','--');
    for a = 1:3
    xcirc2 = a*cos(0:0.01:2*pi); ycirc2 = a*sin(0:0.01:2*pi);
    line(xcirc2,ycirc2,'Color','k');
    end
    xlim([-3.5 3.5]); ylim([-3.5 3.5]);
    pbaspect([1 1 1]);

    
    plot1 = scatter(cal.phis_athw,cal.phis_along,20,TS_beam,'filled');
    %plot1.MarkerFaceAlpha = 0.3;
    colorbar;
    xlabel('Athwartship angle (\circ)');
    ylabel('Alongship angle (\circ)');




end
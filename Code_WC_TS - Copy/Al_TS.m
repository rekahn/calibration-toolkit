% TS for aluminum

function [TS,f]=Al_TS(D,plot_flag);
% D = Al sphere diameter in mm

echo on

%addpath c:\Code_WC_TS\toolboxes\seawater
%addpath c:\Code_WC_TS\sphere_scat

ka_min=0.01; 
ka_max=50; %set to 200 for channels other than ES70 (for 70, set to 50)

R=D./2.*1e-3;% Sphere Radius
%R=20e-3./2;% Sphere Radius
proc_flag=1;%ka
scale_flag=1;%linear
out_flag=3;%scaled scat amplitude
%cw=sw_svel(31.3,20.4,1);%cw=1490;
cw = sw_svel(35,5,sw_pres(700,39)); % salinity, temperature, pressure (depth, latitude)
%g=14900./1000;% Density of Sphere
g= 2700./1000;
%hL=6853/cw;hL=6864./cw; %Longitudinal Sound Speed of Sphere
%hT=4171/cw;hT=4161.2./cw; %Transverse Sound Speed of Sphere
hL = 6420/cw; %6420
hT = 3145/cw;%3150/cw; %3140
para_flag=[3000,ka_min,ka_max,g,hL,hT,180];
[outx,outy]=elastic_fs(proc_flag,scale_flag,out_flag,para_flag);
f=cw.*outx./2./pi./R;
TS=10.*log10(outy./2.*R.^2);

if plot_flag==1;
    %   figure(1)
    %   plot(outx,10.*log10(outy./pi./2),'k')
    %   set(gca,'linewidth',[2],'fontsize',[12])
    %   xlabel('ka','fontsize',[12])
    %   ylabel('REDUCED TARGET STRENGTH (dB)','fontsize',[12])
    %   ylim([-40 0])
    figure(2)
    hold on
    plot(f/1000,TS,'k','linewidth',2)
    set(gca,'linewidth',[2],'fontsize',[12])
    xlabel('FREQUENCY (kHz)','fontsize',[12])
    ylabel('TARGET STRENGTH (dB)','fontsize',[12])
    title('Tungsten Carbide','fontsize',[12])
    ylim([-50 -20])
    grid on
end
echo off
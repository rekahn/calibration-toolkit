% two phase conductivity

clear
n=1000;
phi0=linspace(0,1.0,n);
phi1=1-phi0;
c1=4.291;		% water
c2=0.625;		% animal
ca=c1+phi0./(1./(c2-c1)+phi1./(3*c1));          % more water
cb=c2+phi1./(1./(c1-c2)+phi0./(3*c2));		% more animal
cc=c1*(1-phi0)+phi0*c2;
%disp([phi0(:) ca(:) cb(:) cc(:) ca(:)./cb(:)  cc(:)./ca(:)])
err1=0.5*(cc-ca)./ca;
err2=0.5*(cc-cb)./cb;
plot(phi0,err1,phi0,err2)

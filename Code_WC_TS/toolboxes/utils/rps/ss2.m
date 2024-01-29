function [start,stop]=ss2(jd)
%SS2:  Gregorian start and stop of Julian day variable
%  Usage:  [start,stop]=ss2(jd)
start=gregorian(jd(1));
stop=gregorian(last(jd));
if(nargout==0),
  start(2,:)=stop
end

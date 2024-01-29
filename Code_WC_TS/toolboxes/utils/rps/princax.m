function [theta,maj,min,wr]=princax(w)
% PRINCAX Principal axis, rotation angle, principal ellipse
%
%   [theta,maj,min,wr]=princax(w) 
%
%   Input:  w   = complex vector time series (u+i*v)
%
%   Output: theta = angle of maximum variance, math notation (east == 0, north=90)
%           maj   = major axis of principal ellipse
%           min   = minor axis of principal ellipse
%           wr    = rotated time series, where real(wr) is aligned with 
%                   the major axis.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version 1.0 (12/4/1996) Rich Signell (rsignell@usgs.gov)
% Version 1.1 (4/21/1999) Rich Signell (rsignell@usgs.gov) 
%     fixed bug that sometimes caused the imaginary part 
%     of the rotated time series to be aligned with major axis. 
%     Also simplified the code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w=denan(w(:));   % remove bad values

% find covariance matrix
cv=cov([real(w(:)) imag(w(:))]); 

% find direction of maximum variance
theta=0.5*atan2(2.*cv(2,1),(cv(1,1)-cv(2,2)) );

% rotate into principal ellipse orientation
wr=w.*exp(-i*theta);
c=cov([real(wr(:)) imag(wr(:))]);
maj=sqrt(c(1,1));
min=sqrt(c(2,2));
theta=theta*180./pi;

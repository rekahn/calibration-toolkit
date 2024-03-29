function n2 = sw_bfrq(S,T,P,LAT)

% SW_BFRQ    Brunt-Vaisala Frequency Squared (N^2)
%===========================================================================
% SW_BFRQ  $Revision: 1.11 $   $Date: 1994/10/18 22:26:38 $
%          Copyright (C) CSIRO, Phil Morgan  1993. 
%
% USAGE:  bfrq = sw_bfrq(S,T,P,{LAT}) 
%
% DESCRIPTION:
%    Calculates Brunt-Vaisala Frequency squared (N^2) at the mid depths
%    from the equation,
%
%               -g      d(pdens)
%         N2 =  ----- x --------
%               pdens     d(z)
%
% INPUT:  (all must have same dimensions MxN)
%   S   = salinity    [psu      (PSS-78) ]
%   T   = temperature [degree C (IPTS-68)]
%   P   = pressure    [db]   
%
%   OPTIONAL:
%      LAT     = Latitude in decimal degrees north [-90..+90]
%                May have dimensions 1x1 or 1xn where S(mxn).
%                (Will use sw_g instead of the default g=9.8 m^2/s)
%                (Will also calc d(z) instead of d(p) in numerator)
% OUTPUT:
%   bfrq = Brunt-Vaisala Frequency squared (M-1xN)  [s^-2]
%
% AUTHOR:  Phil Morgan 93-06-24  (morgan@ml.csiro.au)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.  
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCES:
%   A.E. Gill 1982. p.54  eqn 3.7.15
%   "Atmosphere-Ocean Dynamics"
%   Academic Press: New York.  ISBN: 0-12-283522-0
%
%   Jackett, D.R. and McDougall, T.J. 1994.
%   Minimal adjustment of hydrographic properties to achieve static
%   stability.  submitted J.Atmos.Ocean.Tech.
%=========================================================================

% CALLER:  general purpose
% CALLEE:  sw_dens.m sw_pden.m

%$Id: sw_bfrq.M,v 1.11 1994/10/18 22:26:38 morgan Exp $

%-------------
% CHECK INPUTS
%-------------
if ~(nargin==3 | nargin==4) 
   error('sw_bfrq.m: Must pass 3 or 4 parameters ')
end %if

if nargin == 3
  LAT = [];
end %if

% CHECK S,T,P dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);
[mp,np] = size(P);

  
% CHECK THAT S & T HAVE SAME SHAPE
if (ms~=mt) | (ns~=nt)
   error('check_stp: S & T must have same dimensions')
end %if

% CHECK OPTIONAL SHAPES FOR P
if     mp==1  & np==1      % P is a scalar.  Fill to size of S
   P = P(1)*ones(ms,ns);
elseif np==ns & mp==1      % P is row vector with same cols as S
   P = P( ones(1,ms), : ); %   Copy down each column.
elseif mp==ms & np==1      % P is column vector
   P = P( :, ones(1,ns) ); %   Copy across each row
elseif mp==ms & np==ns     % PR is a matrix size(S)
   % shape ok 
else
   error('check_stp: P has wrong dimensions')
end %if
[mp,np] = size(P);
 

  
% IF ALL ROW VECTORS ARE PASSED THEN LET US PRESERVE SHAPE ON RETURN.
Transpose = 0;
if mp == 1  % row vector
   P       =  P(:);
   T       =  T(:);
   S       =  S(:);   

   Transpose = 1;
end %if
%***check_stp

% IF LAT PASSED THEN VERIFY DIMENSIONS
if ~isempty(LAT)
   [mL,nL] = size(LAT);
   if mL==1 & nL==1
      LAT = LAT*ones(size(S));
   end %if  

   if (ms~=mL) | (ns~=nL)              % S & LAT are not the same shape
       if (ns==nL) & (mL==1)           % copy LATS down each column
          LAT = LAT( ones(1,ms), : );  % s.t. dim(S)==dim(LAT)    
       else
          error('sw_bfrq.m:  Inputs arguments have wrong dimensions')
       end %if
   end %if
end %if

   

%------
% BEGIN
%------
if ~isempty(LAT)
   % note that sw_g expects height as argument
   Z = sw_dpth(P,LAT);
   g = sw_g(LAT,-Z);
else
   g = 9.8*ones(size(P));
end %if

[m,n] = size(P);
iup   = 1:m-1;
ilo   = 2:m;
p_ave = (P(iup,:)+P(ilo,:) )/2;
pden_up = sw_pden(S(iup,:),T(iup,:),P(iup,:),p_ave);
pden_lo = sw_pden(S(ilo,:),T(ilo,:),P(ilo,:),p_ave);
 
mid_pden = (pden_up + pden_lo )/2;
dif_pden = pden_up - pden_lo;
mid_g    = (g(iup,:)+g(ilo,:))/2;
if isempty(LAT)
   dif_z = diff(P);
else
   dif_z = diff(Z);
end %if
n2       = -mid_g .* dif_pden ./ (dif_z .* mid_pden);

if Transpose
  n2 = n2';
end %if
return
%-------------------------------------------------------------------


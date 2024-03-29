function [BETA] = sw_beta(S, PTMP, P)

% SW_BETA    Saline contraction coefficient (beta)
%========================================================================
% SW_BETA  $Revision: 1.3 $   $Date: 1994/10/10 04:33:31 $
%   %      Copyright (C) CSIRO, Nathan Bindoff 1993.
%
% USAGE:  [BETA] = sw_beta(S, PTMP, P)
%
% DESCRIPTION
%   The saline contraction coefficient as defined by T.J. McDougall.
%
% INPUT:  (all must have same dimensions)
%   S    = salinity              [psu      (PSS-78) ]
%   PTMP = potential temperature [degree C (IPTS-68)]
%   P    = pressure              [db]
%          (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
%
% OUTPUT
%   BETA = Saline Contraction Coefficient  [psu.^-1]
%
% AUTHOR:   N.L. Bindoff  1993
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.  
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCE: 
%    McDougall, T.J. 1987. "Neutral Surfaces"
%    Journal of Physical Oceanography vol 17 pages 1950-1964, 
%
% CHECK VALUE: 
%    beta=0.72088e-3 psu.^-1 at S=40.0 psu, ptmp = 10.0 C, p=4000 db
%========================================================================

% Modifications
% 93-04-22. Phil Morgan,  Help display modified to suit library
% 93-04-23. Phil Morgan,  Input argument checking

% CHECK INPUT ARGUMENTS
if nargin~=3
  error('sw_beta.m: requires 3 input arguments')
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

% BEGIN

        c1=fliplr([ 0.785567e-3, -0.301985e-5 ...
      0.555579e-7, -0.415613e-9]);
 c2=fliplr([ -0.356603e-6, 0.788212e-8]);
 c3=fliplr([0.0 0.408195e-10, -0.602281e-15]);
 c4=[0.515032e-8];
 c5=fliplr([-0.121555e-7, 0.192867e-9, -0.213127e-11]);
        c6=fliplr([0.176621e-12 -0.175379e-14]);
 c7=[0.121551e-17];
%
% Now calaculate the thermal expansion saline contraction ratio adb
%
 [m,n] = size(S);
 sm35  = S-35*ones(m,n);
 BETA  = polyval(c1,PTMP) + sm35.*(polyval(c2,PTMP) + ...
         polyval(c3,P)) + c4*(sm35.^2) + ...
         P.*polyval(c5,PTMP) + (P.^2).*polyval(c6,PTMP) ...
                +c7*( P.^3);

return
%------------------------------------------------------------------------

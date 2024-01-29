function [ALPHA] = sw_alpha(S, PTMP, P)

% SW_ALPHA   Thermal expansion coefficient (alpha)
%================================================================
% SW_ALPHA  $Revision: 1.4 $   $Date: 1994/10/10 04:14:56 $
%           Copyright (C) CSIRO, Nathan Bindoff 1993.
%
% USAGE:  [ALPHA] = alpha(S, PTMP, P)
%
% DESCRIPTION:
%    A function to calculate the thermal expansion coefficient.
%
% INPUT:
%   S    = salinity              [psu      (PSS-78) ]
%   PTMP = potential temperature [degree C (IPTS-68)]
%   P    = pressure              [db]
%          (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
%
% OUTPUT:
%   ALPHA = Thermal expansion coeff (alpha) [degree_C.^-1]
%
% AUTHOR:   N.L. Bindoff  1993
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.  
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCE: 
%    McDougall, T.J. 1987.  "Neutral Surfaces"
%    Journal of Physical Oceanography vol 17 pages 1950-1964, 
%
% CHECK VALUE:
%    See sw_beta.m amd sw_aonb.m
%================================================================

% Modifications
% 93-04-22. Phil Morgan,  Help display modified to suit library
% 93-04-23. Phil Morgan,  Input argument checking

  
% CHECK INPUT ARGUMENTS
if nargin~=3
  error('sw_alpha.m: requires 3 input arguments')
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

ALPHA = sw_aonb(S,PTMP,P).*sw_beta(S,PTMP,P);

return
%------------------------------------------------------------------------

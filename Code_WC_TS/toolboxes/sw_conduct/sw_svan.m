function svan = sw_svan(S,T,P)

% SW_SVAN    Specific volume anomaly
%=========================================================================
% SW_SVAN  $Revision: 1.3 $  $Date: 1994/10/10 05:52:20 $
%          Copyright (C) CSIRO,  Phil Morgan 1992.
%
% USAGE:  svan = sw_svan(S,T,P)
%
% DESCRIPTION:
%   Specific Volume Anomaly calculated as 
%        svan = 1/sw_dens(s,t,p) - 1/sw_dens(35,0,p)
%   Note that it is often quoted in literature as 1e8*units
%
% INPUT:  (all must have same dimensions)
%   S = salinity    [psu      (PSS-78) ]
%   T = temperature [degree C (IPTS-68)]
%   P = Pressure    [db]
%       (alternatively, may have dimensions 1*1 or 1*n where n is columns in S)
%
% OUTPUT:
%  svan = Specific Volume Anomaly  [m^3 kg^-1]
%
% AUTHOR:  Phil Morgan 92-11-05  (morgan@ml.csiro.au)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.  
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCE:
%     Fofonoff, N.P. and Millard, R.C. Jr 
%     Unesco 1983. Algorithms for computation of fundamental properties of 
%     seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
%     Eqn (9) p.15. 
%
%     S. Pond & G.Pickard  2nd Edition 1986
%     Introductory Dynamical Oceanogrpahy
%     Pergamon Press Sydney.  ISBN 0-08-028728-X
%=========================================================================

%
% CALLER: general purpose
% CALLEE: sw_dens.m

%----------------------
% CHECK INPUT ARGUMENTS
%----------------------
if nargin ~=3
   error('sw_svan.m: Must pass 3 parameters')
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

 
% -----
% BEGIN
% -----
svan = (  ones(size(S)) ./ sw_dens(S,T,P)) - ...
         (ones(size(S)) ./ sw_dens(35*ones(size(S)),zeros(size(S)),P) );

if Transpose
   svan = svan';
end %if

return
%--------------------------------------------------------------------


function PT = sw_temp(S,T,P,PR)

% SW_TEMP    Temperature from potential temperature
%===========================================================================
% TEMP  $Revision: 1.3 $  $Date: 1994/10/10 05:53:39 $
%       Copyright (C) CSIRO, Phil Morgan  1992. 
%
% USAGE:  temp = sw_temp(S,PTMP,P,PR) 
%
% DESCRIPTION:
%    Calculates temperature from potential temperature at the reference
%    pressure PR and in-situ pressure P.
%   
% INPUT:  (all must have same dimensions)
%   S     = salinity              [psu      (PSS-78) ]
%   PTMP  = potential temperature [degree C (IPTS-68)]
%   P     = pressure              [db]
%   PR    = Reference pressure    [db]
%           (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
%
% OUTPUT:
%   temp = temperature [degree C (IPTS-68)]
%
% AUTHOR:  Phil Morgan 92-04-06  (morgan@ml.csiro.au)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.  
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCES:
%    Fofonoff, P. and Millard, R.C. Jr
%    Unesco 1983. Algorithms for computation of fundamental properties of 
%    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
%    Eqn.(31) p.39
%
%    Bryden, H. 1973.
%    "New Polynomials for thermal expansion, adiabatic temperature gradient
%    and potential temperature of sea water."
%    DEEP-SEA RES., 1973, Vol20,401-408.
%=========================================================================

% CALLER:  general purpose
% CALLEE:  sw_ptmp.m

%-------------
% CHECK INPUTS
%-------------
if nargin ~= 4
   error('sw_temp.m: Must pass 4 parameters ')
end %if
% LET sw_ptmp.m DO DIMENSION CHECKING

% CARRY OUT INVERSE CALCULATION BY SWAPPING P0 & PR.
PT = sw_ptmp(S,T,PR,P);

return      
%=========================================================================

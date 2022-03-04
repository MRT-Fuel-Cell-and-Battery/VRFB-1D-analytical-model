%#Distributed under GNU Affero General Public License v3.0#
function [j0to0] = Overpotential(phi,j0)

Parameters

wo        = ko*a*exp(phi/bo);
wr        = ko*a*exp(-phi/br+E0*F/R/T);
EPS0      = (k*a/d_BL)/(wo+wr+k*a/d_BL);
wD        = D/d_el^2;
Pe        = h*d_el/D;
lamx0     = sqrt(EPS0*(wo+wr)/wD);
A2        = (SOC_ch-wr/(wo+wr))/(cosh(lamx0)+lamx0/Pe*sinh(lamx0));
j0_GUESS  = F*d_el*cvan*EPS0*(A2*sinh(lamx0)*(wo+wr)/lamx0);  
j0to0     = j0_GUESS/j0-1;
end
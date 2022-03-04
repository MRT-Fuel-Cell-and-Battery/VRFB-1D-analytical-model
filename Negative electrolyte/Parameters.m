%#Distributed under GNU Affero General Public License v3.0#
%%% DIMENSIONAL PARAMETERS (...)
R_OHM   =    0.5;           %   [Ohm/cm^2] High frequency resistance of EIS
%   Physical Constants
F       = 96485;            %   [C/mol] Faraday constant
R       = 8.314;            %   [J/mol/K] Universal ideal gas constant
%   Operating Conditions
P       = 101325;           %   [Pa] Pressure 
T       = 300;              %   [K] Temperture
SOC_ch  = 0.5;              %   [-] Electrolyte state of charge at channel
cvan    = 1.6*1e-3;         %   [mol/cm3] Vanadium total concentration
%   Kinetics
alpha   = 0.5;              %   [-] Charge transfer coefficient
ko      = 1.45*1e-2;        %   [cm/s] Kinetic rate constant
bo      = R*T/alpha/F;      %   [1/V] Tafel slope
br      = R*T/(1-alpha)/F;  %   [1/V] Tafel slope
%   Electrode
D       = 1.55*1e-5;        %   [cm^2/s] Global mass transport coefficient in electrode
a       = 3.5*1e2;          %   [1/cm] Specific surface area
EPS     = 0.8588;              %   [-] Porosity
d_el    = 0.8*280*1e-4;     %   [cm] Electrode thickness
C_dl    = 1.2*1e-2;         %   [F cm^-3] Double layer Capacitance
ACC     = 0.2;              %   [-] Accumulation term in electrode void space
A       = 25;               %   [cm^2] Electrode geometric area
%   Pore
d_BL    = 4*1e-4;           %   [cm] Pore radius
k       = 0.0018*d_BL;      %   [cm^2/s] Global mass transport coefficient in pores
%   Channel
h       = 0.09;             %   [cm^2/s] Global mass transport coefficient at channel/electrode interface
w       = 0.1;              %   [cm] Channel width
L       = A /2*w;           %   [cm] Channel length
H       = 0.1;              %   [cm] Channel height

%   Thermodynamics
E0      = -0.255;           %   [V] Equilibrium potential
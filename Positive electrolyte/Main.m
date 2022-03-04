%#Distributed under GNU Affero General Public License v3.0#
function Main(j0)

%% Variables Initialization
Parameters
jf      = 0.8; % Max simulated current
ji      = 0.005; % Min simulated current
step    = 0.005; % Simulated current step
v1      = [-jf:step:-ji];
v2      = [ji:step:jf];
j0vec   = [v1 v2];
omega   = logspace(15,-10,250);
f       = omega/2/pi;
f_MODEL = f;

%% Polarisation

kk = 1;
for jpola       = j0vec
    options     = optimset('Display', 'off');
    [phi]       = fzero( @(t) Overpotential(t,jpola), 0.002, options );
    phivec(kk)  = phi;
    kk          = kk+1;
end
I_MODEL = j0vec;
V_MODEL = -2*phivec+2*E0 - R_OHM*j0vec;

%% EIS

options = optimset('Display', 'off');
[phi]   = fzero( @(t) Overpotential(t,j0), 0.002, options );
wo      = ko*a*exp(phi/bo);
wr      = ko*a*exp(-phi/br+E0*F/R/T);
EPS0    = (k*a/d_BL)/(wo+wr+k*a/d_BL);
wD      = D/d_el^2;
Pe      = h*d_el/D;
lamr    = d_BL*sqrt(1i*omega./k);
Zw_BL   = tanh(lamr)./(lamr);
EPS1phi = 1 - (wo+wr)*(Zw_BL*d_BL/a/k)./(1+(wo+wr)*Zw_BL*d_BL/a/k);
EPS1soc = (1./cosh(lamr))./(1+(wo+wr)*d_BL/a/k*Zw_BL);
lamx0   = sqrt(EPS0*(wo+wr)/wD);
A2      = (SOC_ch-wr/(wo+wr))./(cosh(lamx0)+lamx0/Pe*sinh(lamx0));
lamx1   = sqrt(1./wD.*( EPS1soc*(wo+wr)./cosh(lamr) + lamr.*tanh(lamr)*k*a./d_BL + 1i*omega*ACC ));
C3      = EPS0*EPS1soc./wD./cosh(lamr).*A2*(wo/bo-wr/br)./(lamx0.^2-lamx1.^2);
LAMBDA  = -( EPS1phi./cosh(lamr).*(wr/(wo+wr+k*a/d_BL)*(wo/bo-wr/br)+wr/br) + EPS0*EPS1phi./cosh(lamr)*wr/(wr+wo)*(wo/bo-wr/br) )./(lamx1.^2*wD);
A3      = (-C3.*sinh(lamx0)*lamx0/Pe-C3.*cosh(lamx0)-LAMBDA)./(tanh(lamx1).*lamx1/Pe+1);
G       = F*d_el*cvan*EPS0*( EPS1soc*(wo+wr).*( A3.*tanh(lamx1)./lamx1+C3.*sinh(lamx0)./lamx0+LAMBDA ) + EPS1phi.*( (A2*sinh(lamx0)/lamx0+wr/(wo+wr))*EPS0*(wo/bo-wr/br)+(wo/bo-wr/br)*wr/(k*a/d_BL+wo+wr)+wr/br ) )+ 1i*omega*C_dl*d_el;
Z1      = G.^(-1);
Ztot    = 2*Z1 + R_OHM;

%% Plots

figure(1)
grid ON
hold on
xlabel('current density [A/cm2]')
ylabel('cell potential [V]')
axis equal
title('Polarisation curve')
ylim([-0.5;0.5])
POLA = plot(I_MODEL,V_MODEL, 'MarkerFaceColor', 'w', 'LineWidth', 2);hold on
set(gca,'fontsize', 12)
set(gca,'YTick',[-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6])
set(gca,'XTick',[-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6])
figure(2)
grid ON
hold on
xlabel('Z_{TOT,re} [Ohm cm^2]')
ylabel('Z_{TOT,im} [Ohm cm^2]')
title('Nyquist plot')
EIS = plot( real(Ztot), -imag(Ztot), 'MarkerFaceColor', 'w', 'LineWidth', 1.5);hold on
axis equal
set(gca,'fontsize', 12)
figure(3)
BODEimag = semilogx( f_MODEL, -imag(Ztot),'MarkerFaceColor', 'w',  'LineWidth', 1.5);hold on
ylabel('Z_{TOT,im} [Ohm cm^2]')
xlabel('\omega [Hz]')
title('Bode plot Imaginary')
grid ON
ylim([0,+inf])
xlim([1e-6,1e6])
set(gca,'fontsize', 12)
end
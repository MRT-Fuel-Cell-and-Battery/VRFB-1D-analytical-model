%#Distributed under GNU Affero General Public License v3.0#
function Concentrations(j0)

Parameters
options = optimset('Display', 'off');
x       = linspace (0,1,100); % non dimensional x coordinate (through electrode)
r       = linspace (0,1,100); % non dimensional r coordinate (accross boundary layer)
SOCmat0 = zeros(100); SOCmat1 = zeros(100);    % memory allocation for SOC matrix
dphi1   = j0;
omega   = 1;
[phi]   = fzero( @(t) Overpotential(t,j0), 0.002, options );
wo      = ko*a*exp(phi/bo);
wr      = ko*a*exp(-phi/br+E0*F/R/T);
EPS0    = (k*a/d_BL)/(wo+wr+k*a/d_BL);
wD      = D/d_el^2;
Pe      = h*d_el/D;

%% Through electrode (D)

% DC
lamx0   = sqrt(EPS0*(wo+wr)/wD);
A2      = (SOC_ch-wr/(wo+wr))./(cosh(lamx0)+lamx0/Pe*sinh(lamx0));
C2      = wr/(wo+wr);
SOC0x   = A2*cosh(lamx0.*x) + C2;

% AC
lamr    = d_BL*sqrt(1i*omega./k);
Zw_BL   = tanh(lamr)./(lamr);
EPS1phi = 1 - (wo+wr)*(Zw_BL*d_BL/a/k)./(1+(wo+wr)*Zw_BL*d_BL/a/k);
EPS1soc = (1./cosh(lamr))./(1+(wo+wr)*d_BL/a/k*Zw_BL);
lamx1   = sqrt(1./wD.*( EPS1soc*(wo+wr)./cosh(lamr) + lamr.*tanh(lamr)*k*a./d_BL + 1i*omega*ACC ));
C3      = EPS0*EPS1phi./wD./cosh(lamr).*A2*(wo/bo-wr/br)./(lamx0.^2-lamx1.^2);
LAMBDA  = -( EPS1phi./cosh(lamr).*(wr/(wo+wr+k*a/d_BL)*(wo/bo-wr/br)+wr/br) + EPS0*EPS1phi./cosh(lamr)*wr/(wr+wo)*(wo/bo-wr/br) )./(lamx1.^2*wD);
A3      = (-C3.*sinh(lamx0)*lamx0/Pe-C3.*cosh(lamx0)-LAMBDA)./(tanh(lamx1).*lamx1/Pe+1);
SOC1x   = A3*cosh(lamx1.*x)+C3*cosh(lamx0.*x)+LAMBDA;

%% Inside pores (k)

%DC
for kk=1:1:100
    r0BV           = cvan*EPS0*( SOC0x(kk)*(wr+wo)-wr );
    A              = r0BV/a/cvan*d_BL/k;
    B              = SOC0x(kk);
    SOCmat0(:,kk)  = -A*r+B;
end

% AC
for kk=1:1:100
    r1BV           = cvan*(EPS1soc.*SOC1x(kk)*(wo+wr)+EPS1phi*dphi1*(SOC0x(kk)*EPS0*(wo/bo-wr/br)+(wo/bo-wr/br)*wr/(k*a/d_BL+wo+wr)+wr/br));
    B1             = - r1BV/cvan*d_BL/a/k./lamr./cosh(lamr) - SOC1x(kk)*tanh(lamr);
    A1             = SOC1x(kk);
    SOCmat1(:,kk)  = A1*cosh(lamr.*r)+B1*sinh(lamr.*r);
end

%% Plot

figure (1)
plot(1-x,SOC0x,'MarkerFaceColor', 'b', 'LineWidth', 2); hold on
title('SOC along electrode thickness')
xlabel('$$\tilde{x}$$', 'Interpreter', 'LaTeX')
ylabel('SoC [-]')
grid ON
set(gca,'fontsize', 12)
figure (2)
plot(r,SOCmat0(:,end),'MarkerFaceColor', 'b', 'LineWidth', 2); hold on
title('SOC along pore radius at channel/electrode interface')
xlabel('$$\tilde{r}$$', 'Interpreter', 'LaTeX')
ylabel('SoC [-]')
grid ON
set(gca,'fontsize', 12)
figure (3)
plot(r,SOCmat0(:,1),'MarkerFaceColor', 'b', 'LineWidth', 2); hold on
xlabel('$$\tilde{r}$$', 'Interpreter', 'LaTeX')
ylabel('SoC [-]')
grid ON
title('SOC along pore radius at electrode/membrane interface')
set(gca,'fontsize', 12)
end 
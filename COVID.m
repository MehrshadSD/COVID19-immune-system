% Hancioglu et al. 2007

function f = hancioglu (t,x,pars)

% parameters
% dV/dt
gammaV    = 0.5* 510;
gammaVA   = 0.5* 619.2;
gammaVH   = 0.5* 1.02;
alphaV    = 0.5* 1.7;
aV1       = 0.5* 100;
aV2       = 23000;
% dH/dt
bHD       = 0.5* 4;
aR        = 0.5* 1;
gammaHV   = 0.5* 0.34;
bHF       = 0.5* 0.01;
% dI/dt
bIE       = 0.5* 0.066;
aI        = 0.5* 1.5;
% dM/dt
bMD       = 0.5* 1;
bMV       = 0.5* 0.0037;
aM        = 0.5* 1;
% dF/dt
bF        = 0.5* 250000;
cF        = 0.5* 2000;
bFH       = 0.5* 17;
aF        = 0.5* 8;
% dR/dt
% dE/dt
bEM       = 0.5* 8.3;
bEI       = 0.5* 2.72;
aE        = 0.5* 0.4;
% dP/dt
bPM       = 0.5* 11.5;
aP        = 0.5* 0.4;
% dA/dt
bA        = 0.5* 0.043;
gammaAV   = 0.5* 146.2;
aA        = 0.5* 0.043;
% dS/dt
r         = 0.5* 3e-5;

% extract variables
V = x(1);
H = x(2);
L = x(3);
I = x(4);
M = x(5);
F = x(6);
R = x(7);
E = x(8);
P = x(9);
A = x(10);
S = x(11);
D = 1-H-R-I-L;

% rate equations
dVdt = gammaV*I - gammaVA*S*A*V - gammaVH*H*V - alphaV*V - aV1*V/(1+aV2*V);
dHdt = bHD*D*(H+R) + aR*R - gammaHV*V*H - bHF*F*H;
dLdt = gammaHV*V*H - 6*L;
dIdt = 6*L - bIE*E*I - aI*I;
dMdt = (bMD*D + bMV*V)*(1-M) - aM*M;
dFdt = bF*M + cF*I - bFH*H*F - aF*F;
dRdt = bHF*F*H - aR*R;
dEdt = bEM*M*E - bEI*I*E + aE*(1-E);
dPdt = bPM*M*P + aP*(1-P);
dPdt = bPM*M*P + aP*(1-P);
dAdt = bA*P - gammaAV*S*A*V - aA*A;
dSdt = r*P*(1-S);

f = [dVdt; dHdt; dLdt; dIdt; dMdt; dFdt; dRdt; dEdt; dPdt; dAdt; dSdt];
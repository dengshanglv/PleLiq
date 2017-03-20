clear;
cd ~/Desktop/PleLiq;
syms x h l q k positive;
syms I positive;
alpha = 0.33;
beta = 0.95;
pi = 0.037;

r = 1/beta - 1;
irate = 1/beta*(1+pi) - 1;
deltah = 0.1;
deltak = 0.1;
taul = 0.3;
taua = 0.4;

psi = 0.86;
zeta = 0.38;
xi = 0.3;
mub = 0.45;
muk = 0.17;
muh = 0.06;
B = 0.08;
theta = 0.68;
eta = 0.63;
upsi = 1.39;
sigma = 0.67;

ucm = log(x)+psi*log(h)-zeta*l;
f = k^alpha*l^(1-alpha);
udm = upsi*q^(1-eta)/(1-eta);
gi = I^(1+xi)/(1+xi);
cdm = q;
zq = (1-theta)*udm + theta*cdm;

ucmx = diff(ucm, x);
ucmh = diff(ucm, h);
gi1 = diff(gi, I);
fk = diff(f, k);
fl = diff(f, l);
udm1 = diff(udm, q);
zq1 = diff(zq, q);
Lq = udm1/zq1 - 1;
qstar = vpasolve(udm1-zq1==0, q);
qi = vpasolve(sigma*Lq==irate, q);

%% 1. liquid monetary eqm.
eq1 = ucmx - zeta/(1-taul)/fl;
eq2 = 1/beta - 1 - (fk - deltak)*(1 - taua);
eq3 = subs(gi1,I, deltah*h) - beta*fl*(1-taul)*ucmh/(zeta*(1-beta*(1-deltah)));
eq4 = x+subs(gi,I,deltah*h)+deltak*k-f;
eq = [eq1, eq2, eq3, eq4];
[sx,sh,sl,sk] = vpasolve(eq, [x,h,l,k], [ 1.6295337485524106480388388713796, 9.3642916995384207462149818383245, 2.1464840778441239745193076445164, 4.9820353238622916142009358789843]);
s = [sx,sh,sl,sk];
dstar1 = subs(zq*(1-taul)*fl/zeta, [q,k,l],[qstar,sk,sl]);

phib = beta*(1-taua)/(1-beta*taua);
phih = subs(gi1, deltah*sh);
Dbar = subs((1-(1-phib)*taua)*mub*B + (1+(fk-deltak)*(1-taua))*muk*k + ...
    (1-deltah)*muh*phih*h,[h,l,k],[sh,sl,sk]);
% In this setting, dstar1 is larger than Dbar, so no case 1.

%% 2. 

%% 3. 
equ1 = beta*fl*(1-taul)*ucmh/zeta+beta*(1-deltah)*(1+muh*irate) - subs(gi1,I,deltah*h); 
equ2 = (1+(fk-deltak)*(1-taua))*(1+muk*irate) - 1/beta;
equ3 = f - subs(gi, I, deltah*h) - deltak*k - x;
equ4 = ucmx - zeta/(1-taul)/fl;
equ = [equ1, equ2, equ3, equ4];
[sx3,sh3,sl3,sk3] = vpasolve(equ, [ 1.758490013597512746354443843519, 10.07620519731924598101579499508, 2.2388501820413840041916892821563, 6.5453989773017054818494474287658]);
s3 = [sx3,sh3,sl3,sk3];

phib3 = (1-taua)*(1+mub*irate)/(1+r-taua*(1+mub*irate));
phih3 = subs(gi1, deltah*sh3);
Dbar3 = subs((1-(1-phib3)*taua)*mub*B + (1+(fk-deltak)*(1-taua))*muk*k + ...
    (1-deltah)*muh*phih3*h,[h,l,k],[sh3,sl3,sk3]);

di3 = subs(zq*(1-taul)*fl/zeta, [q,k,l],[qi,sk3,sl3]);



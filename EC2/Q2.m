clear
nusp = 13683786;

mnp = rem(nusp,1000)/100;
Z0 = 75; u = 2e8;
Rg = 120 + mnp;
ZL = Rg;
rhol=(ZL-Z0)/(ZL+Z0);
f=100e6;
lambda=u/f;
l=4.4;
Eg = 15/sqrt(2); % Vef

rhoent=rhol*exp(-1j*4*pi*l/lambda);
Zent=Z0*(1+rhoent)/(1-rhoent);
Vent=Eg*Zent/(Zent+Rg);
% item a
disp '|Vent| (V) ='; disp(abs(Vent)*sqrt(2));
VL=Vent*(1+rhol)/(1+rhoent);
disp '|VL| (V) ='; disp(abs(VL)*sqrt(2));
% item c
disp 'Zent (ohm) ='; disp(Zent);
disp '|Zent| (ohm) ='; disp(abs(Zent));
disp 'arg(Zent) (rad) ='; disp(arg(Zent)*180/pi);

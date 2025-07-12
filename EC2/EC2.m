clear
nusp = 786;

format long;

disp 'Questão 1';

mnp = rem(nusp,1000)/10;
ZL = 400 + mnp;
Z0 = 75; u = 2e8;
rhol=(ZL-Z0)/(ZL+Z0);
f=100e6;
lambda=u/f;
l=50;
Eg = 15; % Vef
Rg = 75;

% item a
disp ' '; disp 'item a';
Pd = Eg^2/(4*Rg);
disp 'Pd (W) ='; disp (Pd);

%itens b, c
d = 0.18*lambda; % dist min, para aprox. 227 ohms
yd = Z0/ZL; % inicial = yL = YL*Z0
% yL = Z0/ZL;

while (real(yd) < 1)
  d = d + .000001;
  rhod=rhol*exp(-1j*4*pi*d/lambda);
  %Zd=Z0*(1+rhod)/(1-rhod);
  %yd = Z0/Zd;
  yd = (1-rhod)/(1+rhod);
  if (d > 0.22*lambda) % dist max, para aprox. 667 ohms
    break;
  end
end
disp ' '; disp 'item b';
disp 'd (m) ='; disp(d);
disp ' '; disp 'item c';
disp 'yd (S/S (adimensional)) ='; disp(yd);
disp 'b (S/S (adimensional)) ='; disp(imag(yd));

% item d
disp ' '; disp 'item d';
b = imag(yd);
arg_rho_t = angle((1/(-1j*b) - 1)/(1/(-1j*b) + 1)); % *180/pi
lt = lambda/4*(1-arg_rho_t/pi);
disp 'lt (m) ='; disp(lt);

% item e
f=800:1200; f=f*1e5;
lambda=u./f;
% e-i
rhoent=rhol*exp(-1j*4*pi*l./lambda);
Zent=Z0.*(1+rhoent)./(1-rhoent);
Ient=Eg./(Zent+Rg); Ient=abs(Ient);
Pent=real(Zent).*Ient.^2; % Sem perdas, P_ent = P_L
figure(2);
subplot(2,2,1);
plot(f/1e6,10*log10(Pent/Pd));
axis([-inf inf -4 -3]) % se o primeiro gráfico n aparecer nada, comente essa linha
xlabel('f (MHz)');
ylabel('P_L (dB)');
title('Sem toco');
disp ' '; disp 'item e i: Sem toco';
disp 'P_max (W) ='; disp(max(Pent));
disp 'BandWidth sem toco (MHz) = 0';

% e-ii
rhod=rhol*exp(-1j*4*pi*d./lambda);
rhot=-1*exp(-1j*4*pi*lt./lambda);
b=(1-rhot)./(1+rhot);
Yd = ((1-rhod)./(1+rhod))./Z0 + b./Z0;
rhod=(1./Yd-Z0)./(1./Yd+Z0);
rhoent=rhod.*exp(-1j*4*pi.*(l-d)./lambda);
Zent=Z0.*(1+rhoent)./(1-rhoent);
Ient=Eg./(Zent+Rg); Ient=abs(Ient);
Pent=real(Zent).*Ient.^2; % Carga casada, P_ent = P_L
ganho=10*log10(Pent./Pd);
subplot(2,2,2);
plot(f/1e6,ganho);
xlabel('f (MHz)');
ylabel('P_L (dB)');
title('Com toco, d_{min}');
disp ' '; disp 'item e ii: Com toco';
disp 'P_max (W) ='; disp(max(Pent));
disp 'BW com toco (MHz, MHz) =';
BW = zeros(1,2);
j=1;
for k=1:size(f)(2)
  if (ganho(k) > -2)
      if (j==1)
          BW(j) = f(k); j = j+1;
      else if (k==size(f)(2))
          BW(2)=f(k);
          disp 'Fora da faixa pedida, aumentar a banda de frequencias (linha 52)';
          end end
  else
      if (j==2)
        BW(j) = f(k); j = j+1; end
  end
end
disp (BW*1e-6); disp 'BandWidth (MHz) ='; disp ((BW(2) - BW(1))*1e-6);

% e iii
d = d + 3;
rhod=rhol*exp(-1j*4*pi*d./lambda);
Yd = ((1-rhod)./(1+rhod))./Z0 + b./Z0;
rhod=(1./Yd-Z0)./(1./Yd+Z0);
rhoent=rhod.*exp(-1j*4*pi.*(l-d)./lambda);
Zent=Z0.*(1+rhoent)./(1-rhoent);
Ient=Eg./(Zent+Rg); Ient=abs(Ient);
Pent=real(Zent).*Ient.^2; % Carga casada, P_ent = P_L
ganho=10*log10(Pent./Pd);
subplot(2,2,3);
plot(f/1e6,ganho);
xlabel('f (MHz)');
ylabel('P_L (dB)');
title('Com toco, d_{min} + 3m');
BW2 = zeros(1,2);
j=1;
for k=1:400
  if (ganho(k) > -2)
      if (j==1)
          BW2(j) = f(k); j = j+1; end
  else
      if (j==2)
          BW2(j) = f(k); j = j+1; end
  end
end
disp ' '; disp 'item e iii: Com toco + 3m';
disp 'P_max (W) ='; disp(max(Pent));
disp 'BW com toco + 3m (MHz, MHz) ='; disp (BW2*1e-6); disp 'BandWidth (MHz) ='; disp ((BW2(2) - BW2(1))*1e-6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Questão 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp ' '; disp ' ';
disp 'Questão 2';
clearvars -except nusp

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
disp ' '; disp 'item a, valores de pico';
disp '|Vent| (V) ='; disp(abs(Vent)*sqrt(2));
VL=Vent*(1+rhol)/(1+rhoent);
disp '|VL| (V) ='; disp(abs(VL)*sqrt(2));
% item c
disp ' '; disp 'item c';
disp 'Zent (ohm) ='; disp(Zent);
disp '|Zent| (ohm) ='; disp(abs(Zent));
disp 'arg(Zent) (graus) ='; disp(arg(Zent)*180/pi);

clear
nusp=13683786;

% PTC3314 - Ondas e Linhas: Exercicio-programa 2
% Guilherme Fortunato Miranda - NUSP: 13683786
%
% pr�-defini��es
epsilon0 = 8.854e-12;
mi0 = 4e-7*pi;
mnp=mod(nusp,1000)/1000;
%
% c�lculos (item a)
miar = mi0;                                                                             % dado
epsilonar = epsilon0;                                                                   % dado
mividro = mi0;                                                                          % dado
epsilonvidro = (3+mnp)*epsilon0;                                                          % dado
theta1 = [0:(pi/1800):pi/2];                                                            % alterado para melhorar a precis�o


for aux = [1:901]
    if (sin(theta1(aux)))*sqrt(epsilonvidro/epsilonar) <= 1
        theta2(aux) = asin((sin(theta1(aux)))*sqrt(epsilonvidro/epsilonar));            % theta2 = asin(sin(theta1)*sqrt(epsilon1/epsilon2))
    else
        theta2(aux) = pi/2 + i*acosh((sin(theta1(aux)))*sqrt(epsilonvidro/epsilonar));  % theta2 = pi/2 + i*acosh((sin(theta1))*sqrt(epsilon1/epsilon2))
    end
end
x = theta1*180/pi;                                                                      % vetor abcissas
y1 = real(theta2*180/pi);                                                               % vetor ordenadas 1
y2 = imag(theta2*180/pi);                                                               % vetor ordenadas 2
subplot(2,2,1), plot(x,y1,x,y2)
title('\theta_2 \times \theta_1')
ylabel('Ângulo da onda transmitida (°)')
xlabel('Ângulo da onda incidente (°)')
legend('real','imaginário')
grid
%
% c�lculos (item b)
etaar = sqrt(miar/epsilonar);                                                           % eta = sqrt(mi/epsilon) (meio sem perdas)
etavidro = sqrt(mividro/epsilonvidro);                                                  % eta = sqrt(mi/epsilon) (meio sem perdas)
%1) Hy2(x,z) = (1-ro)*(Hy+)*exp(-i*(betax*x+betaz2*z)) (z>0)
%2) Ex2(x,z) = eta2*cos(theta2)*(1-ro)*(Hy+)*exp(-i*(betax*x+betaz2*z)) (z>0)
%-> Ex2/Hy2 = eta2*cos(theta2) = ZL
ZL = etaar*cos(theta2);                                                                 % ZL = eta2*cos(theta2)
Zz1 = etavidro*cos(theta1);                                                             % Zz1 = eta1*cos(theta1)
ro0 = (ZL - Zz1)./(ZL + Zz1);                                                           % ro0 = (ZL - Zz1)/(ZL + Zz1)
x = theta1*180/pi;                                                                      % vetor abcissas
y1 = real(ZL);                                                                          % vetor ordenadas 1
y2 = imag(ZL);                                                                          % vetor ordenadas 2
subplot(2,2,2), plot(x,y1,x,y2)
title('Impedância Z_L \times \theta_1')
ylabel('Impedância Z_L (ohms)')
xlabel('Ângulo da onda incidente (°)')
legend('real','imaginário')
grid
%
% c�lculos (item c)
%1) Hy1(x,z) = (H1y1+)*exp(-i*betax1*x)*(exp(-i*betaz1*z)-ro0*exp(i*betaz1*z)) (z<0)
%2) Ex1(x,z) = eta1*cos(theta1)*(H1y1+)*exp(-i*betax1*x)*(exp(-i*betaz1*z)+ro0*exp(i*betaz1*z)) (z<0)
%3) Ez1(x,z) = -eta1*sin(theta1)*(H1y1+)*exp(-i*betax1*x)*(exp(-i*betaz1*z)-ro0*exp(i*betaz1*z)) (z<0)
%4) Hy2(x,z) = (1-ro)*(Hy1+)*exp(-i*(betax2*x+betaz2*z)) (z>0)
%5) Ex2(x,z) = eta2*cos(theta2)*(1-ro)*(Hy1+)*exp(-i*(betax2*x+betaz2*z)) (z>0)
%6) Ez2(x,z) = -eta2*sin(theta2)*(1-ro)*(Hy1+)*exp(-i*(betax2*x+betaz2*z)) (z>0)
%-> Nz1+ = Re((Ex1+) x (Hy1*+))/2 = (eta1*cos(theta1)/2)*(abs(Hy1+))^2
%-> Nz1- = Re((Ex1-) x (Hy1*+))/2 = -(eta1*cos(theta1)/2)*((abs(Hy1+))^2)*(abs(ro0))^2
%-> Nz2+ = Re(Ex2 x Hy2*)/2 = (eta2/2)*Re(cos(theta2))*((abs(1-ro0))^2)*((abs(Hy1+))^2)*(abs(exp(-i*betaz2*z)))^2
%-> |Nz2+|/|Nz1+| = (eta2/eta1)*(abs(Re(cos(theta2)))*((abs(1-ro0))^2)/(cos(theta1))
%-> |Nz1-|/|Nz1+| = -abs(ro0)^2
N2pdN1p = ((etaar/etavidro)*real(cos(theta2)).*(abs(ones(1,901)-ro0).^2))./(cos(theta1));
N1ndN1p = -abs(ro0).^2;
x = theta1*180/pi;                                                                      % vetor abcissas
y1 = abs(N2pdN1p);                                                                      % vetor ordenadas 1
y2 = abs(N1ndN1p);                                                                      % vetor ordenadas 2
subplot(2,2,3), plot(x,y1,x,y2)
title('Relações |N_{z2}^+| / |N_{z1}^+| e |N_{z1}^-| / |N_{z1}^+|')
ylabel('|N_{z2}^+| / |N_{z1}^+| e |N_{z1}^-| / |N_{z1}^+|')
xlabel('Ângulo da onda incidente (°)')
legend('|N_{z2}^+| / |N_{z1}^+|','|N_{z1}^-| / |N_{z1}^+|')
grid
%
% c�lculos (item d)
kar = 2*pi*1e+9*sqrt(miar*epsilonar);                                                   % k = omega*sqrt(mi*epsilon) (meio sem perdas)
kvidro = 2*pi*1e+9*sqrt(mividro*epsilonvidro);                                          % k = omega*sqrt(mi*epsilon) (meio sem perdas)
betaxar = kar*sin(theta2);                                                              % betax1 = k1*sin(theta1)
betazar = kar*cos(theta2);                                                              % betaz1 = k1*cos(theta1)
betaxvidro = kvidro*sin(theta1);                                                        % betax2 = k2*sin(theta2)
betazvidro = kvidro*cos(theta1);                                                        % betaz2 = k2*cos(theta2)
z = [-0.2:0.001:0.1];                                                                   % dado
% Hy1(x,z) = (Hy1+)*exp(-i*betax1*x)*(exp(-i*betaz1*z)-ro0*exp(i*betaz1*z)) (z<0)
% Hy2(x,z) = (1-ro)*(Hy1+)*exp(-i*(betax1*x+betaz2*z)) (z>0)
% |Hy1| = |Hy1+|*|1-ro0*exp(i*2*betaz1*z)| (z<0)
% |Hy2| = |Hy1+|*|1-ro0|*|exp(-i*betaz2*z)| (z>0)
for aux = [1:301]
    if z(aux)<=0
        ModHy(aux) = 1*abs(1-ro0(601)*exp(i*2*betazvidro(601)*z(aux)));
    else
        ModHy(aux) = 1*abs(1-ro0(601))*abs(exp(-i*betazar(601)*z(aux)));
    end
end
x = z*100;                                                                                  % vetor abcissas
y = ModHy;                                                                              % vetor ordenadas
subplot(2,2,4), plot(x,y)
title('|H| \times z')
ylabel('módulo de H (mA_{ef}/m)')
xlabel('distância ao plano onde a onda incide (cm)')
grid

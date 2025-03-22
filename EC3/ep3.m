clear
nusp=13683786;
% PTC3314 - Ondas e Linhas: Exercicio-programa 3
% Guilherme Fortunato Miranda - NUSP: 13683786
% Prof. Luis Cezar Trintinalia
%
% Parte 1:
%
%  pr�-defini��es
mnp=mod(nusp,1000)/1000;
epsilon0 = 8.854e-12;
mi0 = 4e-7*pi;
%
%  c�lculos (item a)
disp 'item a';
f0 = 1+mnp; f0=f0*1e+9;                                                                 % dado
mi1 = mi0;                                                                     % considerado
mi2 = mi0;                                                                     % considerado
mi3 = mi0;                                                                     % considerado
epsilon1 = epsilon0;                                                           % dado
epsilon3 = 16*epsilon0;                                                         % dado
eta1 = sqrt(mi1/epsilon1);                                                     % eta = sqrt(mi/epsilon) (meio sem perdas)
eta3 = sqrt(mi3/epsilon3);                                                     % eta = sqrt(mi/epsilon) (meio sem perdas)
eta2 = sqrt(eta1*eta3);                                                        % dado
epsilon2 = mi2/(eta2^2);                                                        % epsilon = mi/(eta^2) (meio sem perdas)
disp 'epsilon2 (epsilon0) ='; disp (epsilon2/epsilon0);

k1 = 2*pi*f0*sqrt(mi1*epsilon1);                                               % k = omega*sqrt(mi*epsilon) (meio sem perdas)
k2 = 2*pi*f0*sqrt(mi2*epsilon2);                                               % k = omega*sqrt(mi*epsilon) (meio sem perdas)
k3 = 2*pi*f0*sqrt(mi3*epsilon3);                                               % k = omega*sqrt(mi*epsilon) (meio sem perdas)
lambda1 = 2*pi/k1;                                                             % lambda = 2*pi/k
lambda2 = 2*pi/k2;                                                             % lambda = 2*pi/k
lambda3 = 2*pi/k3;                                                             % lambda = 2*pi/k
d = lambda2/4;                                                                  % coef. de reflex�o nulo: eta adequado, d = (2*n + 1)*lambda/4
disp 'd (mm) ='; disp (d*10^3);

%
%  c�lculos (item b)
disp ' '; disp 'item b';
a = f0 - 0.5e+9;                                                               % freq��ncia inicial
b = f0 + 0.5e+9;                                                               % freq��ncia final
p = 0.5e+6;                                                                     % passo
x = [a:p:b];                                                                   % vetor abcissas
k1 = 2*pi*x*sqrt(mi1*epsilon1);                                                % k = omega*sqrt(mi*epsilon) (meio sem perdas)
k2 = 2*pi*x*sqrt(mi2*epsilon2);                                                % k = omega*sqrt(mi*epsilon) (meio sem perdas)
k3 = 2*pi*x*sqrt(mi3*epsilon3);                                                % k = omega*sqrt(mi*epsilon) (meio sem perdas)
Z3 = eta3;                                                                     % ZN = etaN (�ltimo meio)
gama23 = (Z3 - eta2)/(Z3 + eta2);                                              % gamaN,N+1 = (ZN+1 - etaN)/(ZN+1 + etaN)
Z2 = eta2*(1 + gama23*exp(-i*2*k2*(-d)))./(1 - gama23*exp(-i*2*k2*(-d)));      % ZN(z) = etaN*(1 + gamaN,N+1(z))/(1 - gamaN,N+1(z))
gama12 = (Z2 - eta1)./(Z2 + eta1);                                             % gamaN,N+1 = (ZN+1 - etaN)/(ZN+1 + etaN)
mgamaq23 = (abs(gama23))^2;                                                    % m�dulo ao quadrado de gama
mgamaq12 = (abs(gama12)).^2;                                                   % m�dulo ao quadrado de gama
y = mgamaq12;                                                                  % vetor ordenadas
figure 1;
plot(x*1e-9,y)
grid on;
axis([1.2 2.41 -0.00054 0.101]);
title('|\rho_l|^2 \times f para casamento com camada única')
ylabel('módulo ao quadrado coeficiente de reflexão')
xlabel('frequência (GHz)')
disp '|rho_l|^2 (f0-500MHz) ='; disp (y(1));
disp '|rho_l|^2 (f0) ='; disp (min(y));
disp '|rho_l|^2 (f0+500MHz) ='; disp (y(end));
%
%  c�lculos (item c)
disp ' '; disp 'item c';
n1 = 1001;                                                                       % ponto de freq��ncia 1,6786Hz
while mgamaq12(n1) <= 0.001,
    n1 = n1-2;                                                                 % ponto de freq��ncia m�nima da banda - 1MHz (2*p)
end
n2 = 1001;
while mgamaq12(n2) <= 0.001,
    n2 = n2+2;                                                                 % ponto de freq��ncia m�xima da banda + 1 MHz (2*p)
end
B = [x(n1+2),x(n2-2)];
deltaB = x(n2-2) - x(n1+2);
disp 'B (GHz) ='; disp (B*1e-9); disp 'deltaB (MHz) ='; disp (deltaB*1e-6);


%
% Parte 2:
%
clearvars -except mnp
disp ' '; disp 'Parte 2';
%  pr�-defini��es
epsilon0 = 8.854e-12;
mi0 = 4e-7*pi;
%
%  c�lculos (item a)
disp 'item a';
f0 = 1+mnp; f0=f0*1e+9;                                                                 % dado
mi1 = mi0;                                                                     % considerado
mi2 = mi0;                                                                     % considerado
mi3 = mi0;                                                                     % considerado
mi4 = mi0;                                                                     % considerado
epsilon1 = epsilon0;                                                           % dado
epsilon4 = 16*epsilon0;                                                         % dado
eta1 = sqrt(mi1/epsilon1);                                                     % eta = sqrt(mi/epsilon) (meio sem perdas)
eta4 = sqrt(mi4/epsilon4);                                                     % eta = sqrt(mi/epsilon) (meio sem perdas)
eta2 = sqrt(sqrt((eta1^3)*eta4));                                              % dado
eta3 = sqrt(sqrt(eta1*(eta4^3)));                                              % dado
epsilon2 = mi2/(eta2^2);                                                        % epsilon = mi/(eta^2) (meio sem perdas)
epsilon3 = mi3/(eta3^2);                                                        % epsilon = mi/(eta^2) (meio sem perdas)
k1 = 2*pi*f0*sqrt(mi1*epsilon1);                                               % k = omega*sqrt(mi*epsilon) (meio sem perdas)
k2 = 2*pi*f0*sqrt(mi2*epsilon2);                                               % k = omega*sqrt(mi*epsilon) (meio sem perdas)
k3 = 2*pi*f0*sqrt(mi3*epsilon3);                                               % k = omega*sqrt(mi*epsilon) (meio sem perdas)
k4 = 2*pi*f0*sqrt(mi4*epsilon4);                                               % k = omega*sqrt(mi*epsilon) (meio sem perdas)
lambda1 = 2*pi/k1;                                                             % lambda = 2*pi/k
lambda2 = 2*pi/k2;                                                             % lambda = 2*pi/k
lambda3 = 2*pi/k3;                                                             % lambda = 2*pi/k
lambda4 = 2*pi/k4;                                                             % lambda = 2*pi/k
d2 = lambda2/4;                                                                 % coef. de reflex�o nulo: eta adequado, d = (2*n + 1)*lambda/4
d3 = lambda3/4;                                                                 % coef. de reflex�o nulo: eta adequado, d = (2*n + 1)*lambda/4
disp 'epsilon2 (epsilon0) ='; disp (epsilon2/epsilon0);
disp 'epsilon3 (epsilon0) ='; disp (epsilon3/epsilon0);
disp 'd2 (mm) ='; disp (d2*1e3);
disp 'd3 (mm) ='; disp (d3*1e3);
%
%  c�lculos (item b)
disp ' '; disp 'item b';
a = f0 - 0.5e+9;                                                               % freq��ncia inicial
b = f0 + 0.5e+9;                                                               % freq��ncia final
p = 0.5e+6;                                                                     % passo
x = [a:p:b];                                                                   % vetor abcissas
k1 = 2*pi*x*sqrt(mi1*epsilon1);                                                % k = omega*sqrt(mi*epsilon) (meio sem perdas)
k2 = 2*pi*x*sqrt(mi2*epsilon2);                                                % k = omega*sqrt(mi*epsilon) (meio sem perdas)
k3 = 2*pi*x*sqrt(mi3*epsilon3);                                                % k = omega*sqrt(mi*epsilon) (meio sem perdas)
k4 = 2*pi*x*sqrt(mi4*epsilon4);                                                % k = omega*sqrt(mi*epsilon) (meio sem perdas)
Z4 = eta4;                                                                     % ZN = etaN (�ltimo meio)
gama34 = (Z4 - eta3)/(Z4 + eta3);                                              % gamaN,N+1 = (ZN+1 - etaN)/(ZN+1 + etaN)
Z3 = eta3*(1 + gama34*exp(-i*2*k3*(-d3)))./(1 - gama34*exp(-i*2*k3*(-d3)));    % ZN(z) = etaN*(1 + gamaN,N+1(z))/(1 - gamaN,N+1(z))
gama23 = (Z3 - eta2)./(Z3 + eta2);                                             % gamaN,N+1 = (ZN+1 - etaN)/(ZN+1 + etaN)
Z2 = eta2*(1 + gama23.*exp(-i*2*k2*(-d2)))./(1 - gama23.*exp(-i*2*k2*(-d2)));  % ZN(z) = etaN*(1 + gamaN,N+1(z))/(1 - gamaN,N+1(z))
gama12 = (Z2 - eta1)./(Z2 + eta1);                                             % gamaN,N+1 = (ZN+1 - etaN)/(ZN+1 + etaN)
mgamaq34 = (abs(gama34))^2;                                                    % m�dulo ao quadrado de gama
mgamaq23 = (abs(gama23)).^2;                                                   % m�dulo ao quadrado de gama
mgamaq12 = (abs(gama12)).^2;                                                   % m�dulo ao quadrado de gama
y = mgamaq12;                                                                  % vetor ordenadas
figure(2);
plot(x*1e-9,y)
grid on;
axis([1.2 2.41 -0.00054 0.021]);
title('|\rho_l|^2 \times f para casamento com duas camadas')
ylabel('módulo ao quadrado coeficiente de reflexão')
xlabel('frequência (GHz)')
disp '|rho_l|^2 (f0-500MHz) ='; disp (y(1));
disp '|rho_l|^2 (f0) ='; disp (min(y));
disp '|rho_l|^2 (f0+500MHz) ='; disp (y(end));
%
%  c�lculos (item c)
disp ' '; disp 'item c';
n1 = 1001;                                                                       % ponto de freq��ncia 1,645GHz
while mgamaq12(n1) <= 0.001,
    n1 = n1-2;                                                                 % ponto de freq��ncia m�nima da banda - 1
end
n2 = 1001;
while mgamaq12(n2) <= 0.001,
    n2 = n2+2;                                                                 % ponto de freq��ncia m�xima da banda + 1
end
B = [x(n1+2),x(n2-2)];
deltaB = x(n2-2) - x(n1+2);
disp 'B (GHz) ='; disp (B*1e-9); disp 'deltaB (MHz) ='; disp (deltaB*1e-6);


% spin current for NM-QD-MI system
% 2020-11-04
% by ligy

%clc;
close all;
clear all;

% constants
k_B = physconst('Boltzman');
charge_e = 1.602176634e-19;
meV = 1.0e-3.*charge_e;

global omegac alpha E0_up E0_down Gamma0

omegac = 80;                            % cutoff frequency, unit: meV
alpha = 0.2;                            % disspation strength, dimensionless
T0 = 300;                               % average temperature, unit:K
mu0 = 20;                               % average spin baias, unit: meV
E0_up = 3;                              % QD up level, unit: meV
E0_down = -2;                           % QD down level, unit: meV
Gamma0=1;                               % effective coupling, unit:meV

% integral limits of E
E_limit = 5e2;
E_lower = -1.*E_limit;
E_upper = E_limit;

%==================================================================================
% % output
% display('Warning! the integral limit is [-2k_B*T0, 2k_B*T0], [0, 80]');
%==================================================================================
%calculate current
dT = linspace(-2*T0, 2*T0, 1000);
d_mu = linspace(-2*mu0, 2*mu0, 1000);
data = zeros(100, 3);
fileID = fopen('current.txt','w');
count = 1;
for deltaT = dT
    for delta_mu = d_mu                 % difference of spin voltage bias, unit: meV
        T_L = T0 + deltaT./2.0;               % left lead temperature, unit: K
        T_R = T0 - deltaT./2.0;              % right lead temperature, unit: K
        beta_L = 1./(k_B.*T_L./meV);    % beta of left metal lead, unit: meV^-1
        beta_R = 1./(k_B.*T_R./meV);      % beta of right MI lead, unit: meV^-1

        mu_up = mu0 + delta_mu./2.0;    % spin-up chemical, unit: meV
        mu_down = mu0 - delta_mu./2.0;  % spin-down chemical, unit: meV

        % Bosonic distribution
        N_L = @(omega) 1./(exp(beta_L.*(omega + delta_mu)) - 1);
        N_R = @(omega) 1./(exp(beta_R.*(omega + delta_mu)) - 1);
        
        % Fermionic distribution
        f_L_up = @(E) 1./(exp(beta_L.*(E-mu_up))+1);
        f_L_down = @(E) 1./(exp(beta_L.*(E-mu_down))+1);

        % define integrant
        my_integrant2 = @(E, omega) rho(omega) .* (N_R(omega) - N_L(omega)) .* (f_L_up(E) - f_L_down(E+omega)) .* A(E, omega);

        currt = quad2d(my_integrant2, E_lower, E_upper, 0, omegac);
        count = count + 1;
        fprintf(fileID, '%-15.10g%-15.10g%-15.10g\n',deltaT, delta_mu, currt);
    end
end
display(count);
%==================================================================================
% % plot contour, save to file
% fig = figure;
% axes;
% set(fig, 'InvertHardcopy', 'off');
% contour(deltaT, delta_mu, current, 'ShowText','on');
% saveas(fig, 'current.pdf')
%==================================================================================
% x = linspace(-1000, 1000, 10000);
% y = linspace(0, omegac, 10000);
% plot3(x, y, my_integrant2(x, y));
% xlabel('E')
% ylabel('omega')
%==================================================================================
%functions

% % integrant
% function out = calc_current(E_lower, E_upper)
%    out = integral2(integrant2, E_lower, E_upper, 0, omegac);
% end

% integrant
function out = integrant1(E, omega)
   out = rho(omega) .* (N_R(omega) - N_L(omega)) .* (f_L_up(E) - f_L_down(E+omega)) .* A(E, omega);
end

% matrix A
function out = A(E, omega)
    Gamma_R = @(omega) rho(omega).*2.*pi;
    out = DL_up(E).*DL_down(E+omega).*Gamma_R(omega);
end

function out = rho(omega)
    % density of states for right MI, functoin.
    global omegac alpha
    out = 0.5*alpha.*omega.*exp(-1.*omega./omegac);
end

function out = DL_up(E)
    global E0_up
    out = 1./(E-E0_up+1i.*Gamma_L(E)./2) .* Gamma_L(E) .* 1./(E-E0_up-1i.*Gamma_L(E)./2);
end

function out = DL_down(E)
    global E0_down
    out = 1./(E-E0_down+1i.*Gamma_L(E)./2) .* Gamma_L(E) .* 1./(E-E0_down-1i.*Gamma_L(E)./2);
end

function out = Gamma_L(E)
    global Gamma0
    out = Gamma0;
end
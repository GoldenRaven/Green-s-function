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

global omegac
omegac = 80;

T0 = 300;                               % unit: K
DeltaT = 50;
T_L = T0 + DeltaT./2.0;               % left lead temperature, unit: K
T_R = T0 - DeltaT./2.0;              % right lead temperature, unit: K
beta_L = 1./(k_B.*T_L/meV);    % beta of left metal lead, unit: meV^-1
beta_R = 1./(k_B.*T_R/meV);      % beta of right MI lead, unit: meV^-1

mu0 = 20;
delta_mu = 10;% difference of spin voltage bias, unit: meV
mu_up = mu0 + delta_mu./2.0;
mu_down = mu0 - delta_mu./2.0;

% Bosonic distribution
N_L = @(omega) 1./(exp(beta_L.*(omega + delta_mu)) - 1);
N_R = @(omega) 1./(exp(beta_R.*(omega + delta_mu)) - 1);

% Fermionic distribution
f_L_up = @(E) 1./(exp(beta_L.*(E-mu_up))+1);
f_L_down = @(E) 1./(exp(beta_L.*(E-mu_down))+1);

% define integrant
my_integrant2 = @(E, omega) rho(omega) .* (N_R(omega) - N_L(omega)) .* (f_L_up(E) - f_L_down(E+omega)) .* A(E, omega);

% integral limits of E
E_limit = 5e2;
E_lower = -1.*E_limit;
E_upper = E_limit;

%==================================================================================
% % output
% display('Warning! the integral limit is [-2k_B*T0, 2k_B*T0], [0, 80]');
% display('');
% out = quad2d(my_integrant2, E_lower, E_upper, 0, omegac);
% % out = integral2(my_integrant2, E_lower, E_upper, 0, omegac);
% display(out);
%==================================================================================
dT = linspace(-2*T0, 2*T0, 100);
d_mu = linspace(-2*mu0, 2*mu0, 100);
data = zeros(100, 3);
count = 1
for deltaT = dT
    for delta_mu = d_mu
        out = quad2d(my_integrant2, E_lower, E_upper, 0, omegac);
        data(count, :) = [deltaT, delta_mu, out];
        count = count + 1;
        %display(out);
    end
end
display(count);
writematrix(data, 'current.txt')
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
    E0_up = 3;
    E0_down = -2;
    
    Gamma_R = rho(omega).*2.*pi;
    out = DL_up(E, E0_up).*DL_down(E+omega, E0_down).*Gamma_R;
end

function out = rho(omega)
    % density of states for right MI, functoin.
    global omegac
    alpha = 0.2;                             % dissipation strength, dimensionless
    out = 0.5.*alpha.*omega.*exp(-1.*omega/omegac);
end

function out = DL_up(E, e0_up)
    out = 1./(E-e0_up+1i.*Gamma_L(E)./2) .* Gamma_L(E) .* 1./(E-e0_up-1i.*Gamma_L(E)./2);
end

function out = DL_down(E, e0_down)
    out = 1./(E-e0_down+1i.*Gamma_L(E)./2) .* Gamma_L(E) .* 1./(E-e0_down-1i.*Gamma_L(E)./2);
end

function out = Gamma_L(E)
    Gamma0 = 0.1;
    out = Gamma0;
end
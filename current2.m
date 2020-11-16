% spin current for NM-QD-MI system
% 2020-11-04
% by ligy

%clc;
%close all;
%clear all;
%tic;

% constants
k_B = physconst('Boltzman');
charge_e = 1.602176634e-19;
meV = 1.0e-3.*charge_e;

global omegac alpha E0_up E0_down Gamma0 W delta_mu beta_L beta_R mu_up mu_down

omegac = 80;                            % cutoff frequency, unit: meV
alpha = 0.2;                            % disspation strength, dimensionless
T0 = 300;                               % average temperature, unit:K
mu0 = -10;                         % average spin baias, unit: meV
E0_up = -5;                              % QD up level, unit: meV
E0_down = 5;                           % QD down level, unit: meV
Gamma0=1;                               % effective coupling, unit:meV
W=80;                                   % bandwidth of left metal lead Lorentz spectral

%==================================================================================
% % output
% display('Warning! the integral limit is [-2k_B*T0, 2k_B*T0], [0, 80]');
%==================================================================================
%calculate current
fileID = fopen('current.txt','w');
%fileID2 = fopen('heat.txt','w');
% d_mu = 0;
% dT = 0;
dT = linspace(-1.999*T0, 1.999*T0, 30);
d_mu = linspace(-80, 80, 30);

for deltaT = dT
    for delta_mu = d_mu                 % difference of spin voltage bias, unit: meV
        T_L = T0 + deltaT./2.0;               % left lead temperature, unit: K
        T_R = T0 - deltaT./2.0;               % left lead temperature, unit: K
        beta_L = 1./(k_B.*T_L./meV);    % beta of left metal lead, unit: meV^-1
        beta_R = 1./(k_B.*T_R./meV);      % beta of right MI lead, unit: meV^-1
        mu_up = mu0 - delta_mu./2.0;    % spin-up chemical, unit: meV
        mu_down = mu0 + delta_mu./2.0;    % spin-down chemical, unit: meV
 
        % define integrant
        current_integrant2 = @(E, omega) rho(omega) .* w(E, omega) .* A(E, omega);
        % current_integrant3 = @(E, omega) rho(omega) .* (N_R(omega) - N_L(omega)) .* (f_L_up(E) - f_L_down(E+omega)) .* A(E, omega);

        % heat_integrant2 = @(E, omega) omega.*rho(omega) .* (N_R(omega) - N_L(omega)) .* (f_L_up(E) - f_L_down(E+omega)) .* A(E, omega);

        currt = quad2d(current_integrant2, -100, 100, 0, 100, 'Singular', true, 'MaxFunEvals', 10000);
        % currt2 = quad2d(current_integrant3, -100, 100, 0, 100, 'Singular', true, 'MaxFunEvals', 10000);
        %heat = quad2d(heat_integrant2, -100, 100, -100, 100, 'Singular', true, 'MaxFunEvals', 10000);
        fprintf(fileID, '%-15.10g%-15.10g%-15.10g\n',deltaT, delta_mu, -currt);
        %fprintf(fileID2, '%-15.10g%-15.10g%-15.10g\n',deltaT, delta_mu, heat);
        %currt

        %==================================================================================
        % % plot the integrant function
        % x = linspace(-200, 200, 10000);
        % y = linspace(-3*omegac, 3*omegac, 10000);
        % plot3(x, y, current_integrant2(x, y));
        % xlabel('E');
        % ylabel('omega');
        % text(100, 100, 0.04, 'dT=-500,d\_mu=-40,E0\_up=30,E0\_down=35,Gamma0=4,mu0=32.5,\{-200,200\}');
        %==================================================================================

    end
end
%toc;
%display('Time used: ', num2str(toc));
%==================================================================================
% % plot contour, save to file
% fig = figure;
% axes;
% set(fig, 'InvertHardcopy', 'off');
% contour(deltaT, delta_mu, current, 'ShowText','on');
% saveas(fig, 'current.pdf')
%==================================================================================
% x = linspace(-500, 500, 10000);
% y = linspace(-10*omegac, 10*omegac, 10000);
% plot3(x, y, current_integrant2(x, y));
% xlabel('E')
% ylabel('omega')
% text(100, 100, 0.04, 'dT=-500,d\_mu=-40,E0\_up=30,E0\_down=35,Gamma0=4,mu0=32.5,\{-200,200\}');
%==================================================================================
%functions

% % integrant
% function out = calc_current(E_lower, E_upper)
%    out = integral2(integrant2, E_lower, E_upper, 0, omegac);
% end

% integrant
% function out = integrant1(E, omega)
%    out = rho(omega) .* (N_R(omega) - N_L(delta_mu)) .* (f_L_up(E) - f_L_down(E+omega)) .* A(E, omega);
% end

% matrix A
function out = A(E, omega)
    out = DL_up(E).*DL_down(E+omega);
    %out = 1.0;                          % Ren Jie's case
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
    global Gamma0 W
    out = Gamma0;
    %out = Gamma0.*W.*W./(2*pi*(E.*E+W.*W));
end

function out = w(E, omega)
    out = f_L_up(E).*(1-f_L_down(E+omega)).*N_R(omega) - f_L_down(E+omega).*(1-f_L_up(E)).*(1+N_R(omega));
end

function out = N_L(omega)
    global beta_L delta_mu
    out = 1./(exp(beta_L.*(omega-delta_mu)) - 1);
end

function out = N_R(omega)
    global beta_R
    out = 1./(exp(beta_R.*(omega)) - 1);
end

function out = f_L_up(E)
    global beta_L mu_up
    out = 1./(exp(beta_L.*(E-mu_up)) + 1);
end

function out = f_L_down(E)
    global beta_L mu_down
    out = 1./(exp(beta_L.*(E-mu_down)) + 1);
end

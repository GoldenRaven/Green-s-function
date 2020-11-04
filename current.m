% spin current for NM-QD-MI system
% 2020-11-04
% by ligy

% constants
k_B =                                   % unit:

alpha = 0.2;                             % dissipation strength, dimensionless
omegac = 80;                             % cut off frequency, unit: meV
betaL = 1./(k_B.*T_L);                                % temperature of left metal lead, unit: meV
betaR = 1./(k_B.*T_R);                                % temperature of right MI lead, unit: meV
mu_up =                                  % spin voltage bias, unit: meV
mu_down =                                  % spin voltage bias, unit: meV
delta_mu = mu_up - mu_down                  % difference of spin voltage bias, unit: meV

% density of states for right MI, functoin.
rho = @(omega) pi.*alpha.*omega.*exp(-1.*omega/omegac);

% Bosonic distribution
N_L = @(omega) 1./(exp(betaL.*(omega + delta_mu)) + 1)
N_R = @(omega) 1./(exp(betaR.*(omega + delta_mu)) + 1)

% Fermionic distribution
f_L_up = @(E) 1./(exp(beta_L.*(E-mu_up)+1))
f_L_down = @(E) 1./(exp(beta_L.*(E-mu_down)+1))

% matrix A

A(E, omega) = @()

out = rho(1).*rho(2);
display(out)
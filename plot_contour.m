% plot contour, save to file
data = importdata('current.txt');
deltaT = data(:,1);
delta_mu = data(:,2);
currt = data(:,3);

T_R = 300;                               % average temperature, unit:K
mu_down = 32.5;                         % average spin baias, unit: meV

dT = linspace(-T_R, T_R, 80);
d_mu = linspace(-mu_down, mu_down, 80);
[x, y] = meshgrid(dT, d_mu);
z = griddata(deltaT, delta_mu, currt, x, y);

fig = figure;
axes;
surf(x, y, z);
set(fig, 'InvertHardcopy', 'off');
contourf(x, y, z, 'ShowText','on');
xlabel('deltaT');
ylabel('delta mu');
colorbar()
export_fig current.pdf

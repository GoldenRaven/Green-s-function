% plot contour, save to file
data = importdata('current.txt');
deltaT = data(:,1);
delta_mu = data(:,2);
currt = data(:,3);

T0 = 300;                               % average temperature, unit:K
mu0 = 32.5;                               % average spin baias, unit: meV
dT = linspace(-1.99*T0, 1.99*T0, 50);
d_mu = linspace(-1.99*mu0, 1.99*mu0, 50);
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

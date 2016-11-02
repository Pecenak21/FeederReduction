%% PV Power fitting
load('MAY2_2013_CLEAR_PV.mat') % Load whatever PV data you have
target = siDeployment('PointLoma');

%% Smooth PV data
a = 0.01;
pfilt = filter(a, [1 a-1], p); % Low pass filter

figure; plot(t, p, 'linewidth', 2)
hold on; plot(t, pfilt, '--r')

%% Plot efficiency
pv.pos.longitude = -117.222; pv.pos.latitude = 32.88; pv.pos.altitude = 100;
csk = clearSkyIrradiance( pv.pos , t, target.tilt, target.azimuth );

SF = (max(p)/max(csk.gi));

p_mod = csk.gi*SF;
eff = p./p_mod;
p_mod(eff < 0) = [];
eff(eff < 0) = [];
p_mod(eff > 1.5) = [];
eff(eff > 1.5) = [];

figure; plot(p_mod, eff); hold on; plot([0 2.5e4], [1 1], '--k')

%% Fit
eff_coeff = polyfit(p_mod, eff, 4);
figure; plot(p_mod, eff); hold on; plot([0 2.5e4], [1 1], '--k'); plot(p_mod, polyval(eff_coeff, p_mod), '--r')

% Test
figure; plot(t, p, 'linewidth', 2, 'color', 'g')
hold on; plot(t, polyval(eff_coeff, csk.gi*SF).*csk.gi*SF, '--r', 'linewidth', 2)
legend('measured', 'modeled', 'location', 'northeast'); datetick
clearvars;
close all;

path="/mnt/sdb/yyamaguchi/psdata2matlab/simulation/rawsignal/bubblerand2d/data/solid_liquid_reflector_linear1800.mat";
p=load(path);
t=p.kgrid.t_array;
p=p.sensor_data.p;
sz = size(p);
plt_idx = round(sz(1)/2);
start_idx = round(sz(2)/2);
p = p(plt_idx,start_idx:end);
t = t(start_idx:end);

figure(1);
plot(t*1e6,p*1e-3);
xlim([53 58]);
% ylim([-20 20]);
xlabel('Time [us]');
ylabel('Pressure [kPa]');
fontsize(18,"points");
saveas(gcf,"/mnt/sdb/yyamaguchi/psdata2matlab/simulation/processed/bubblerand2d/sandbox.png");


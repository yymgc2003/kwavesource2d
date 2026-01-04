clearvars;
close all;

config=jsondecode(fileread("/home/user01/Document/yyamaguchi/documents/simulationB4/sandbox/config.json"));
save_data_path = config.save_data_path;
save_pic_path = config.save_pic_path;

Nx = config.grid.Nx;
Ny = config.grid.Ny;
t_end = config.simulation.t_end;



mat_struct = load(fullfile(save_data_path, 'sandbox/glass_lossyballnew_8.mat'));
p1 = mat_struct.p;
kgrid1 = mat_struct.kgrid;
disp(size(kgrid1.t_array));

mat_struct = load(fullfile(save_data_path, 'sandbox/glass_lossyballnew_16.mat'));
p2 = mat_struct.p;
kgrid2 = mat_struct.kgrid;
disp(size(kgrid2.t_array));

mat_struct = load(fullfile(save_data_path, 'sandbox/glass_lossyballnew_32.mat'));
p3 = mat_struct.p;
kgrid3 = mat_struct.kgrid;

mat_struct = load(fullfile(save_data_path, 'sandbox/glass_lossyballnew_64.mat'));
p4 = mat_struct.p;
kgrid4 = mat_struct.kgrid;

figure('Position', [0,0,600,200]);
plot(1e6*kgrid1.t_array, p1*1e-3,'Color','magenta');
hold on;
plot(1e6*kgrid2.t_array, p2*1e-3,'Color','b');
hold on;
plot(1e6*kgrid3.t_array, p3*1e-3,'Color','r');
hold on;
plot(1e6*kgrid4.t_array, p4*1e-3,'Color','black');
xlim([27,32]);
ylim([-100,100]);
xlabel("Time [μs]");
ylabel("Pressure [kPa]")
fontsize(16,"points");
% legend("linear", "nonlinear");
legend("ppw=8", "ppw=16", "ppw=32", "ppw=64");
saveas(gcf, fullfile(save_pic_path, ['sandbox/glass_lossyballnew_ppw.png']));

an1=angle(hilbert(p1));
an2=angle(hilbert(p2));
an3=angle(hilbert(p3));
an4=angle(hilbert(p4));
figure('Position', [0,0,600,200]);
plot(1e6*kgrid1.t_array, an1,'Color','magenta');
hold on;
plot(1e6*kgrid2.t_array, an2,'Color','b');
hold on;
plot(1e6*kgrid3.t_array, an3,'Color','r');
hold on;
plot(1e6*kgrid4.t_array, an4,'Color','black');
xlim([27,32]);
ylim([-4,4]);
xlabel("Time [μs]");
ylabel("Phase [rad]")
fontsize(16,"points");
% legend("linear", "nonlinear");
legend("ppw=8", "ppw=16", "ppw=32", "ppw=64");
saveas(gcf, fullfile(save_pic_path, ['sandbox/glass_lossyballnew_angle_ppw.png']));

dif1 = min([abs(an1-an4); abs(an1+an4)],[],1);
disp(size(dif1));
dif2 = min([abs(an2-an4); abs(an2+an4)],[],1);
dif3 = min([abs(an3-an4); abs(an3+an4)],[],1);
figure('Position', [0,0,600,200]);
plot(1e6*kgrid4.t_array,dif1(1,:)/2/pi,'Color','magenta');
hold on;
plot(1e6*kgrid4.t_array,dif2(1,:)/2/pi,'Color','b');
hold on;
plot(1e6*kgrid4.t_array,dif3(1,:)/2/pi,'Color','r');
yline(0.1, '--', 'Color', 'black');
xlim([27,32]);
ylim([0,1]);
xlabel("Time [μs]");
ylabel("Phase [-]")
fontsize(16,"points");
% legend("linear", "nonlinear");
legend("ppw=8", "ppw=16", "ppw=32");
saveas(gcf, fullfile(save_pic_path, ['sandbox/glass_lossyballnew_angle_dif_ppw.png']));

mat_struct = load(fullfile(save_data_path, 'sandbox/air_lossyballnew_4.mat'));
p1 = mat_struct.p;
kgrid1 = mat_struct.kgrid;

mat_struct = load(fullfile(save_data_path, 'sandbox/air_lossyballnew_8.mat'));
p2 = mat_struct.p;
kgrid2 = mat_struct.kgrid;

mat_struct = load(fullfile(save_data_path, 'sandbox/air_lossyballnew_16.mat'));
p3 = mat_struct.p;
kgrid3 = mat_struct.kgrid;

% mat_struct = load(fullfile(save_data_path, 'sandbox/air_lossyball_c1200_10.mat'));
% p4 = mat_struct.p;
% kgrid4 = mat_struct.kgrid;

figure('Position', [0,0,600,200]);
plot(1e6*kgrid1.t_array, p1*1e-3);
hold on;
fontsize(16,"points");
plot(1e6*kgrid2.t_array, p2*1e-3);
hold on;
plot(1e6*kgrid3.t_array, p3*1e-3);
% hold on;
% plot(1e6*kgrid4.t_array, p4*1e-3);
xlim([28.5,30.5]);
ylim([-100,100]);
xlabel("Time [μs]");
ylabel("Pressure [kPa]")
% legend("linear", "nonlinear");
legend("ppw=4", "ppw=8", "ppw=16");
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossyballnew_ppw.png']));

an1=angle(hilbert(p1));
an2=angle(hilbert(p2));
an3=angle(hilbert(p3));
figure('Position', [0,0,600,200]);
plot(1e6*kgrid1.t_array, an1);
hold on;
plot(1e6*kgrid2.t_array, an2);
hold on;
plot(1e6*kgrid3.t_array, an3);
xlim([28.5,30.5]);
ylim([-4,4]);
xlabel("Time [μs]");
ylabel("Phase [rad]")
% legend("linear", "nonlinear");
legend("ppw=4", "ppw=8", "ppw=16");
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossyballnew_angle_ppw.png']));

dif1 = min([abs(an1-an3); abs(an1+an3)],[],1);
disp(size(dif1));
dif2 = min([abs(an2-an3); abs(an2+an3)],[],1);
figure('Position', [0,0,600,200]);
plot(1e6*kgrid3.t_array,dif1(1,:)/2/pi);
hold on;
plot(1e6*kgrid3.t_array,dif2(1,:)/2/pi);
yline(0.1, '--', 'Color', 'black');
xlim([28.5,30.5]);
ylim([0,1]);
xlabel("Time [μs]");
ylabel("Phase [-]")
fontsize(16,"points");
% legend("linear", "nonlinear");
legend("ppw=4", "ppw=8", "ppw=16");
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossyballnew_angle_dif_ppw.png']));

% mat_struct = load(fullfile(save_data_path, 'sandbox/glass_lossyballnew_4.mat'));
% p1 = mat_struct.p;
% kgrid1 = mat_struct.kgrid;

% mat_struct = load(fullfile(save_data_path, 'sandbox/glass_lossyballnew_8.mat'));
% p2 = mat_struct.p;
% kgrid2 = mat_struct.kgrid;

% mat_struct = load(fullfile(save_data_path, 'sandbox/glass_lossyballnew_16.mat'));
% p3 = mat_struct.p;
% kgrid3 = mat_struct.kgrid;

% mat_struct = load(fullfile(save_data_path, 'sandbox/air_lossyball8.mat'));
% p4 = mat_struct.p;
% kgrid4 = mat_struct.kgrid;

% figure('Position', [0,0,600,300]);
% plot(1e6*kgrid1.t_array, p1*1e-3);
% hold on;
% plot(1e6*kgrid2.t_array, p2*1e-3);
% hold on;
% plot(1e6*kgrid3.t_array, p3*1e-3);
% % hold on;
% % plot(1e6*kgrid4.t_array, p4*1e-3);
% xlim([26,30]);
% ylim([-100,100]);
% xlabel("Time [μs]");
% ylabel("Pressure [kPa]");
% fontsize(16,"points");
% % legend("linear", "nonlinear");
% legend("ppw=4", "ppw=8","ppw=16");
% saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossyballnew_ppw.png']));

mat_struct = load(fullfile(save_data_path, 'sandbox/glass_lossystrong1024_cfl0.05.mat'));
sensor_data = mat_struct.sensor_data;
kgrid1 = mat_struct.kgrid;

p1 = sensor_data.p;
p1 = mean(p1);

mat_struct = load(fullfile(save_data_path, 'sandbox/glass_nonlinear1024_cfl0.05.mat'));
sensor_data = mat_struct.sensor_data;
kgrid2 = mat_struct.kgrid;

p2 = sensor_data.p;
p2 = mean(p2);

figure('Position', [0,0,600,300]);
plot(1e6*kgrid1.t_array, p1*1e-3, "Color", "b");
hold on;
fontsize(16,"points");
plot(1e6*kgrid2.t_array, p2*1e-3, "Color", "r");
xlim([40,50]);
ylim([-30,30]);
xlabel("Time [μs]");
ylabel("Pressure [kPa]")
% legend("linear", "nonlinear");
legend("lossy", "lossless");
saveas(gcf, fullfile(save_pic_path, ['sandbox/glass_lossy_nonlinear.png']));

p1=p1(round(size(p1,2)/2):end);
p2=p2(round(size(p2,2)/2):end);
pf1 = fft(p1);
pf2 = fft(p2);
l1 = size(pf1,2);
l2 = size(pf2,2);
fs1 = 1/kgrid1.dt/l1*(0:l1-1);
fs2 = 1/kgrid2.dt/l2*(0:l2-1);

figure('Position', [0,0,600,300]);
plot(1e-6*fs1, abs(pf1)*2, "Color", "b");
hold on;
fontsize(16,"points");
plot(1e-6*fs2, abs(pf2), "Color", "r");
xlim([0,20]);
% ylim([-20,20]);
xlabel("Frequency [MHz]");
ylabel("|FFT(P)|")
% legend("linear", "nonlinear");
legend("lossy", "lossless");
saveas(gcf, fullfile(save_pic_path, ['sandbox/glass_lossy_nonlinear_fft.png']));
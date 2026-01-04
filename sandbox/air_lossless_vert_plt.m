clearvars;
close all;

config=jsondecode(fileread("/home/user01/Document/yyamaguchi/documents/simulationB4/sandbox/config.json"));
save_data_path = config.save_data_path;
save_pic_path = config.save_pic_path;

Nx = config.grid.Nx;
Ny = config.grid.Ny;
t_end = config.simulation.t_end;



mat_struct = load(fullfile(save_data_path, 'sandbox/air_lossy_vert_refined343_1_3_1536.mat'));
kgrid1 = mat_struct.kgrid;

p1 = mat_struct.p1;
t1 = mat_struct.p2;
% p1 = mean(p1);
% t1 = mean(t1);
A1 = max(p1(1:round(end*0.375)))-min(p1(1:round(end*0.375)));
A2 = max(p1(round(end*0.375):end))-min(p1(round(end*0.375):end));
raito1 = A2/A1

mat_struct = load(fullfile(save_data_path, 'sandbox/air_lossy_vert_refined1200_10_3_1536.mat'));
kgrid2 = mat_struct.kgrid;

p2 = mat_struct.p1;
t2 = mat_struct.p2;
A1 = max(p2(1:round(end*0.375)))-min(p2(1:round(end*0.375)));
A2 = max(p2(round(end*0.375):end))-min(p2(round(end*0.375):end));
raito2 = A2/A1

figure('Position', [0,0,600,300]);
plot(1e6*kgrid1.t_array, p1*1e-3, "Color", "b");
hold on;
plot(1e6*kgrid2.t_array, p2*1e-3, "Color", "r");
fontsize(20,"points");
xlim([10,25]);
ylim([-400,400]);
xlabel("Time [μs]");
ylabel("Pressure [kPa]")
legend("Genuine properties", "Fake properties");
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossy_properties2.png']));

figure('Position', [0,0,600,300]);
plot(1e6*kgrid1.t_array, t1*1e-3, "Color", "b");
hold on;
plot(1e6*kgrid2.t_array, t2*1e-3, "Color", "r");
fontsize(20,"points");
xlim([10,25]);
ylim([-100,100]);
xlabel("Time [μs]");
ylabel("Pressure [kPa]")
legend("Genuine properties", "Fake properties");
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossy_trans_properties2.png']));

% figure('Position', [0,0,600,300]);
% plot(1e6*kgrid1t.t_array, t1*1e-3, "Color", "b");
% hold on;
% fontsize(20,"points");
% plot(1e6*kgrid2t.t_array, t2*1e-3, "Color", "r");
% hold on;
% % plot(1e6*kgrid3.t_array, p3*1e-3, "Color", "g");
% % hold on;
% % plot(1e6*kgrid4.t_array, p4*1e-3);
% xlim([10,25]);
% ylim([-40,40]);
% xlabel("Time [μs]");
% ylabel("Pressure [kPa]")
% % legend("linear", "nonlinear");
% legend("Genuine properties", "Fake properties");
% saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossy_properties2_trans.png']));



% mat_struct = load(fullfile(save_data_path, 'sandbox/air_lossy343_1_cfl0.03.mat'));
% sensor_data = mat_struct.sensor_data;
% kgrid1 = mat_struct.kgrid;

% p1 = sensor_data.p;
% p1 = mean(p1);
% A1 = max(p1(round(end*0.375):end))-min(p1(round(end*0.375):end));
% A1 =A1/2

% mat_struct = load(fullfile(save_data_path, 'sandbox/air_lossy600_10_cfl0.03.mat'));
% sensor_data = mat_struct.sensor_data;
% kgrid2 = mat_struct.kgrid;

% p2 = sensor_data.p;
% p2 = mean(p2);
% A2 = max(p2(round(end*0.375):end))-min(p2(round(end*0.375):end));
% A2 = A2/2

% figure('Position', [0,0,600,300]);
% plot(1e6*kgrid1.t_array, p1*1e-3, "Color", "b");
% hold on;
% fontsize(20,"points");
% plot(1e6*kgrid2.t_array, p2*1e-3, "Color", "r");
% hold on;
% % plot(1e6*kgrid3.t_array, p3*1e-3, "Color", "g");
% % hold on;
% % plot(1e6*kgrid4.t_array, p4*1e-3);
% xlim([25,35]);
% ylim([-100,100]);
% xlabel("Time [μs]");
% ylabel("Pressure [kPa]")
% % legend("linear", "nonlinear");
% legend("Genuine properties", "Fake properties");
% saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossy_properties2.png']));
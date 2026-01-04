clearvars;
close all;

config=jsondecode(fileread("/home/user01/Document/yyamaguchi/documents/simulationB4/bubble_rand/config.json"));
save_dir = config.save_path;
[filepath, filename, ext] = fileparts(mfilename('fullpath'));

dataset_name = strsplit(filepath, '/');
dataset_name = string(dataset_name(end));
% dataset_name = glass_rand
save_dir = fullfile(save_dir, dataset_name);

save_pic_path=fullfile(save_dir,"pic");

t_end = config.simulation.t_end;



mat_struct = load(fullfile(save_dir, 'case64/data/sim_data1.mat'));
kgrid1 = mat_struct.kgrid;

p1 = mat_struct.p1;
p1 = mean(p1);

mat_struct = load(fullfile(save_dir, 'case64/data/sim_data2.mat'));
kgrid2 = mat_struct.kgrid;

p2 = mat_struct.p1;
p2 = mean(p2);

mat_struct = load(fullfile(save_dir, 'case64/data/sim_data3.mat'));
kgrid3 = mat_struct.kgrid;

p3 = mat_struct.p1;
p3 = mean(p3);

figure('Position', [0,0,600,300]);
plot(1e6*kgrid1.t_array, p1*1e-3, "Color", "b");
hold on;
plot(1e6*kgrid2.t_array, p2*1e-3, "Color", "r");
% hold on;
% plot(1e6*kgrid3.t_array, p3*1e-3);
fontsize(20,"points");
xlim([20,50]);
ylim([-100,100]);
xlabel("Time [μs]");
ylabel("Pressure [kPa]")
legend("dx=0.02 mm", "dx=0.04 mm");
saveas(gcf, fullfile(save_pic_path, 'dx_c343_rho10.png'));
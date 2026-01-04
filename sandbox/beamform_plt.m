clearvars;
close all;

config=jsondecode(fileread("/home/user01/Document/yyamaguchi/documents/simulationB4/sandbox/config.json"));
save_data_path = config.save_data_path;
save_pic_path = config.save_pic_path;

Nx = config.grid.Nx;
Ny = config.grid.Ny;
dx = config.grid.dx;

sensor_data = load(fullfile(save_data_path, 'sandbox/beamform2048.mat'));
sensor_data = sensor_data.sensor_data;
p_matrix = reshape(sensor_data.p_max, Nx, Ny);
p_matrix = transpose(p_matrix);
figure('Position', [0, 0, 720, 240]);
imagesc(dx:dx*Nx, dx:dx*Ny, p_matrix);
colorbar;
daspect([1,1,1]);
saveas(gcf, fullfile(save_pic_path, 'sandbox/beamform2048.png'));

figure;
p_axis = p_matrix(Ny/2,:);
plot(21:Nx-21, p_axis(21:Nx-21));
saveas(gcf, fullfile(save_pic_path, 'sandbox/beamform_axis2048.png'));
% Read CSV file as table (with header support)
config = jsondecode(fileread('../config.json'));
sample_table = readtable(fullfile(config.location_seedfiles_path, 'location2.csv'));
x = sample_table.Var1;
y = sample_table.Var2;
z = sample_table.Var3;
% 半径rの球を各点(x, y, z)を中心に描画する
r = 0.1; % 半径を適宜設定してください

% 球面の座標を生成
[phi, theta] = meshgrid(linspace(0, pi, 20), linspace(0, 2*pi, 40));
xs = sin(phi) .* cos(theta);
ys = sin(phi) .* sin(theta);
zs = cos(phi);


%saveas(gcf, fullfile(config.location_seedfiles_path, 'sample_spheres.png'));

% 2D scatter plot in XY plane with unit circle
figure;
scatter(x, y, 'filled');
hold on;
theta = linspace(0, 2*pi, 200);
plot(cos(theta), sin(theta), 'r-', 'LineWidth', 2);
hold off;
xlabel('X');
ylabel('Y');
title('Sample points and unit circle in XY plane');
grid on;
axis equal;
saveas(gcf, fullfile(config.location_seedfiles_path, 'pipeplot.png'));

% Overlay a unit cylinder on the previous 3D scatter plot
% (Assumes hold is already on from previous plotting)
n_cylinder = 100;
[theta_cyl, z_cyl] = meshgrid(linspace(0, 2*pi, n_cylinder), linspace(0, 1, n_cylinder));
x_cyl = cos(theta_cyl);
y_cyl = sin(theta_cyl);
% Draw the surface of the unit cylinder
surf(x_cyl, y_cyl, z_cyl, 'FaceAlpha', 0.8, 'EdgeColor', 'none', 'FaceColor', [0.8 0.2 0.2]);
% Draw the top and bottom edges of the cylinder
plot3(cos(theta_cyl(1,:)), sin(theta_cyl(1,:)), ones(1, n_cylinder), 'r-', 'LineWidth', 1.5);   % Top edge
plot3(cos(theta_cyl(1,:)), sin(theta_cyl(1,:)), zeros(1, n_cylinder), 'r-', 'LineWidth', 1.5);  % Bottom edge
% No need to call hold off or create a new figure, as this overlays on the existing plot
% Add legend or annotation if needed

% XZ plane projection (z from 0 to 1, square from 0 to 1)
figure;
scatter(x, z, 'filled');
hold on;
x_square = [-1 1 1 -1 -1];
z_square = [0 0 1 1 0];
plot(x_square, z_square, 'r-', 'LineWidth', 2);
hold off;
xlabel('X');
ylabel('Z');
title('Sample points and [0,1] square in XZ plane');
grid on;
axis equal;
saveas(gcf, fullfile(config.location_seedfiles_path, 'pipeplot_xz.png'));

% YZ plane projection (z from 0 to 1, square from 0 to 1)
figure;
scatter(y, z, 'filled');
hold on;
y_square = [-1 1 1 -1 -1];
z_square = [0 0 1 1 0];
plot(y_square, z_square, 'r-', 'LineWidth', 2);
hold off;
xlabel('Y');
ylabel('Z');
title('Sample points and [0,1] square in YZ plane');
grid on;
axis equal;
saveas(gcf, fullfile(config.location_seedfiles_path, 'pipeplot_yz.png'));
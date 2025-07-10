% CSVファイルをテーブルとして読み込む（ヘッダー対応）
config = jsondecode(fileread('../config.json'));
sample_table = readtable(fullfile(config.save_path, 'sample.csv'));

% Plot 2D scatter of X and Y columns and save as image
figure;
scatter(sample_table.X, sample_table.Y, 'filled');
hold on;
% Plot a circle of radius 1 centered at (0,0)
theta = linspace(0, 2*pi, 200);
plot(cos(theta), sin(theta), 'r-', 'LineWidth', 2);
hold off;
xlabel('X');
ylabel('Y');
title('Sample points and unit circle in XY plane');
grid on;
axis equal;
saveas(gcf, fullfile(config.save_path, 'pipeplot.png'));
% 3D scatter plot of sample points and a cylinder (centered at 0, radius 1, height from -1 to 1)
figure;
scatter3(sample_table.X, sample_table.Y, sample_table.Z, 36, 'filled'); % 3D scatter of samples
hold on;

% Create cylinder data (centered at (0,0), radius 1, height from -1 to 1)
n_cylinder = 100; % Number of points for smoothness
[theta_cyl, z_cyl] = meshgrid(linspace(0, 2*pi, n_cylinder), linspace(-1, 1, n_cylinder));
x_cyl = cos(theta_cyl);
y_cyl = sin(theta_cyl);

% Plot the cylinder surface
surf(x_cyl, y_cyl, z_cyl, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', [0.8 0.2 0.2]);

% Plot the top and bottom circles of the cylinder
plot3(cos(theta_cyl(1,:)), sin(theta_cyl(1,:)), ones(1, n_cylinder), 'r-', 'LineWidth', 1.5);    % Top
plot3(cos(theta_cyl(1,:)), sin(theta_cyl(1,:)), -ones(1, n_cylinder), 'r-', 'LineWidth', 1.5);   % Bottom

hold off;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Sample points and unit cylinder in 3D');
grid on;
axis equal;
view(3);
saveas(gcf, fullfile(config.save_path, 'pipeplot3d.png'));

% Plot XZ plane projection (reference: XY plane plot)
figure;
scatter(sample_table.X, sample_table.Z, 'filled'); % XZ scatter plot
hold on;
% Draw the projection of the unit square (Y=0) in XZ plane
x_square = [-1 1 1 -1 -1];
z_square = [-1 -1 1 1 -1];
plot(x_square, z_square, 'r-', 'LineWidth', 2); % XZ square projection
hold off;
xlabel('X');
ylabel('Z');
title('Sample points and unit square in XZ plane');
grid on;
axis equal;
saveas(gcf, fullfile(config.save_path, 'pipeplot_xz.png'));

% Plot YZ plane projection (reference: XY plane plot)
figure;
scatter(sample_table.Y, sample_table.Z, 'filled'); % YZ scatter plot
hold on;
% Draw the projection of the unit square (X=0) in YZ plane
y_square = [-1 1 1 -1 -1];
z_square = [-1 -1 1 1 -1];
plot(y_square, z_square, 'r-', 'LineWidth', 2); % YZ square projection
hold off;
xlabel('Y');
ylabel('Z');
title('Sample points and unit circle in YZ plane');
grid on;
axis equal;
saveas(gcf, fullfile(config.save_path, 'pipeplot_yz.png'));
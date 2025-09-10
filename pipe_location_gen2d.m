

% --- Script to generate multiple location CSV files ---
addpath('..');
config = jsondecode(fileread('config2d.json'));
num_repeat = config.simulation.num_dataset; % Number of times to repeat
m = config.simulation.num_particles;         % Number of samples per file
d = config.simulation.distance_particles;
for i = 1:num_repeat
    samples = glass_location_gen(m);
    config_file = 'config2d.json';
    if ~exist(config_file, 'file')
        error('Configuration file not found: %s', config_file);
    end
    config = jsondecode(fileread(config_file));
    save_path = config.location_seedfiles_path;
    if ~exist(save_path, 'dir')
        mkdir(save_path);
    end
    filename = fullfile(save_path, sprintf('location%d.csv', i));
    writematrix(samples', filename); % Save as CSV (transpose to get m rows)
end
% XY plane scatter plot with unit circle
figure;
scatter(samples(1,:), samples(2,:), 'filled');
hold on;
theta = linspace(0, 2*pi, 200);
plot(cos(theta), sin(theta), 'r-', 'LineWidth', 2);
hold off;
xlabel('X');
ylabel('Y');
title('Sample points and unit circle in XY plane');
grid on;
axis equal;
saveas(gcf, fullfile(save_path, 'pipeplot.png'));

% XZ plane projection (z from 0 to 1, square from -1 to 1 in X)
figure;
scatter(samples(1,:), samples(3,:), 'filled');
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
saveas(gcf, fullfile(save_path, 'pipeplot_xz.png'));

% YZ plane projection (z from 0 to 1, square from -1 to 1 in Y)
figure;
scatter(samples(2,:), samples(3,:), 'filled');
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
saveas(gcf, fullfile(save_path, 'pipeplot_yz.png'));

% 3D scatter plot with unit cylinder overlay
figure;
scatter3(samples(1,:), samples(2,:), samples(3,:), 36, 'filled');
hold on;
n_cylinder = 100;
[theta_cyl, z_cyl] = meshgrid(linspace(0, 2*pi, n_cylinder), linspace(0, 1, n_cylinder));
x_cyl = cos(theta_cyl);
y_cyl = sin(theta_cyl);
% Draw the surface of the unit cylinder
surf(x_cyl, y_cyl, z_cyl, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', [0.8 0.2 0.2]);
% Draw the top and bottom edges of the cylinder
plot3(cos(theta_cyl(1,:)), sin(theta_cyl(1,:)), ones(1, n_cylinder), 'r-', 'LineWidth', 1.5);   % Top edge
plot3(cos(theta_cyl(1,:)), sin(theta_cyl(1,:)), zeros(1, n_cylinder), 'r-', 'LineWidth', 1.5);  % Bottom edge
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D sample points and unit cylinder');
grid on;
axis equal;
hold off;
saveas(gcf, fullfile(save_path, 'pipeplot_3d.png'));
fprintf('plot saved to %s\n', save_path);
% --- 関数定義部分（同じファイルの末尾、または別ファイルに分離）---
function samples = glass_location_gen(m)
    % SAMPLING - Sample from 3D Gaussian distribution
    % 
    % Inputs:
    %   m - Number of samples to generate
    %
    % Outputs:
    %   samples - 3xm array containing the sampled points
    
    % Read configuration file
    config_file = 'config.json';
    if ~exist(config_file, 'file')
        error('Configuration file not found: %s', config_file);
    end
    
    % Read JSON configuration
    fid = fopen(config_file, 'r');
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);
    config = jsondecode(str);
    
    % Extract save path from configuration
    if ~isfield(config, 'save_path')
        error('save_path not found in configuration file');
    end
    save_path = config.save_path;
    
    % Create save directory if it doesn't exist
    if ~exist(save_path, 'dir')
        mkdir(save_path);
    end
    
    % Sample from 3D Gaussian (mean 0, var 1), keep only those inside unit circle in XY plane,
    % and ensure each sample is at least a certain distance away from all previous samples

    mu = [0, 0];
    sigma = eye(2);

    min_dist = config.simulation.distance_particles; % Minimum allowed Euclidean distance between any two samples (change as needed)
    samples = zeros(2, m); % Storage for samples
    count = 0;
    max_attempts = 100000; % Prevent infinite loop
    attempts = 0;
    while count < m && attempts < max_attempts
        % x, yはガウス分布、zは[-1,1]の一様分布からサンプリング
        xy = mvnrnd([0, 0], eye(2), 1)'; % 2x1ベクトル
        z = (1 - min_dist) * rand(1,1) + min_dist/2;   % zを[min_dist/2, 1-min_dist/2]の範囲で一様分布から生成
        candidate = [xy; z];             % 3x1 vec
        % Check if (X,Y) is inside unit circle
        if (candidate(1))^2 + (candidate(2))^2 <= (1-min_dist)^2
            % If this is the first sample, always accept
            if count == 0
                count = count + 1;
                samples(:, count) = candidate;
            else
                % Compute Euclidean distances to all previous samples
                dists = sqrt(sum((samples(:,1:count) - candidate).^2, 1));
                % Accept only if all distances are greater than or equal to min_dist
                if all(dists >= min_dist)
                    count = count + 1;
                    samples(:, count) = candidate;
                end
            end
        end
        attempts = attempts + 1;
    end
    if count == m
        % Calculate the Euclidean distance between all pairs of samples and display the minimum value
        D = pdist(samples'); % pdist computes the pairwise distances between rows of the matrix
        min_dist_val = min(D);
        fprintf('Minimum Euclidean distance between samples: %.6f\n', min_dist_val);
        disp(['Successfully generated ' num2str(m) ' samples! ']);
    end
    if count < m
        error('Could not generate enough samples in %d attempts. Try reducing min_dist or increasing max_attempts.', max_attempts);
    end
    
    % Save samples to CSV file
    csv_file = fullfile(save_path, 'sample.csv');
    sample_table = array2table(samples', 'VariableNames', {'X', 'Y', 'Z'});
    writetable(sample_table, csv_file);

    fprintf('Generated %d samples from 3D Gaussian distribution\n', m);
    fprintf('Samples saved to: %s\n', csv_file);
    fprintf('\nSample Statistics:\n');
    fprintf('Mean X: %.4f, Y: %.4f, Z: %.4f\n', mean(samples(1,:)), mean(samples(2,:)), mean(samples(3,:)));
    fprintf('Std  X: %.4f, Y: %.4f, Z: %.4f\n', std(samples(1,:)), std(samples(2,:)), std(samples(3,:)));

    % --- Plotting section: replicate the plots from tutorials/sampleplot.m ---

    % % XY plane scatter plot with unit circle
    % figure;
    % scatter(samples(1,:), samples(2,:), 'filled');
    % hold on;
    % theta = linspace(0, 2*pi, 200);
    % plot(cos(theta), sin(theta), 'r-', 'LineWidth', 2);
    % hold off;
    % xlabel('X');
    % ylabel('Y');
    % title('Sample points and unit circle in XY plane');
    % grid on;
    % axis equal;
    % saveas(gcf, fullfile(save_path, 'pipeplot.png'));

    % % XZ plane projection (z from 0 to 1, square from -1 to 1 in X)
    % figure;
    % scatter(samples(1,:), samples(3,:), 'filled');
    % hold on;
    % x_square = [-1 1 1 -1 -1];
    % z_square = [0 0 1 1 0];
    % plot(x_square, z_square, 'r-', 'LineWidth', 2);
    % hold off;
    % xlabel('X');
    % ylabel('Z');
    % title('Sample points and [0,1] square in XZ plane');
    % grid on;
    % axis equal;
    % saveas(gcf, fullfile(save_path, 'pipeplot_xz.png'));

    % % YZ plane projection (z from 0 to 1, square from -1 to 1 in Y)
    % figure;
    % scatter(samples(2,:), samples(3,:), 'filled');
    % hold on;
    % y_square = [-1 1 1 -1 -1];
    % z_square = [0 0 1 1 0];
    % plot(y_square, z_square, 'r-', 'LineWidth', 2);
    % hold off;
    % xlabel('Y');
    % ylabel('Z');
    % title('Sample points and [0,1] square in YZ plane');
    % grid on;
    % axis equal;
    % saveas(gcf, fullfile(save_path, 'pipeplot_yz.png'));

    % % 3D scatter plot with unit cylinder overlay
    % figure;
    % scatter3(samples(1,:), samples(2,:), samples(3,:), 36, 'filled');
    % hold on;
    % n_cylinder = 100;
    % [theta_cyl, z_cyl] = meshgrid(linspace(0, 2*pi, n_cylinder), linspace(0, 1, n_cylinder));
    % x_cyl = cos(theta_cyl);
    % y_cyl = sin(theta_cyl);
    % % Draw the surface of the unit cylinder
    % surf(x_cyl, y_cyl, z_cyl, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', [0.8 0.2 0.2]);
    % % Draw the top and bottom edges of the cylinder
    % plot3(cos(theta_cyl(1,:)), sin(theta_cyl(1,:)), ones(1, n_cylinder), 'r-', 'LineWidth', 1.5);   % Top edge
    % plot3(cos(theta_cyl(1,:)), sin(theta_cyl(1,:)), zeros(1, n_cylinder), 'r-', 'LineWidth', 1.5);  % Bottom edge
    % xlabel('X');
    % ylabel('Y');
    % zlabel('Z');
    % title('3D sample points and unit cylinder');
    % grid on;
    % axis equal;
    % hold off;
    % saveas(gcf, fullfile(save_path, 'pipeplot_3d.png'));
    % fprintf('plot saved to %s\n', save_path);
    % If you want to run the plotting script directly (optional, if available and path is correct):
    % try
    %     run(fullfile(fileparts(mfilename('fullpath')), 'tutorials', 'sampleplot.m'));
    % catch ME
    %     warning('Could not run sampleplot.m: %s', ME.message);
    % end

end
% --- 英語での説明 ---
% In MATLAB, you cannot define a function inside a script file and execute the script as a function.
% To fix this, separate the script part (variable assignment and function call) from the function definition.
% Place the function definition at the end of the file or in a separate file.
% Then, run the script part to execute your code.

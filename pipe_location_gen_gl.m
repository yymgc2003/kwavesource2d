% --- Script to generate multiple location CSV files ---
addpath('..');
config = jsondecode(fileread('config_gl.json'));
num_repeat = config.simulation.num_dataset; % Number of times to repeat
for i = 1:num_repeat
    samples = gas_location_gen();
    config_file = 'config_gl.json';
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

% --- 関数定義部分（同じファイルの末尾、または別ファイルに分離）---
function samples = gas_location_gen()
    % SAMPLING - Sample from 3D Gaussian distribution
    % 
    % Inputs:
    %   m - Number of samples to generate
    %
    % Outputs:
    %   samples - 3xm array containing the sampled points
    
    % Read configuration file
    config_file = 'config_gl.json';
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

    num_gas = 0;
    cur_gas_fraction = 0;
    %max_gas_fraction = config.simulation.max_gas_fraction;
    %min_gas_fraction = config.simulation.min_gas_fraction;
    %target_gas_fraction = (max_gas_fraction - min_gas_fraction)*rand + min_gas_fraction;
    target_gas_fraction = config.simulation.target_gas_fraction;
    cutoff_gas_fraction = config.simulation.cutoff_gas_fraction;
    inner_radius = config.pipe.inner_radius;
    min_gas_radius = config.simulation.min_gas_radius;
    min_gas_radius = min_gas_radius/inner_radius;
    min_dist = config.simulation.distance_gas;
    min_dist = min_dist / inner_radius;
    attempts_radius = 0;
    max_attempts_radius = 1000;

    while attempts_radius < max_attempts_radius
        num_gas = 0;
        cur_gas_fraction = 0;
        max_gas_radius = sqrt(target_gas_fraction);
        gas_radius_arr = zeros(1000);
        
        % Generate radius first.
        if target_gas_fraction > cutoff_gas_fraction
            num_gas = num_gas + 1;
            cur_min_gas_radius = 0.9*max_gas_radius;
            cur_gas_radius = (max_gas_radius - cur_min_gas_radius)*rand + cur_min_gas_radius;
            cur_gas_fraction = cur_gas_fraction + cur_gas_radius^2;
            gas_radius_arr(num_gas) = cur_gas_radius;
            max_gas_radius = max(0.5 - cur_gas_radius*0.5 - min_dist, min_gas_radius*2);
        end
        while cur_gas_fraction < target_gas_fraction
            num_gas = num_gas + 1;
            cur_gas_radius = (max_gas_radius - min_gas_radius)*rand + min_gas_radius;
            cur_gas_fraction = cur_gas_fraction + cur_gas_radius^2;
            gas_radius_arr(num_gas) = cur_gas_radius;
            max_gas_radius_candidate = min(0.5 - gas_radius_arr(1)*0.5 - min_dist, ...
                                        1.1*sqrt(target_gas_fraction - cur_gas_fraction));
            max_gas_radius = max(max_gas_radius_candidate, min_gas_radius*2);
        end
        gas_radius_arr = gas_radius_arr(1:num_gas);
        gas_radius_arr = sort(gas_radius_arr, 'descend');
        samples = zeros(3, num_gas);
        count = 0;
        max_attempts = 1000; % Prevent infinite loop
        attempts = 0;



        while count < num_gas && attempts < max_attempts
            % x, yはガウス分布、zは[-1,1]の一様分布からサンプリング
            %xy = mvnrnd([0, 0], eye(2), 1)'; % 2x1ベクトル
            xy = [2*rand() - 1; 2*rand() - 1];
            cur_gas_radius = gas_radius_arr(count+1);
            cur_min_dist = min_dist + cur_gas_radius;
            %z = (1 - cur_min_dist) * rand(1,1) + min_dist/2;   % zを[min_dist/2, 1-min_dist/2]の範囲で一様分布から生成
            candidate = [xy; cur_gas_radius];             % 3x1 vec
            % Check if (X,Y) is inside unit circle
                % If this is the first sample, always accept
            if count ~= 0
                if (candidate(1))^2 + (candidate(2))^2 <= (1-cur_min_dist)^2
                    % Compute Euclidean distances to all previous samples
                    dists = -samples(3, 1:count) + sqrt(sum((samples(1:2,1:count) - candidate(1:2)).^2, 1));
                    % Accept only if all distances are greater than or equal to min_dist
                    if all(dists >= cur_min_dist)
                        count = count + 1;
                        samples(:, count) = candidate;
                        attempts = 0;
                    end
                end
            else
                if (candidate(1))^2 + (candidate(2))^2 <= ((1-cur_gas_radius)/3)^2
                    count = count + 1;
                    samples(:, count) = candidate;
                end
            end
            attempts = attempts + 1;
        end
        attempts_radius = attempts_radius + 1;
        if count == num_gas
            % Calculate the Euclidean distance between all pairs of samples and display the minimum value
            D = pdist(samples'); % pdist computes the pairwise distances between rows of the matrix
            min_dist_val = min(D);
            fprintf('Minimum Euclidean distance between samples: %.6f\n', min_dist_val);
            disp(['Successfully generated ' num2str(num_gas) ' samples! ']);
            fprintf("Target gas fraction: %.4f\n", target_gas_fraction);
            fprintf("Current gas fraction: %.4f\n", cur_gas_fraction);
            break;
        end
        if count < num_gas
            attempts_radius = attempts_radius + 1;
        end
        if attempts_radius == max_attempts_radius
            error('Could not generate enough samples in %.4f fraction. Try reducing min_dist or increasing max_attempts.', target_gas_fraction);
        end
    end
    
    % Save samples to CSV file
    csv_file = fullfile(save_path, 'sample.csv');
    sample_table = array2table(samples', 'VariableNames', {'X', 'Y', 'Z'});
    writetable(sample_table, csv_file);

    fprintf('Generated %d samples from 3D Gaussian distribution\n', num_gas);
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
% --- Script to generate multiple location CSV files ---
addpath('..');
config = jsondecode(fileread('config3d.json'));
num_repeat = config.simulation.num_dataset; % Number of times to repeat
num_bubble = config.simulation.num_bubble;
slug_length = config.simulation.slug_length;
slug_range = config.simulation.slug_range;
for i = 1:num_repeat
    if config.simulation.flow_pattern == "bubble"
        samples = gas_location_gen3d(num_bubble);
    end
    if config.simulation.flow_pattern == "slug"
        samples = slug_location_gen3d(slug_length, slug_range);
    end
    config_file = 'config3d.json';
    if ~exist(config_file, 'file')
        error('Configuration file not found: %s', config_file);
    end
    config = jsondecode(fileread(config_file));
    save_path = config.location_seedfiles_path;
    if ~exist(save_path, 'dir')
        mkdir(save_path);
    end
    filename = fullfile(save_path, sprintf('location%d.csv', i));
    writematrix(transpose(samples), filename); % Save as CSV (transpose to get m rows)
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
function samples = gas_location_gen3d(num_bubble)
    % SAMPLING - Sample from 3D Gaussian distribution
    % 
    % Inputs:
    %   m - Number of samples to generate
    %
    % Outputs:
    %   samples - 3xm array containing the sampled points
    
    % Read configuration file
    config_file = 'config3d.json';
    if ~exist(config_file, 'file')
        error('Configuration file not found: %s', config_file);
    end
    
    % Read JSON configuration
    fid = fopen(config_file, 'r');
    raw = fread(fid, inf);
    str = char(transpose(raw));
    fclose(fid);
    config = jsondecode(str);
    
    % Extract save path from configuration
    if ~isfield(config, 'save_path')
        error('save_path not found in configuration file');
    end
    save_path = config.save_path;
    
    % Create save directory if it doesnt exist
    %if ~exist(save_path, 'dir')
    %    mkdir(save_path);
    %end
    
    % Sample from 3D Gaussian (mean 0, var 1), keep only those inside unit circle in XY plane,
    % and ensure each sample is at least a certain distance away from all previous samples

    mu = [0, 0];
    sigma = eye(2);

    k_gamma=5;
    theta_gamma=0.5;

    cur_gas_fraction = 0;
    inner_radius = config.pipe.inner_radius;
    min_diameter_bubble = config.simulation.min_diameter_bubble;
    min_diameter_bubble = min_diameter_bubble/inner_radius;
    min_dist = config.simulation.distance_bubble;
    min_dist = min_dist / inner_radius;
    z_range = config.grid.Nz*(config.grid.dz*1e3) / inner_radius;
    attempts_radius = 0;
    max_attempts_radius = 1000;
    diameter_bubble = zeros(num_bubble);
    pd = makedist('Gamma','a',k_gamma,'b',theta_gamma);
    samples = zeros(8,num_bubble); %4番目はバブル長径、5番目は短径、6~8はx,y,z周りの回転角


    while attempts_radius < max_attempts_radius
        count = 0;
        cur_gas_fraction = 0;
        diameter_bubble = zeros(num_bubble);
        
        while count < num_bubble
            cur_diameter_bubble = random(pd)/inner_radius;
            if cur_diameter_bubble < 0.4 && cur_diameter_bubble> min_diameter_bubble
                count = count+1;
                diameter_bubble(count)=cur_diameter_bubble;
            end
        end
        diameter_bubble = sort(diameter_bubble, 'descend');
        count = 0;
        max_attempts = 1000; % Prevent infinite loop
        attempts = 0;
        samples = zeros(8,count);



        while count < num_bubble && attempts < max_attempts
            cur_diameter_bubble = diameter_bubble(count+1);
            cur_min_dist = min_dist + cur_diameter_bubble/2;
            % x, yはガウス分布、zは[-1,1]の一様分布からサンプリング
            xy = transpose(mvnrnd([0, 0], eye(2)./(cur_diameter_bubble*inner_radius), 1)); % 2x1ベクトル
            %xy = [2*rand-1; 2*rand-1];
            z = (z_range - cur_diameter_bubble) * rand + cur_diameter_bubble/2;   % zを[min_dist/2, 1-min_dist/2]の範囲で一様分布から生成
            candidate = [xy; z];             % 3x1 vec
            % Check if (X,Y) is inside unit circle
                % If this is the first sample, always accept
            if count ~= 0
                if (candidate(1))^2 + (candidate(2))^2 <= (1-cur_min_dist)^2
                    % Compute Euclidean distances to all previous samples
                    dists = -samples(4, 1:count)./2 + sqrt(sum((samples(1:3,1:count) - candidate(1:3)).^2, 1));
                    % Accept only if all distances are greater than or equal to min_dist
                    if all(dists >= cur_min_dist)
                        count = count + 1;
                        samples(1:3, count) = candidate;
                        samples(4, count) = cur_diameter_bubble;
                        samples(5, count) = cur_diameter_bubble / (rand + 1);
                        samples(6:8, count) = [2*pi*rand, 2*pi*rand, 2*pi*rand];
                        cur_gas_fraction = cur_gas_fraction + samples(5, count)*cur_diameter_bubble^2/z_range/6;
                        attempts = 0;
                        %fprintf('Generated %d bubbles at %.4f %.4f %.4f diameter %.4f\n', count, samples(1,count),samples(2,count),samples(3,count), samples(4,count));
                    end
                end
            else
                if (candidate(1))^2 + (candidate(2))^2 <= (1-cur_min_dist)^2
                    count = count + 1;
                    samples(1:3, count) = candidate;
                    samples(4, count) = cur_diameter_bubble;
                    samples(5, count) = cur_diameter_bubble / (rand*0.5 + 1);
                    samples(6:8, count) = [2*pi*rand, 2*pi*rand, 2*pi*rand];
                end
            end
            attempts = attempts + 1;
        end
        if count == num_bubble
            % Calculate the Euclidean distance between all pairs of samples and display the minimum value
            D = pdist(transpose(samples(1:3,:))); % pdist computes the pairwise distances between rows of the matrix
            min_dist_val = min(D);
            fprintf('Minimum Euclidean distance between samples: %.6f\n', min_dist_val);
            disp(['Successfully generated ' num2str(num_bubble) ' samples! ']);
            fprintf("Current gas fraction: %.4f\n", cur_gas_fraction);
            break;
        end
        if count < num_bubble
            attempts_radius = attempts_radius + 1;
            %fprintf('Failed\n');
        end
        if attempts_radius == max_attempts_radius
            error('Could not generate enough samples in %.4f fraction. Try reducing min_dist or increasing max_attempts.', num_bubble);
        end
    end
    
    % Save samples to CSV file
    csv_file = fullfile(save_path, 'sample.csv');
    sample_table = array2table(transpose(samples), 'VariableNames', {'X', 'Y', 'Z', 'A', 'B', 'RL', 'PC', 'YW'});
    writetable(sample_table, csv_file);

    fprintf('Generated %d samples from 3D Gaussian distribution\n', num_bubble);
    fprintf('Samples saved to: %s\n', csv_file);
    fprintf('\nSample Statistics:\n');
    fprintf('Mean X: %.4f, Y: %.4f, Z: %.4f, R: %.4f\n', mean(samples(1,:)), mean(samples(2,:)), mean(samples(3,:)), mean(samples(4,:)));
    fprintf('Std  X: %.4f, Y: %.4f, Z: %.4f, R: %.4f\n', std(samples(1,:)), std(samples(2,:)), std(samples(3,:)), std(samples(4,:)));

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
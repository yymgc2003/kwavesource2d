% --- Script to generate multiple location CSV files ---
%anu(環状流)において、このファイルでは円がどこを通るかを指定する
%i=1では、必ず真の円になる
addpath('..');
config = jsondecode(fileread('config_anu.json'));
num_repeat = config.simulation.num_dataset; % Number of times to repeat
for i = 1:num_repeat
    samples = gas_location_gen(i);
    config_file = 'config_anu.json';
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

% --- 関数定義部分（同じファイルの末尾、または別ファイルに分離）---
function samples = gas_location_gen(i)
    % SAMPLING - Sample from 3D Gaussian distribution
    % 
    % Inputs:
    %   m - Number of samples to generate
    %
    % Outputs:
    %   samples - 3xm array containing the sampled points
    
    % Read configuration file
    config_file = 'config_anu.json';
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

    num_gas = 0;
    cur_gas_fraction = 0;
    %max_gas_fraction = config.simulation.max_gas_fraction;
    %min_gas_fraction = config.simulation.min_gas_fraction;
    %target_gas_fraction = (max_gas_fraction - min_gas_fraction)*rand + min_gas_fraction;
    inner_radius = config.pipe.inner_radius;
    min_dist = config.simulation.distance_gas;
    min_dist = min_dist / inner_radius;
    spline_point_num = config.simulation.annular_spline_point_num;
    mu = config.simulation.annular_radius_mean;
    mu = mu/inner_radius;
    sigma = config.simulation.annular_radius_std;
    sigma = sigma / inner_radius;
    count = 0;
    attempts = 0;
    max_attempts = 10000;
    samples = zeros(1, spline_point_num+1);
    while attempts < max_attempts
        samples_tmp = zeros(1, spline_point_num);
        while count < spline_point_num
            radius_candidate = randn;
            if abs(radius_candidate) < 2
                radius_candidate = radius_candidate*sigma + mu;
                count = count + 1;
                samples_tmp(1, count) = radius_candidate;
            end
        end
        samples_tmp(1, spline_point_num+1) = samples_tmp(1,1);
        theta = linspace(0,2*pi,spline_point_num+1);
        cs = spline(theta, [0 samples_tmp 0]);
        if max(ppval(cs,theta))<1-min_dist
            samples = samples_tmp;
            break;
        end
        attempts = attempts + 1;
    end
    if i==1
        samples(1, :) = mu;
    end

    % Save samples to CSV file
    csv_file = fullfile(save_path, 'sample.csv');
    sample_table = array2table(samples', 'VariableNames', {'X'});
    writetable(sample_table, csv_file);

    fprintf('Samples saved to: %s\n', csv_file);
    fprintf('\nSample Statistics:\n');
    fprintf('Mean X: %.4f, Std  X: %.4f\n', mean(samples(:)), std(samples(:)));
end
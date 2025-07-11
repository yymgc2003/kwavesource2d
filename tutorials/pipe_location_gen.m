% MATLABでは、スクリプトファイル（.mファイル）内で関数定義を行う場合、そのファイル全体を関数ファイルとして書き直す必要があります。
% しかし、スクリプトとして実行したい場合は、関数定義を別ファイルに分けるか、関数呼び出し部分と関数定義部分を分離してください。
% 例えば、以下のようにスクリプト部分と関数部分を分けて記述します。

% --- Script to generate multiple location CSV files ---
num_repeat = 5; % Number of times to repeat
m = 10;         % Number of samples per file

for i = 1:num_repeat
    samples = glass_location_gen(m);
    config_file = '../config.json';
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
function samples = glass_location_gen(m)
    % SAMPLING - Sample from 3D Gaussian distribution
    % 
    % Inputs:
    %   m - Number of samples to generate
    %
    % Outputs:
    %   samples - 3xm array containing the sampled points
    
    % Read configuration file
    config_file = '../config.json';
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

    mu = [0, 0, 0];
    sigma = eye(3);

    min_dist = 0.4; % Minimum allowed Euclidean distance between any two samples (change as needed)
    samples = zeros(3, m); % Storage for samples
    count = 0;
    max_attempts = 100000; % Prevent infinite loop
    attempts = 0;
    while count < m && attempts < max_attempts
        % x, yはガウス分布、zは[-1,1]の一様分布からサンプリング
        xy = mvnrnd([0, 0], eye(2), 1)'; % 2x1ベクトル
        z = -1 + 2 * rand(1,1);          % [-1,1]の一様分布
        candidate = [xy; z];             % 3x1ベクトル
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
config = jsondecode(fileread('../config.json'));
end
% --- 英語での説明 ---
% In MATLAB, you cannot define a function inside a script file and execute the script as a function.
% To fix this, separate the script part (variable assignment and function call) from the function definition.
% Place the function definition at the end of the file or in a separate file.
% Then, run the script part to execute your code.

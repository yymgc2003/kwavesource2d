% =========================================================================
% 2D k-Wave simulation with a focusing arc transducer and vinyl pipe
% =========================================================================
clearvars;
close all;
DATA_CAST='gpuArray-single';
% -------------------------------------------------------------------------
% 1) 設定ファイルの読み込み
% -------------------------------------------------------------------------
config = jsondecode(fileread('/home/user01/Document/yyamaguchi/documents/kwavesource2d/config.json'));

num_bubble = config.simulation.num_bubble;
mu = [0, 0];
sigma = eye(2);

k_gamma=5;
theta_gamma=0.5;

cur_gas_fraction = 0;
inner_radius = config.pipe.inner_radius;
min_diameter_bubble = config.simulation.min_diameter_bubble;
min_dist = config.simulation.distance_bubble;
min_dist_wall = config.simulation.min_dist_wall;
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
        cur_diameter_bubble = random(pd);
        if cur_diameter_bubble < 6 && cur_diameter_bubble> min_diameter_bubble
            count = count+1;
            diameter_bubble(count)=cur_diameter_bubble;
        end
    end
    diameter_bubble = sort(diameter_bubble, 'descend');
    count = 0;
    max_attempts = 1000; % Prevent infinite loop
    attempts = 0;
    samples = zeros(6,count);



    while count < num_bubble && attempts < max_attempts
        cur_diameter_bubble = diameter_bubble(count+1);
        cur_min_dist = min_dist + cur_diameter_bubble/2;
        cur_min_dist_wall = min_dist_wall + cur_diameter_bubble/2;
        % x, yはガウス分布、zは[-1,1]の一様分布からサンプリング
        if cur_diameter_bubble >= 5.5e-3
            xy = transpose(mvnrnd([0, 0], 0.1*eye(2), 1)); % 2x1ベクトル
        end
        if cur_diameter_bubble < 5.5e-3 & cur_diameter_bubble >= 3.5e-3
            R = normrnd(0.7, 0.2)*inner_radius;
            theta = 2*pi*rand;
            xy = [R*cos(theta); R*sin(theta)];
        end
        if cur_diameter_bubble < 3.5e-3
            R = normrnd(0.8, 0.1)*inner_radius;
            theta = 2*pi*rand;
            xy = [R*cos(theta); R*sin(theta)];
        end
        % xy = transpose(mvnrnd([0, 0], eye(2)./(cur_diameter_bubble), 1)); % 2x1ベクトル
        %xy = [2*rand-1; 2*rand-1];
        % z = (z_range - cur_diameter_bubble) * rand - z_range/2  + cur_diameter_bubble/2;   % zを[min_dist/2, 1-min_dist/2]の範囲で一様分布から生成
        candidate = xy;            % 3x1 vec
        % Check if (X,Y) is inside unit circle
            % If this is the first sample, always accept
        if count ~= 0
            if (candidate(1))^2 + (candidate(2))^2 <= (1-cur_min_dist_wall)^2
                % Compute Euclidean distances to all previous samples
                dists = -samples(3, 1:count)./2 + sqrt(sum((samples(1:2,1:count) - candidate(1:2)).^2, 1));
                % Accept only if all distances are greater than or equal to min_dist
                if all(dists >= cur_min_dist)
                    count = count + 1;
                    samples(1:2, count) = candidate;
                    samples(3, count) = cur_diameter_bubble;
                    samples(4, count) = cur_diameter_bubble/(rand+1);
                    % samples(4, count) = cur_diameter_bubble;
                    samples(5:6, count) = [2*pi*rand, 2*pi*rand];
                    cur_gas_fraction = cur_gas_fraction + samples(4, count)*cur_diameter_bubble^2;
                    attempts = 0;
                    %fprintf('Generated %d bubbles at %.4f %.4f %.4f diameter %.4f\n', count, samples(1,count),samples(2,count),samples(3,count), samples(4,count));
                end
            end
        else
            if (candidate(1))^2 + (candidate(2))^2 <= (1-cur_min_dist)^2
                count = count + 1;
                samples(1:2, count) = candidate;
                samples(3, count) = cur_diameter_bubble;
                % samples(4, count) = cur_diameter_bubble;
                samples(4, count) = cur_diameter_bubble/(rand+1);
                samples(5:6, count) = [2*pi*rand, 2*pi*rand];
            end
        end
        attempts = attempts + 1;
    end
    if count == num_bubble
        % Calculate the Euclidean distance between all pairs of samples and display the minimum value
        D = pdist(transpose(samples(1:2,:))); % pdist computes the pairwise distances between rows of the matrix
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
save_path = "/mnt/sdb/yyamaguchi/kwavesource2d/bubblerand/loc_seed";

csv_file = fullfile(save_path, 'sample.csv');
name = ["center x","center y","major axis length","minor axis length","raw","pitch"];
T = array2table(transpose(samples), 'VariableNames', name);
writetable(T, csv_file);

fprintf('Generated %d samples from 3D Gaussian distribution\n', num_bubble);
fprintf('Samples saved to: %s\n', csv_file);
fprintf('\nSample Statistics:\n');
fprintf('Mean X: %.4f, Y: %.4f, Z: %.4f, R: %.4f\n', mean(samples(1,:)), mean(samples(2,:)), mean(samples(3,:)), mean(samples(4,:)));
fprintf('Std  X: %.4f, Y: %.4f, Z: %.4f, R: %.4f\n', std(samples(1,:)), std(samples(2,:)), std(samples(3,:)), std(samples(4,:)));

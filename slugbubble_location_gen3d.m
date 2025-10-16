function samples_tuple = slugbubble_location_gen3d(num_bubble, slug_length, slug_range)
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
    min_dist_wall = config.simulation.min_dist_wall;
    min_dist_wall = min_dist_wall/ inner_radius;
    %z_range = config.grid.Nz*(config.grid.dz*1e3) / inner_radius;
    z_range =  2;
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
            cur_min_dist_wall = min_dist_wall + cur_diameter_bubble/2;
            % x, yはガウス分布、zは[-1,1]の一様分布からサンプリング
            xy = transpose(mvnrnd([0, 0], eye(2)./(cur_diameter_bubble*inner_radius), 1)); % 2x1ベクトル
            %xy = [2*rand-1; 2*rand-1];
            z = (z_range/2 - cur_diameter_bubble) * rand - z_range/2  + cur_diameter_bubble/2;   % zを[min_dist/2, 1-min_dist/2]の範囲で一様分布から生成
            candidate = [xy; z];             % 3x1 vec
            % Check if (X,Y) is inside unit circle
                % If this is the first sample, always accept
            if count ~= 0
                if (candidate(1))^2 + (candidate(2))^2 <= (1-cur_min_dist_wall)^2
                    % Compute Euclidean distances to all previous samples
                    dists = -samples(4, 1:count)./2 + sqrt(sum((samples(1:3,1:count) - candidate(1:3)).^2, 1));
                    % Accept only if all distances are greater than or equal to min_dist
                    if all(dists >= cur_min_dist)
                        count = count + 1;
                        samples(1:3, count) = candidate;
                        samples(4, count) = cur_diameter_bubble;
                        samples(5, count) = cur_diameter_bubble/(rand+1);
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
                    samples(5, count) = cur_diameter_bubble;
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
    name = ["center x","center y","center z","major axis length","minor axis length","raw","pitch","yaw"];
    T = array2table(transpose(samples), 'VariableNames', name);
    writetable(T, csv_file);

    fprintf('Generated %d samples from 3D Gaussian distribution\n', num_bubble);
    fprintf('Samples saved to: %s\n', csv_file);
    fprintf('\nSample Statistics:\n');
    fprintf('Mean X: %.4f, Y: %.4f, Z: %.4f, R: %.4f\n', mean(samples(1,:)), mean(samples(2,:)), mean(samples(3,:)), mean(samples(4,:)));
    fprintf('Std  X: %.4f, Y: %.4f, Z: %.4f, R: %.4f\n', std(samples(1,:)), std(samples(2,:)), std(samples(3,:)), std(samples(4,:)));
    
    min_liquid_thickness = config.simulation.min_dist_wall;
    min_liquid_thickness = min_liquid_thickness/inner_radius;
    slug_range = slug_range/inner_radius;
    %z_range = config.grid.Nz*config.grid.dz / inner_radius*1e3;
    z_range = 2;
    slug_pow_num = config.simulation.slug_pow_num;

    major_axis_length = slug_length/inner_radius;
    minor_axis_length = 1 - min_liquid_thickness;

    %ここで楕円の中心をどこに持ってくるか決める
    %上限は、円柱の下面と楕円の原点をそろえる場合か、slug_rangeの部分
    %下限は、円柱の中心面と楕円の頂点をそろえる場合
    %center_z = -rand*(major_axis_length - z_range/2) - z_range/2;
    center_z = 0;
    samples_slug = [center_z; major_axis_length; minor_axis_length; slug_pow_num];

    csv_file_slug = fullfile(save_path, 'sample_slug.csv');
    sample_table_slug = array2table(transpose(samples_slug), 'VariableNames', {'center z', 'major axis', 'minor axis', 'power number'});
    writetable(sample_table_slug, csv_file_slug);

    samples_tuple = {samples, samples_slug};
    csv_file_tuple = {csv_file, csv_file_slug};

    save_path = config.location_seedfiles_path;
    if ~exist(save_path, 'dir')
        mkdir(save_path);
    end
    plot_gen3d_gl(config_file, csv_file_tuple, '1', save_path);
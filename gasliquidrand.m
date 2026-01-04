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

% -------------------------------------------------------------------------
% 2) シミュレーション用グリッドの定義
% -------------------------------------------------------------------------
Nx = config.grid.Nx*2;
Ny = config.grid.Ny*2;
dx = config.grid.dx/2;
dy = config.grid.dy/2;
kgrid = kWaveGrid(Nx, dx, Ny, dy);
save_path = config.save_path;
t_end = config.simulation.t_end;
CFL = config.simulation.CFL;

% -------------------------------------------------------------------------
% 3) 媒質パラメータ
% -------------------------------------------------------------------------
% 水のパラメータ
medium.sound_speed = config.medium.water.sound_speed;
medium.density     = config.medium.water.density;
medium.alpha_coeff = config.medium.water.alpha_coeff;
medium.alpha_power = config.medium.water.alpha_power;
% medium.BonA = config.medium.water.BonA;

% ガラスのパラメータ
PMMA.sound_speed = config.medium.PMMA.sound_speed;
PMMA.density     = config.medium.PMMA.density;
PMMA.alpha_coeff = config.medium.PMMA.alpha_coeff;

% -------------------------------------------------------------------------
% 4) トランスデューサーの設定
% -------------------------------------------------------------------------
% トランスデューサーの位置と向きを設定
transducer_pos = config.transducer.position;
transducer_rot = config.transducer.rotation;
transducer_focus = config.transducer.focus;

% ガラス円環のマスクを作成
%glass_mask = zeros(Nx, Ny);
outer_radius = config.pipe.outer_radius;  % 外側の半径
inner_radius = config.pipe.inner_radius;   % 内側の半径
thickness = outer_radius - inner_radius;

% -------------------------------------------------------------------------
% 5) 塩ビパイプの設定
% -------------------------------------------------------------------------
% パイプの中心位置
center_x = 30 + round(config.pipe.center_x/dx);
center_y = Ny/2;

% 媒質パラメータのマスクを作成
medium.sound_speed = medium.sound_speed * ones(Nx, Ny);
medium.density = medium.density * ones(Nx, Ny);
medium.alpha_coeff = medium.alpha_coeff * ones(Nx, Ny);
medium.alpha_power = medium.alpha_power;
% medium.BonA = medium.BonA * ones(Nx, Ny);

% 塩ビパイプのパラメータを設定
pipe_mask=zeros(Nx,Ny);
for i = 1:Nx
    for j = 1:Ny
        if dx*sqrt((i-center_x)^2+(j-center_y)^2)<= outer_radius & dx*sqrt((i-center_x)^2+(j-center_y)^2)>= inner_radius
            medium.sound_speed(i, j) = PMMA.sound_speed;
            medium.density(i, j) = PMMA.density;
            medium.alpha_coeff(i, j) = PMMA.alpha_coeff;
            pipe_mask(i,j) = 1;
        end
    end
end

m = config.simulation.num_bubble;
bubble_r = 0.75e-3;
samples=zeros(2,m);
attempts=0;
max_attempts=1000;
count=0;
min_dist=config.simulation.distance_bubble+2*bubble_r;
min_dist_wall=config.simulation.min_dist_wall;

% % ーーーーーーーーーーーここから気泡生成ーーーーーーーーーーー
% while count<m & attempts<max_attempts
%     xy = mvnrnd([0, 0], eye(2), 1)';
%     candidate = xy*inner_radius;             % 3x1 vec
%     % Check if (X,Y) is inside unit circle
%     if (candidate(1))^2 + (candidate(2))^2 <= (inner_radius-min_dist-min_dist_wall)^2
%         % If this is the first sample, always accept
%         if count == 0
%             count = count + 1;
%             samples(:, count) = candidate;
%             disp("ok");
%         else
%             % Compute Euclidean distances to all previous samples
%             dists = sqrt(sum((samples(:,1:count) - candidate).^2, 1));
%             % Accept only if all distances are greater than or equal to min_dist
%             if all(dists >= min_dist)
%                 count = count + 1;
%                 samples(:, count) = candidate;
%                 attempts=0;
%                 disp("ok");
%             end
%         end
%     end
%     attempts = attempts + 1;
% end
% if count == m
%     % Calculate the Euclidean distance between all pairs of samples and display the minimum value
%     D = pdist(samples'); % pdist computes the pairwise distances between rows of the matrix
%     min_dist_val = min(D);
%     fprintf('Minimum Euclidean distance between samples: %.6f\n', min_dist_val);
%     disp(['Successfully generated ' num2str(m) ' samples! ']);
% end
% if count < m
%     error('Could not generate enough samples in %d attempts. Try reducing min_dist or increasing max_attempts.', max_attempts);
% end

% save_path="/mnt/sdb/yyamaguchi/kwavesource2d/bubblerand/loc_seed";
% csv_path=fullfile(save_path,'sample.csv');
% writematrix(samples',csv_path);

save_path="/mnt/sdb/yyamaguchi/kwavesource2d/bubblerand/loc_seed";
csv_path=fullfile(save_path,'sample.csv');
samples=readmatrix(csv_path);
samples=samples';

bubble_r = round(bubble_r/dx);
samples = round(samples./dx);
bubble_mask=zeros(Nx,Ny);

for k = 1:m
    for i = center_x-bubble_r+samples(1,k)-5:center_x+bubble_r+samples(1,k)+5
        for j = center_y-bubble_r+samples(2,k)-5:center_y+bubble_r+samples(2,k)+5
            if sqrt((i-center_x-samples(1,k))^2+(j-center_y-samples(2,k))^2) <= bubble_r
                medium.sound_speed(i, j)= config.medium.air.sound_speed;
                medium.density(i, j) = config.medium.air.density;
                medium.alpha_coeff(i, j) = config.medium.air.alpha_coeff;
                % medium.BonA(i,j) = config.medium.air.BonA;
                bubble_mask(i,j)=1;
            end
        end
    end
end

% -------------------------------------------------------------------------
% 6) シミュレーション時間配列の作成
% -------------------------------------------------------------------------
kgrid.makeTime(config.medium.water.sound_speed, CFL, t_end);

% -------------------------------------------------------------------------
% 7) ソース波形の設定
% -------------------------------------------------------------------------
distance_pipe_source = config.source.distance_pipe_source;
distance_pipe_source = distance_pipe_source/dy;
source_diameter = config.source.diameter;
source_diameter = source_diameter/dx;
source_mask=zeros(Nx,Ny);
source_mask(30,round(Ny/2-source_diameter/2):round(Ny/2+source_diameter/2))= 1;
source.u_mask = zeros(Nx, Ny);
source.u_mask(30,round(Ny/2-source_diameter/2):round(Ny/2+source_diameter/2))= 1;
source_signal = zeros(size(kgrid.t_array));
source_frequency = config.source.tone_burst_freq;
t_array = kgrid.t_array;

% Input signal properties
source_strength = config.source.source_strength;
tone_burst_freq = config.source.tone_burst_freq;
tone_burst_cycles = config.source.tone_burst_cycles;

% Generate input signal
% Read upsampled pulse from mat file
% original_pulse = load('/home/matsubara/Scripts/kwavesource/src/pulse_4000_upsampled.mat');
%input_signal = original_pulse.pulse_upsampled;

input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% Scale by impedance
input_signal = (source_strength / (config.medium.water.sound_speed * config.medium.water.density)) * input_signal(1:end-10);
source.ux = input_signal;

% -------------------------------------------------------------------------
% 8) センサーの設定
% -------------------------------------------------------------------------
sensor.mask = zeros(Nx, Ny);
sensor.mask(30,round(Ny/2-source_diameter/2):round(Ny/2+source_diameter/2))= 1;
sensor.record = {'p'};
save_data_path="/mnt/sdb/yyamaguchi/kwavesource2d/bubblerand/data";
save_pic_path="/mnt/sdb/yyamaguchi/kwavesource2d/bubblerand/pic";
% -------------------------------------------------------------------------
% 9) シミュレーションのオプション設定
% -------------------------------------------------------------------------
input_args = {
    'PMLInside', true, 'PlotPML', true, ...
    'PMLSize', config.simulation.pml_size, ...
    'RecordMovie', false, ...
    'PlotFreq', 50, ...
    'DataCast', DATA_CAST, ...
    };

% Convert to double for visualization
transmit_mask_double = double(source_mask);
pipe_mask_double     = double(pipe_mask);
bubble_mask_double    = double(bubble_mask);

composite_mask = transmit_mask_double * 1 + ...
                 pipe_mask_double * 2 + ...
                 bubble_mask_double * 3;

my_colormap = [
0.0 0.0 1.0;  % 0: 背景 (黒)
0.0 0.0 0.0;  % 1: transmit_mask (青)
0.0 1.0 0.0;  % 2: pipe_mask (緑)
1.0 1.0 1.0   % 3: bubble_mask (赤)
];

% プロット
figure;
imagesc(composite_mask); 
colormap(my_colormap); % 定義したカラーマップを適用
axis equal tight;
            
saveas(gcf, fullfile(save_pic_path, ['experimental_setup.png']));

% -------------------------------------------------------------------------
% 10) シミュレーション実行
% -------------------------------------------------------------------------
sensor_data = kspaceFirstOrder2DG(kgrid, medium, source, sensor, input_args{:});
save(fullfile(save_data_path, ['solid_liquid_reflector_linear' num2str(Nx) '.mat']), 'sensor_data', 'kgrid', '-v7.3');
% -------------------------------------------------------------------------
% 11) 結果の可視化
% -------------------------------------------------------------------------
p = gather(sensor_data.p);
sz = size(p);
plot_idx = round(sz(1)/2);

figure(1);
plot(kgrid.t_array * 1e6, p(plot_idx,:) * 1e-3, 'b-');
xlabel('Time [\mus]');
ylabel('Pressure [kPa]');
ylim([-400 400]);
title('Signal from Transducer transmit');
grid on;
saveas(gcf, fullfile(save_pic_path, ['signal_solid_liquid_reflector_linear' num2str(Nx) '.png']));
clearvars;
close all;
DATA_CAST='gpuArray-single';

config=jsondecode(fileread("/home/user01/Document/yyamaguchi/documents/simulationB4/sandbox/config.json"));

% Nx = config.grid.Nx;
% Ny = config.grid.Ny;
% dx = config.grid.dx;
% dy = config.grid.dy;
ppw=64;
Nx = 128*ppw;
Ny = 64*ppw;
dx = 1.5/4*1e-3/ppw;
dy = dx;
kgrid = kWaveGrid(Nx, dx, Ny, dy);
t_end = config.simulation.t_end;
cfl = config.simulation.CFL;

save_data_path = config.save_data_path;
save_pic_path = config.save_pic_path;

water.sound_speed = config.medium.water.sound_speed;
water.density     = config.medium.water.density;
water.alpha_coeff = config.medium.water.alpha_coeff;
water.alpha_power = config.medium.water.alpha_power;
water.BonA = config.medium.water.BonA;

medium.sound_speed = water.sound_speed * ones(Nx, Ny);
medium.density = water.density * ones(Nx, Ny);
medium.alpha_coeff = water.alpha_coeff * ones(Nx, Ny);
medium.alpha_power = water.alpha_power;
medium.BonA = water.BonA * ones(Nx, Ny);

air.sound_speed = 0;
air.density = 10;
air.alpha_coeff = config.medium.air.alpha_coeff;
air.BonA = config.medium.air.BonA;

glass.sound_speed = config.medium.glass.sound_speed;
glass.density = config.medium.glass.density;
glass.alpha_coeff = config.medium.glass.alpha_coeff;
glass.BonA = config.medium.glass.BonA;

medium.sound_speed_ref = water.sound_speed;

% config.json内のパイプの中心はtdx1からの相対座標！！
% tdx1の位置は(30, Ny/2)で固定！！
distance_tdx = config.pipe.distance_tdx;
inner_radius = config.pipe.inner_radius;
outer_radius = config.pipe.outer_radius;
pcx = 30 + round(distance_tdx/dx); % pipe center x in grid
% pcx=Nx-pcx;
pcy = Ny/2; % pipe center y in grid
outer_radius_sqrt2_grid = round(outer_radius/dx/sqrt(2));

pipe_mask = zeros(Nx, Ny);

bubble_diameter = config.bubble.diameter;
glass_diameter = config.glass.diameter;

% bubble_mask = makeDisc(Nx, Ny, pcx, pcy,...
%                     round(bubble_diameter/2/dx));
% glass_mask = zeros(Nx, Ny);
% glass_mask(pcx:end,:) = 1;

% medium.sound_speed(glass_mask==1) = glass.sound_speed;
% medium.density(glass_mask==1) = glass.density;
% medium.alpha_coeff(glass_mask==1) = glass.alpha_coeff;
% medium.BonA(glass_mask==1) = glass.BonA;

bubble_mask = zeros(Nx, Ny);
bubble_mask(pcx:end,:) = 1;
% bubble_mask(1:pcx,:)=1;

medium.sound_speed(bubble_mask==1) = air.sound_speed;
medium.density(bubble_mask==1) = air.density;
medium.alpha_coeff(bubble_mask==1) = air.alpha_coeff;
medium.BonA(bubble_mask==1) = air.BonA;

kgrid.makeTime(water.sound_speed, cfl, t_end);

% tdx1の位置は(30, Ny/2)で固定！！
% config.jsonのtdx2以上はパイプ中心からの相対座標！！

tone_burst_cycles = config.transducer1.tone_burst_cycles;
tone_burst_freq = config.transducer1.tone_burst_freq;

tdx1_pressure = config.transducer1.pressure;

tdx1_diameter = config.transducer1.diameter/2;
tdx1_radius_grid = round(tdx1_diameter/2/dx);
% tdx1_pos = [30, Ny/2];
% tdx1_focus_distance = config.transducer1.focus_distance;
% tdx1_steering_angle = config.transducer1.steering_angle;
% tdx1_mask = makeArc([Nx, Ny], tdx1_pos, ...
%                     round(2*tdx1_focus_distance/dx), ...
%                     round(tdx1_diameter/dx), ...
%                     [tdx1_pos(1)+round(tdx1_focus_distance/dx), tdx1_pos(2)]);
tdx1_mask = zeros(Nx, Ny);
tdx1_mask(30, 11:Ny-11) = 1;
% tdx1_mask(Nx-30, 11:Ny-11) = 1;
source.u_mask= tdx1_mask;

input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

input_signal = tdx1_pressure/water.sound_speed/water.density*input_signal;
source.ux = input_signal;
% source.ux = -input_signal;

sensor_mask = zeros(Nx, Ny);
sensor_mask(pcx-round(1.0e-3/dx), 11:Ny-11) = 1;
sensor_mask(pcx+round(0.3e-3/dx), 11:Ny-11) = 1;
% sensor_mask(pcx+round(1.0e-3/dx), 11:Ny-11) = 1;
% sensor_mask(pcx-round(0.3e-3/dx), 11:Ny-11) = 1;
% distance_tdx_sqrt2_grid = round(distance_tdx/sqrt(2)/dx);
% tdx2_mask = makeArc([Nx,Ny],...
%                     [pcx-distance_tdx_sqrt2_grid, pcy-distance_tdx_sqrt2_grid],...
%                     inf,2*tdx1_radius_grid+1,[pcx, pcy]);
% sensor_mask(tdx2_mask==1) = 1;
% tdx3_mask = makeArc([Nx,Ny],...
%                     [pcx+distance_tdx_sqrt2_grid, pcy-distance_tdx_sqrt2_grid],...
%                     inf,2*tdx1_radius_grid+1,[pcx, pcy]);
% sensor_mask(tdx3_mask==1) = 1;
% tdx4_mask = zeros(Nx, Ny);
% tdx4_mask(2*pcx-30, Ny/2-tdx1_radius_grid:Ny/2+tdx1_radius_grid) = 1;
sensor.mask = zeros(Nx, Ny);
sensor.mask(sensor_mask==1)=1;
sensor.record = {'p'};
num_sensor=sum(sensor_mask(:));
disp(num_sensor);

input_args = {
    'PMLInside', true, 'PlotPML', true, ...
    'PMLSize', config.simulation.pml_size, ...
    'RecordMovie', false, ...
    'PlotFreq', 50, ...
    'DataCast', DATA_CAST, ...
    };

sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
p=sensor_data.p;
tdx1_idx=1:2:num_sensor-1;
tdx2_idx=2:2:num_sensor;
p2=p(tdx1_idx,:);
p1=p(tdx2_idx,:);
p1=mean(p1);
p2=mean(p2);

save(fullfile(save_data_path, ['sandbox/air_lossy_vert_refined_void' num2str(ppw) '.mat']), 'p1', 'p2', 'kgrid', '-v7.3');

transmit_mask_double = double(tdx1_mask);
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
composite_mask = transpose(composite_mask);
imagesc(composite_mask); 
colormap(my_colormap); % 定義したカラーマップを適用
axis equal tight;
            
saveas(gcf, fullfile(save_pic_path, ['sandbox/experimental_setup_air_lossless_vert_void.png']));

figure(1);
plot(1e6*kgrid.t_array, p1);
hold on;
plot(1e6*kgrid.t_array, p2);
xlim([0,t_end*1e6]);
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossy_vert_refined_void' num2str(ppw) '.png']));
clearvars;
close all;
DATA_CAST='gpuArray-single';

config=jsondecode(fileread("/home/user01/Document/yyamaguchi/documents/simulationB4/sandbox/config.json"));

Nx = config.grid.Nx;
Ny = config.grid.Ny;
dx = config.grid.dx;
dy = config.grid.dy;
kgrid = kWaveGrid(Nx, dx, Ny, dy);
t_end = config.simulation.t_end;
cfl = config.simulation.CFL;

save_data_path = config.save_data_path;
save_pic_path = config.save_pic_path;

water.sound_speed = config.medium.water.sound_speed;
water.sound_speed_shear = config.medium.water.sound_speed_shear;
water.density     = config.medium.water.density;

medium.sound_speed_compression = water.sound_speed * ones(Nx, Ny);
medium.sound_speed_shear = water.sound_speed_shear * ones(Nx, Ny);
medium.density = water.density * ones(Nx, Ny);

air.sound_speed = config.medium.air.sound_speed;
air.density = config.medium.air.density;

glass.sound_speed = config.medium.glass.sound_speed;
glass.sound_speed_shear = config.medium.glass.sound_speed_shear;
glass.density = config.medium.glass.density;

% config.json内のパイプの中心はtdx1からの相対座標！！
% tdx1の位置は(30, Ny/2)で固定！！
distance_tdx = config.pipe.distance_tdx;
pcx = 30 + round(distance_tdx/dx); % pipe center x in grid
pcy = Ny/2; % pipe center y in grid

bubble_diameter = config.bubble.diameter;
glass_diameter = config.glass.diameter;

glass_mask = makeDisc(Nx, Ny, pcx, pcy,...
                    round(glass_diameter/2/dx));

medium.sound_speed_compression(glass_mask==1) = glass.sound_speed;
medium.sound_speed_shear(glass_mask==1) = glass.sound_speed_shear;
medium.density(glass_mask==1) = glass.density;

kgrid.makeTime(medium.sound_speed_compression, cfl, t_end);

% tdx1の位置は(30, Ny/2)で固定！！
% config.jsonのtdx2以上はパイプ中心からの相対座標！！

tone_burst_cycles = config.transducer1.tone_burst_cycles;
tone_burst_freq = config.transducer1.tone_burst_freq;

tdx1_pressure = config.transducer1.pressure;

tdx1_diameter = config.transducer1.diameter;
tdx1_radius_grid = round(tdx1_diameter/2/dx);
% tdx1_pos = [30, Ny/2];
% tdx1_focus_distance = config.transducer1.focus_distance;
% tdx1_steering_angle = config.transducer1.steering_angle;
% tdx1_mask = makeArc([Nx, Ny], tdx1_pos, ...
%                     round(2*tdx1_focus_distance/dx), ...
%                     round(tdx1_diameter/dx), ...
%                     [tdx1_pos(1)+round(tdx1_focus_distance/dx), tdx1_pos(2)]);
tdx1_mask = zeros(Nx, Ny);
tdx1_mask(30, Ny/2-tdx1_radius_grid:Ny/2+tdx1_radius_grid) = 1;
source.u_mask= tdx1_mask;

input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

input_signal = tdx1_pressure/water.sound_speed/water.density*input_signal;
source.ux = input_signal;

sensor_mask = zeros(Nx, Ny);
sensor_mask(tdx1_mask==1) = 1;
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

input_args = {
    'PMLInside', true, 'PlotPML', true, ...
    'PMLSize', config.simulation.pml_size, ...
    'RecordMovie', false, ...
    'PlotFreq', 50, ...
    'DataCast', DATA_CAST, ...
    };

sensor_data = pstdElastic2D(kgrid, medium, source, sensor, input_args{:});
save(fullfile(save_data_path, ['sandbox/glass_elastic' num2str(Nx) '.mat']), 'sensor_data', 'kgrid', '-v7.3');

p = sensor_data.p;
p = mean(p);

figure(1);
plot(1e6*kgrid.t_array, p);
xlim([0,t_end*1e6]);
saveas(gcf, fullfile(save_pic_path, ['sandbox/glass_elastic' num2str(Nx) '.png']));
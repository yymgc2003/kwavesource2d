clearvars;
close all;
DATA_CAST='gpuArray-single';

config=jsondecode(fileread("/home/user01/Document/yyamaguchi/documents/simulationB4/sandbox/config.json"));

% Nx = config.grid.Nx;
% Ny = config.grid.Ny;
% dx = config.grid.dx;
% dy = config.grid.dy;
% ppw=2;
for ppw=2:1:15
    Nx = 64*ppw;
    dx = 1.5/4*1e-3/ppw;
    kgrid = kWaveGrid(Nx, dx);
    t_end = 4e-6;
    cfl = config.simulation.CFL;
    surf = 2e-3;

    save_data_path = config.save_data_path;
    save_pic_path = config.save_pic_path;

    water.sound_speed = config.medium.water.sound_speed;
    water.density     = config.medium.water.density;
    water.alpha_coeff = config.medium.water.alpha_coeff;
    water.alpha_power = config.medium.water.alpha_power;
    water.BonA = config.medium.water.BonA;

    medium.sound_speed = water.sound_speed * ones(Nx,1);
    medium.density = water.density * ones(Nx,1);
    % medium.alpha_coeff = water.alpha_coeff * ones(Nx, Ny);
    % medium.alpha_power = water.alpha_power;
    % medium.BonA = water.BonA * ones(Nx, Ny);

    air.sound_speed = config.medium.air.sound_speed;
    air.density = config.medium.air.density;
    air.alpha_coeff = config.medium.air.alpha_coeff;
    air.BonA = config.medium.air.BonA;

    glass.sound_speed = config.medium.glass.sound_speed;
    glass.density = config.medium.glass.density;
    glass.alpha_coeff = config.medium.glass.alpha_coeff;
    glass.BonA = config.medium.glass.BonA;

    medium.sound_speed_ref = water.sound_speed;

    pcx = 30 + round(surf/dx); % pipe center x in grid
    % glass_mask = zeros(Nx, Ny);
    % glass_mask(pcx:end,:) = 1;

    % medium.sound_speed(glass_mask==1) = glass.sound_speed;
    % medium.density(glass_mask==1) = glass.density;
    % medium.alpha_coeff(glass_mask==1) = glass.alpha_coeff;
    % medium.BonA(glass_mask==1) = glass.BonA;

    bubble_mask = zeros(Nx,1);
    bubble_mask(pcx:end,:) = 1;

    medium.sound_speed(bubble_mask==1) = glass.sound_speed;
    medium.density(bubble_mask==1) = glass.density;

    kgrid.makeTime(water.sound_speed, cfl, t_end);

    tdx1_mask = zeros(Nx,1);
    tdx1_mask(30,:) = 1;
    source.u_mask= tdx1_mask;

    tone_burst_cycles = config.transducer1.tone_burst_cycles;
    tone_burst_freq = config.transducer1.tone_burst_freq;
    tdx1_pressure = config.transducer1.pressure;
    input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

    input_signal = tdx1_pressure/water.sound_speed/water.density*input_signal;
    source.ux = input_signal;

    sensor_mask = zeros(Nx,1);
    sensor_mask(pcx-round(1.0e-3/dx),:) = 1;
    sensor_mask(pcx+round(0.5e-3/dx),:) = 1;
    sensor.mask = zeros(Nx,1);
    sensor.mask(sensor_mask==1)=1;
    sensor.record = {'p'};

    input_args = {
        'PMLInside', true, 'PlotPML', true, ...
        'PMLSize', 20, ...
        'RecordMovie', false, ...
        'PlotFreq', 50, ...
        'DataCast', DATA_CAST, ...
        };

    sensor_data = kspaceFirstOrder1D(kgrid, medium, source, sensor, input_args{:});
    p=sensor_data.p;
    p1=p(1,:);
    p2=p(2,:);

    save(fullfile(save_data_path, ['sandbox/glass_lossless1D_ppw' num2str(ppw) '.mat']), 'p1', 'p2', 'kgrid', '-v7.3');

    figure(1);
    plot(1e6*kgrid.t_array, p1);
    hold on;
    plot(1e6*kgrid.t_array, p2);
    xlim([0,t_end*1e6]);
    saveas(gcf, fullfile(save_pic_path, ['sandbox/glass_lossless1D_ppw' num2str(ppw) '.png']));
    figure(close);
end
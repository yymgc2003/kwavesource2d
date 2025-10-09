function plot_gen(config_file, location_csv, locnum_str, ... 
    device_num, save_full_path)

    config = jsondecode(fileread(config_file));
    save_logs_path = fullfile(save_full_path, 'logs');
    save_data_path = fullfile(save_full_path, 'data');

    location = csvread(location_csv);
    glass_radius = config.simulation.glass_radius;

    figure;
    hold on;
    for i = 1:size(location,1)
        loc_seed = location(i,:);
        loc_seed(1) = config.pipe.inner_radius * loc_seed(1);
        loc_seed(2) = config.pipe.inner_radius * loc_seed(2);
        glass_radius_2d = glass_radius^2 - (Nz*dz*(0.5-loc_seed(3)))^2;
        if glass_radius_2d>0
            glass_radius_2d = sqrt(glass_radius_2d);
            rectangle('Position',[loc_seed(1)-glass_radius_2d,...
                loc_seed(2)-glass_radius_2d,...
                2*glass_radius_2d,2*glass_radius_2d], ...
                'Curvature', [1 1], ...
                'FaceColor', 'b', ...
                'EdgeColor', 'b');
            hold on;
        end
    end
    rectangle('Position',[-16e-3,-16e-3,32e-3,32e-3], ...
        'Curvature', [1 1], ...
        'EdgeColor', 'black');
    hold on;
    rectangle('Position',[-13e-3,-13e-3,26e-3,26e-3], ...
        'Curvature', [1 1], ...
        'EdgeColor', 'black');
    hold on;
    axis equal;
    xlim([-19e-3 19e-3]);
    ylim([-19e-3 19e-3]);
    axis off;
    saveas(gcf, fullfile(save_logs_path, ['experimental_setup' locnum_str '.png']));
end
function plot_gen(config_file, location_csv, locnum_str, ... 
    save_full_path)

    config = jsondecode(fileread(config_file));
    save_logs_path = fullfile(save_full_path, 'logs');
    save_data_path = fullfile(save_full_path, 'data');

    location = csvread(location_csv);
    glass_radius = config.simulation.glass_radius;
    Nz = config.grid.Nz;
    dz = config.grid.dz;

    figure;
    hold on;
    for i = 1:size(location,1)
        loc_seed = location(i,:);
        loc_seed(1) = config.pipe.inner_radius * loc_seed(1);
        loc_seed(2) = config.pipe.inner_radius * loc_seed(2);
        glass_radius = loc_seed(4)*config.pipe.inner_radius/2;
        glass_radius_2d = (glass_radius)^2 - (config.pipe.inner_radius*(loc_seed(3)))^2;
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
    rectangle('Position',[-16,-16,32,32], ...
        'Curvature', [1 1], ...
        'EdgeColor', 'black');
    hold on;
    rectangle('Position',[-13,-13,26,26], ...
        'Curvature', [1 1], ...
        'EdgeColor', 'black');
    hold on;
    axis equal;
    xlim([-19 19]);
    ylim([-19 19]);
    axis off;
    saveas(gcf, fullfile(save_logs_path, ['experimental_setup' locnum_str '.png']));
end
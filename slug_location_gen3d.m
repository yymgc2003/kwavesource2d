function samples = slug_location_gen3d(slug_length, slug_range)
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

    cur_gas_fraction = 0;
    inner_radius = config.pipe.inner_radius;
    min_liquid_thickness = config.simulation.min_dist_wall;
    min_liquid_thickness = min_liquid_thickness/inner_radius;
    slug_range = slug_range/inner_radius;
    z_range = config.grid.Nz*config.grid.dz / inner_radius*1e3;
    slug_pow_num = config.simulation.slug_pow_num;

    major_axis_length = slug_length/inner_radius;
    minor_axis_length = 1 - min_liquid_thickness;

    %ここで楕円の中心をどこに持ってくるか決める
    %上限は、円柱の下面と楕円の原点をそろえる場合か、slug_rangeの部分
    %下限は、円柱の中心面と楕円の頂点をそろえる場合
    center_z = rand*( -slug_range + z_range/inner_radius/2);
    samples = [center_z; major_axis_length; minor_axis_length; slug_pow_num];

    csv_file = fullfile(save_path, 'sample.csv');
    sample_table = array2table(transpose(samples), 'VariableNames', {'CZ', 'MAJ_AXIS', 'MIN_AXIS', 'POW'});
    writetable(sample_table, csv_file);

    fprintf('Samples saved to: %s\n', csv_file);
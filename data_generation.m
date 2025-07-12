function signalgen_module_all()
% This function reads config.json and sequentially runs signalgen_module for all location1.csv, location2.csv, ... in the specified directory.

    % Path to config file
    config_file = 'config.json';

    % Get location_seedfiles_path from config.json
    config = jsondecode(fileread(config_file));
    location_dir = config.location_seedfiles_path;

    % Get all location*.csv files
    files = dir(fullfile(location_dir, 'location*.csv'));
    fprintf('--- %d location*.csv files found. ---\n', length(files));
    if isempty(files)
        error('No location*.csv files found in %s.', location_dir);
    end

    % Run signalgen_module for each CSV file
    for i = 1:length(files)
        location_csv = fullfile(location_dir, files(i).name);
        fprintf('--- Processing %s ---\n', location_csv);
        simulation_execution(config_file, location_csv);
    end
end
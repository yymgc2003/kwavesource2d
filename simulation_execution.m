function simulation_execution(config_file, location_csv, ...
    device_num, save_full_path)
% Wrapper to extract location number and call kwavesim for simulation.

    % Extract location number from file name
    [~, locname, ~] = fileparts(location_csv);
    locnum = regexp(locname, '\d+', 'match');
    if isempty(locnum)
        locnum_str = '';
    else
        locnum_str = locnum{1};
    end
    %save_full_path = config.save_full_path;

    % Call the main simulation function
    kwavesim(config_file, location_csv, locnum_str, ...
    device_num, save_full_path);
end
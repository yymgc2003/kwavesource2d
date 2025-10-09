function simulation_execution3d_gl(config_file, location_csv)
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
        kwavesim3d_gl(config_file, location_csv, locnum_str);
    end
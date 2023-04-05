function [] = launch_simulation(cfg)
    if (ispc)
        bin_name = 'mcxyzn.exe';
    elseif (ismac)
        bin_name = 'mcxyzn.mac';
    elseif (isunix)
        bin_name = 'mcxyzn.linux';
    else
        fprintf('Could not find the appropriate binary \n');
    end
    
    filepath = split(cfg.name,'/');
    name = filepath(end);
    name = split(name,'\');
    name = name(end);
    
    old_folder = cd('./data_files/outputs/mcxyzn/');
    system_command_string = string(bin_name)+" "+name;
    
    [status] = system(system_command_string);
    cd(old_folder);

end

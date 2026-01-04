mat_path = "/mnt/yamaguchi/kwavesource_gl/csv";

mat_files = dir(fullfile(mat_path, '*.mat'));

for t = 1:length(mat_files)
    mat_file = mat_files(t);
    mat_file_path = fullfile(mat_path, mat_file.name)
    data = load(mat_file_path);
    var_name = fieldnames(data);
    disp(var_name);

    p = data.p_cpu;

    Nx = size(p, 1);
    Ny = size(p, 2);
    Nz = size(p, 3);

    disp(Nx);
    disp(Ny);
    disp(Nz);

    dx = 0.042;
    dy = 0.042;
    dz = 0.042;

    x = dx/2:dx:(Nx-1/2)*dx;
    y = dy/2:dy:(Ny-1/2)*dy;
    z = dz/2:dz:(Nz-1/2)*dz;

    save_path = "/mnt/yamaguchi/kwavesource_gl/csv";
    file_id = fopen(fullfile(save_path, sprintf('pressure_xy%d.xy', t)), 'w');
    for i = 1:Nx
        for j = 1:Ny
            fprintf(file_id, '%.3f %.3f %.2f\n', x(i), y(j), p(i, j, Nz/2));
        end
        fprintf(file_id, '\n');
    end

    fclose(file_id);

    file_id = fopen(fullfile(save_path, sprintf('pressure_xz%d.xy', t)), 'w');
    for i = 1:Nx
        for j = 1:Nz
            fprintf(file_id, '%.3f %.3f %.2f\n', x(i), z(j), p(i, Ny/2, j));
        end
        fprintf(file_id, '\n');
    end

    fclose(file_id);
end


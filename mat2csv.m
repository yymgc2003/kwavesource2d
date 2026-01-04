clearvars;
close all;

mat_path = '/mnt/sdb/yyamaguchi/psdata2matlab/simulation/rawsignal/solidliquid/data/solid_liquid_reflector1.mat';

mat_load = load(mat_path);
disp(fieldnames(mat_load.sensor_data));
p = mat_load.sensor_data;
disp(size(mean(p)));
p = mean(p);

save_path = '/mnt/sdb/yyamaguchi/psdata2matlab/simulation/rawsignal/solidliquid/data/solid_liquid_reflector1.txt';
file_id = fopen(save_path, 'w');
for i=1:5000
    fprintf(file_id, '%.2f %.5e\n', i/50, p(1,i*20));
end
fclose(file_id);
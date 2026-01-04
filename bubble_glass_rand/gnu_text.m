config_file = 'config.json';

config = jsondecode(fileread(config_file));

save_dir = config.save_path;
[filepath, filename, ext] = fileparts(mfilename('fullpath'));

dataset_name = strsplit(filepath, '/');
dataset_name = string(dataset_name(end));
% dataset_name = glass_rand
save_dir = fullfile(save_dir, dataset_name);
disp(save_dir);

num_particles = 55;

save_path = sprintf("%s/case_test%d",save_dir,num_particles);
loc_dir =fullfile(save_path,"location_seed");
gnu_dir =fullfile(save_path,"gnu_csv");
locnum =1;

loctable = readtable(sprintf("%s/location%d.csv", loc_dir, locnum));
loc = table2array(loctable);
disp(length(loc));
loc(:,1:4)=loc(:,1:4)*1e3;
loc(:,3:4)=loc(:,3:4)/2;

fid = fopen(sprintf("%s/gnu_text%d.txt", loc_dir, locnum), "w");
fprintf(fid,'plot ');
for ii=1:length(loc)
    x0=loc(ii,1);y0=loc(ii,2);aa=loc(ii,3);bb=loc(ii,4);th=loc(ii,5);
    fprintf(fid,'(%.5g)*((%.5g)*cos(t))-(%.5g)*((%.5g)*sin(t))+(%.5g),\\\n',...
        cos(th),aa,sin(th),bb,x0);
    fprintf(fid,'(%.5g)*((%.5g)*cos(t))+(%.5g)*((%.5g)*sin(t))+(%.5g)\\\n',...
        sin(th),aa,cos(th),bb,y0);
    fprintf(fid,'lw 3 lc "white",\\\n');
end

fclose(fid);

clearvars;
close all;

config=jsondecode(fileread("/home/user01/Document/yyamaguchi/documents/simulationB4/sandbox/config.json"));
save_data_path = config.save_data_path;
save_pic_path = config.save_pic_path;

Z2=343*1.2;
Z1=1500*1000;
rrt=(Z2-Z1)/(Z2+Z1);
trt=2*Z2/(Z2+Z1);

t_array=[];
rratio=[];
tratio=[];
p=[];
t=[];
ppwlist=[];
count=0;
max_len=0;
simtype='c343rho1';

for i=20:-0.25:2
    if exist(fullfile(save_data_path, ['sandbox/air_lossless1D' simtype '_ppw' num2str(i) '.mat']))
        mat_struct = load(fullfile(save_data_path, ['sandbox/air_lossless1D' simtype '_ppw' num2str(i) '.mat']));
        cur_t_array=mat_struct.kgrid.t_array;
        cur_t_array=transpose(cur_t_array);
        cur_len = length(cur_t_array);
        ppwlist = [ppwlist, i];
        p1=mat_struct.p1;
        p1=transpose(p1);
        t1=mat_struct.p2;
        t1=transpose(t1);
        if count~=0
            t_array=[t_array,[cur_t_array;cur_t_array(end)*ones(max_len-cur_len,1)]];
            p=[p,[p1;p1(end)*ones(max_len-cur_len,1)]];
            t=[t,[t1;t1(end)*ones(max_len-cur_len,1)]];
        else
            t_array=[t_array,cur_t_array];
            p=[p,p1];
            t=[t,t1];
            max_len=cur_len;
        end
        A1 = max(p1(1:round(end*0.5)))-min(p1(1:round(end*0.5)));
        A2 = max(p1(round(end*0.5):end))-min(p1(round(end*0.5):end));
        A3 = max(t1)-min(t1);
        rratio = [rratio; -A2/A1];
        tratio = [tratio; A3/A1];
        count = count+1;
    end
end

writematrix([rratio,tratio], fullfile(save_data_path, ['sandbox/air_lossless1D' simtype '_ratio.txt']));

rratio = abs(rratio-rrt);
tratio = abs(tratio-trt);

figure('Position', [0,0,600,300]);
plot(ppwlist, rratio(:,1),'-o','LineWidth',2);
hold on;
plot(ppwlist, tratio(:,1),'-o','LineWidth',2);
fontsize(20,"points");
set(gca, "LineWidth",2);
set(gca, "YScale", "log");
xlim([2,20]);
xlabel("Points per Wavelength");
ylabel("Error");
legend({"Reflection","Transmission"},'Location','southeast');
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossless1D_conv_' simtype '.png']));

mat_struct=load("/mnt/sdb/yyamaguchi/simulationB4/data/sandbox/glass_lossless1D_ppw10.mat");
p1=mat_struct.p1;
p2=mat_struct.p2;
kgrid=mat_struct.kgrid;
figure('Position', [0,0,600,300]);
plot(1e6*kgrid.t_array, p1*1e-3,'LineWidth',2);
hold on;
plot(1e6*kgrid.t_array, p2*1e-3,'LineWidth',1.5);
xlim([0,4]);
ylabel("Pressure (kPa)");
xlabel("Time (µs)");
legend({"Reflection","Transmission"},'Location','southeast');
fontsize(16,"points");
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossless1Dppw10.png']));

Z2=5789*2240;
Z1=1500*1000;
rrt=(Z2-Z1)/(Z2+Z1);
trt=2*Z2/(Z2+Z1);

t_array=[];
rratio=[];
tratio=[];
p=[];
t=[];
ppwlist=[];
count=0;
max_len=0;
simtype='c343rho1';

for i=20:-0.25:2
    if exist(fullfile(save_data_path, ['sandbox/glass_lossless1D_ppw' num2str(i) '.mat']))
        mat_struct = load(fullfile(save_data_path, ['sandbox/glass_lossless1D_ppw' num2str(i) '.mat']));
        cur_t_array=mat_struct.kgrid.t_array;
        cur_t_array=transpose(cur_t_array);
        cur_len = length(cur_t_array);
        ppwlist = [ppwlist, i];
        p1=mat_struct.p1;
        p1=transpose(p1);
        t1=mat_struct.p2;
        t1=transpose(t1);
        if count~=0
            t_array=[t_array,[cur_t_array;cur_t_array(end)*ones(max_len-cur_len,1)]];
            p=[p,[p1;p1(end)*ones(max_len-cur_len,1)]];
            t=[t,[t1;t1(end)*ones(max_len-cur_len,1)]];
        else
            t_array=[t_array,cur_t_array];
            p=[p,p1];
            t=[t,t1];
            max_len=cur_len;
        end
        A1 = max(p1(1:round(end*0.5)))-min(p1(1:round(end*0.5)));
        A2 = max(p1(round(end*0.5):end))-min(p1(round(end*0.5):end));
        A3 = max(t1)-min(t1);
        rratio = [rratio; A2/A1];
        tratio = [tratio; A3/A1];
        count = count+1;
    end
end

writematrix([rratio,tratio], fullfile(save_data_path, ['sandbox/glass_lossless1D_ratio.txt']));

rratio = abs(rratio-rrt);
tratio = abs(tratio-trt);

figure('Position', [0,0,600,300]);
plot(ppwlist, rratio(:,1),'-o','LineWidth',2);
hold on;
plot(ppwlist, tratio(:,1),'-o','LineWidth',2);
fontsize(20,"points");
set(gca, "LineWidth",2);
set(gca, "YScale", "log");
xlim([2,15]);
xlabel("Points per Wavelength");
ylabel("Error");
legend({"Reflection","Transmission"},'Location','northeast');
saveas(gcf, fullfile(save_pic_path, ['sandbox/glass_lossless1D_conv.png']));

Z2=343*1.2;
Z1=5789*2240;
rrt=(Z2-Z1)/(Z2+Z1);
trt=2*Z2/(Z2+Z1);

t_array=[];
rratio=[];
tratio=[];
p=[];
t=[];
ppwlist=[];
count=0;
max_len=0;
simtype='c343rho1';

for i=15:-0.25:2
    if exist(fullfile(save_data_path, ['sandbox/airglass_lossless1D' simtype '_ppw' num2str(i) '.mat']))
        mat_struct = load(fullfile(save_data_path, ['sandbox/airglass_lossless1D' simtype '_ppw' num2str(i) '.mat']));
        cur_t_array=mat_struct.kgrid.t_array;
        cur_t_array=transpose(cur_t_array);
        cur_len = length(cur_t_array);
        ppwlist = [ppwlist, i];
        p1=mat_struct.p1;
        p1=transpose(p1);
        t1=mat_struct.p2;
        t1=transpose(t1);
        if count~=0
            t_array=[t_array,[cur_t_array;cur_t_array(end)*ones(max_len-cur_len,1)]];
            p=[p,[p1;p1(end)*ones(max_len-cur_len,1)]];
            t=[t,[t1;t1(end)*ones(max_len-cur_len,1)]];
        else
            t_array=[t_array,cur_t_array];
            p=[p,p1];
            t=[t,t1];
            max_len=cur_len;
        end
        A1 = max(p1(1:round(end*0.5)))-min(p1(1:round(end*0.5)));
        A2 = max(p1(round(end*0.5):round(end*0.9)))-min(p1(round(end*0.5):round(end*0.9)));
        A3 = max(t1)-min(t1);
        rratio = [rratio; -A2/A1];
        tratio = [tratio; A3/A1];
        count = count+1;
    end
end

writematrix([rratio,tratio], fullfile(save_data_path, ['sandbox/airglass_lossless1D' simtype '_ratio.txt']));

rratio = abs(rratio-rrt);
tratio = abs(tratio-trt);

figure('Position', [0,0,600,300]);
plot(ppwlist, rratio(:,1),'-o','LineWidth',2);
hold on;
plot(ppwlist, tratio(:,1),'-o','LineWidth',2);
fontsize(20,"points");
set(gca, "LineWidth",2);
set(gca, "YScale", "log");
xlim([2,15]);
xlabel("Points per Wavelength");
ylabel("Error");
legend({"Reflection","Transmission"},'Location','east');
saveas(gcf, fullfile(save_pic_path, ['sandbox/airglass_lossless1D_conv_' simtype '.png']));

Z2=343*10;
Z1=1500*1000;
rrt=(Z2-Z1)/(Z2+Z1);
trt=2*Z2/(Z2+Z1);

t_array=[];
rratio=[];
tratio=[];
p=[];
t=[];
ppwlist=[];
count=0;
max_len=0;
simtype='c343rho10';

for i=10:-0.25:2
    if exist(fullfile(save_data_path, ['sandbox/air_lossless1D' simtype '_ppw' num2str(i) '.mat']))
        mat_struct = load(fullfile(save_data_path, ['sandbox/air_lossless1D' simtype '_ppw' num2str(i) '.mat']));
        cur_t_array=mat_struct.kgrid.t_array;
        cur_t_array=transpose(cur_t_array);
        cur_len = length(cur_t_array);
        ppwlist = [ppwlist, i];
        p1=mat_struct.p1;
        p1=transpose(p1);
        t1=mat_struct.p2;
        t1=transpose(t1);
        if count~=0
            t_array=[t_array,[cur_t_array;cur_t_array(end)*ones(max_len-cur_len,1)]];
            p=[p,[p1;p1(end)*ones(max_len-cur_len,1)]];
            t=[t,[t1;t1(end)*ones(max_len-cur_len,1)]];
        else
            t_array=[t_array,cur_t_array];
            p=[p,p1];
            t=[t,t1];
            max_len=cur_len;
        end
        A1 = max(p1(1:round(end*0.5)))-min(p1(1:round(end*0.5)));
        A2 = max(p1(round(end*0.5):end))-min(p1(round(end*0.5):end));
        A3 = max(t1)-min(t1);
        rratio = [rratio; -A2/A1];
        tratio = [tratio; A3/A1];
        count = count+1;
    end
end

writematrix([rratio,tratio], fullfile(save_data_path, ['sandbox/air_lossless1D' simtype '_ratio.txt']));

rratio = abs(rratio-rrt);
tratio = abs(tratio-trt);

figure('Position', [0,0,600,300]);
plot(ppwlist, rratio(:),'-o','LineWidth',2);
hold on;
plot(ppwlist, tratio(:),'-o','LineWidth',2);
fontsize(20,"points");
set(gca, "LineWidth",2);
set(gca, "YScale", "log");
xlim([2,10]);
xlabel("Points per Wavelength");
ylabel("Error");
legend({"Reflection","Transmission"},'Location','southeast');
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossless1D_conv_' simtype '.png']));

Z2=343*100;
Z1=1500*1000;
rrt=(Z2-Z1)/(Z2+Z1);
trt=2*Z2/(Z2+Z1);

t_array=[];
rratio=[];
tratio=[];
p=[];
t=[];
ppwlist=[];
count=0;
max_len=0;
simtype='c343rho100';

for i=10:-0.25:2
    if exist(fullfile(save_data_path, ['sandbox/air_lossless1D' simtype '_ppw' num2str(i) '.mat']))
        mat_struct = load(fullfile(save_data_path, ['sandbox/air_lossless1D' simtype '_ppw' num2str(i) '.mat']));
        cur_t_array=mat_struct.kgrid.t_array;
        cur_t_array=transpose(cur_t_array);
        cur_len = length(cur_t_array);
        ppwlist = [ppwlist, i];
        p1=mat_struct.p1;
        p1=transpose(p1);
        t1=mat_struct.p2;
        t1=transpose(t1);
        if count~=0
            t_array=[t_array,[cur_t_array;cur_t_array(end)*ones(max_len-cur_len,1)]];
            p=[p,[p1;p1(end)*ones(max_len-cur_len,1)]];
            t=[t,[t1;t1(end)*ones(max_len-cur_len,1)]];
        else
            t_array=[t_array,cur_t_array];
            p=[p,p1];
            t=[t,t1];
            max_len=cur_len;
        end
        A1 = max(p1(1:round(end*0.5)))-min(p1(1:round(end*0.5)));
        A2 = max(p1(round(end*0.5):end))-min(p1(round(end*0.5):end));
        A3 = max(t1)-min(t1);
        rratio = [rratio; -A2/A1];
        tratio = [tratio; A3/A1];
        count = count+1;
    end
end

writematrix([rratio,tratio], fullfile(save_data_path, ['sandbox/air_lossless1D' simtype '_ratio.txt']));

rratio = abs(rratio-rrt);
tratio = abs(tratio-trt);

figure('Position', [0,0,600,300]);
plot(ppwlist, rratio(:,1)/10,'-o','LineWidth',2);
hold on;
plot(ppwlist, tratio(:,1),'-o','LineWidth',2);
fontsize(20,"points");
set(gca, "LineWidth",2);
set(gca, "YScale", "log");
xlim([2,10]);
xlabel("Points per Wavelength");
ylabel("Error");
legend({"Reflection","Transmission"},'Location','northeast');
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossless1D_conv_' simtype '.png']));

Z2=1200*10;
Z1=1500*1000;
rrt=(Z2-Z1)/(Z2+Z1);
trt=2*Z2/(Z2+Z1);

t_array=[];
rratio=[];
tratio=[];
p=[];
t=[];
ppwlist=[];
count=0;
max_len=0;
simtype='c1200rho10';

for i=10:-0.25:2
    if exist(fullfile(save_data_path, ['sandbox/air_lossless1D' simtype '_ppw' num2str(i) '.mat']))
        mat_struct = load(fullfile(save_data_path, ['sandbox/air_lossless1D' simtype '_ppw' num2str(i) '.mat']));
        cur_t_array=mat_struct.kgrid.t_array;
        cur_t_array=transpose(cur_t_array);
        cur_len = length(cur_t_array);
        ppwlist = [ppwlist, i];
        p1=mat_struct.p1;
        p1=transpose(p1);
        t1=mat_struct.p2;
        t1=transpose(t1);
        if count~=0
            t_array=[t_array,[cur_t_array;cur_t_array(end)*ones(max_len-cur_len,1)]];
            p=[p,[p1;p1(end)*ones(max_len-cur_len,1)]];
            t=[t,[t1;t1(end)*ones(max_len-cur_len,1)]];
        else
            t_array=[t_array,cur_t_array];
            p=[p,p1];
            t=[t,t1];
            max_len=cur_len;
        end
        A1 = max(p1(1:round(end*0.5)))-min(p1(1:round(end*0.5)));
        A2 = max(p1(round(end*0.5):end))-min(p1(round(end*0.5):end));
        A3 = max(t1)-min(t1);
        rratio = [rratio; -A2/A1];
        tratio = [tratio; A3/A1];
        count = count+1;
    end
end

writematrix([rratio,tratio], fullfile(save_data_path, ['sandbox/air_lossless1D' simtype '_ratio.txt']));

rratio = abs(rratio-rrt);
tratio = abs(tratio-trt);

figure('Position', [0,0,600,300]);
plot(ppwlist, rratio(:,1)/10,'-o','LineWidth',2);
hold on;
plot(ppwlist, tratio(:,1),'-o','LineWidth',2);
fontsize(20,"points");
set(gca, "LineWidth",2);
set(gca, "YScale", "log");
xlim([2,10]);
xlabel("Points per Wavelength");
ylabel("Error");
legend({"Reflection","Transmission"},'Location','northeast');
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossless1D_conv_' simtype '.png']));

Z2=600*10;
Z1=1500*1000;
rrt=(Z2-Z1)/(Z2+Z1);
trt=2*Z2/(Z2+Z1);

t_array=[];
rratio=[];
tratio=[];
p=[];
t=[];
ppwlist=[];
count=0;
max_len=0;
simtype='c600rho10';

for i=10:-0.25:2
    if exist(fullfile(save_data_path, ['sandbox/air_lossless1D' simtype '_ppw' num2str(i) '.mat']))
        mat_struct = load(fullfile(save_data_path, ['sandbox/air_lossless1D' simtype '_ppw' num2str(i) '.mat']));
        cur_t_array=mat_struct.kgrid.t_array;
        cur_t_array=transpose(cur_t_array);
        cur_len = length(cur_t_array);
        ppwlist = [ppwlist, i];
        p1=mat_struct.p1;
        p1=transpose(p1);
        t1=mat_struct.p2;
        t1=transpose(t1);
        if count~=0
            t_array=[t_array,[cur_t_array;cur_t_array(end)*ones(max_len-cur_len,1)]];
            p=[p,[p1;p1(end)*ones(max_len-cur_len,1)]];
            t=[t,[t1;t1(end)*ones(max_len-cur_len,1)]];
        else
            t_array=[t_array,cur_t_array];
            p=[p,p1];
            t=[t,t1];
            max_len=cur_len;
        end
        A1 = max(p1(1:round(end*0.5)))-min(p1(1:round(end*0.5)));
        A2 = max(p1(round(end*0.5):end))-min(p1(round(end*0.5):end));
        A3 = max(t1)-min(t1);
        rratio = [rratio; -A2/A1];
        tratio = [tratio; A3/A1];
        count = count+1;
    end
end

writematrix([rratio,tratio], fullfile(save_data_path, ['sandbox/air_lossless1D' simtype '_ratio.txt']));

rratio = abs(rratio-rrt);
tratio = abs(tratio-trt);

figure('Position', [0,0,600,300]);
plot(ppwlist, rratio(:,1)/10,'-o','LineWidth',2);
hold on;
plot(ppwlist, tratio(:,1),'-o','LineWidth',2);
fontsize(20,"points");
set(gca, "LineWidth",2);
set(gca, "YScale", "log");
xlim([2,10]);
xlabel("Points per Wavelength");
ylabel("Error");
legend({"Reflection","Transmission"},'Location','northeast');
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossless1D_conv_' simtype '.png']));

Z2=5789*2240;
Z1=1500*1000;
rrt=(Z2-Z1)/(Z2+Z1);
trt=2*Z2/(Z2+Z1);

t_end = config.simulation.t_end;
t_array=[];
rratio=[];
tratio=[];
p=[];
t=[];
count=0;
max_len=0;

for i=10:-1:2
    mat_struct = load(fullfile(save_data_path, ['sandbox/glass_lossy_vert_refined' num2str(i) '.mat']));
    cur_t_array=mat_struct.kgrid.t_array;
    cur_t_array=transpose(cur_t_array);
    cur_len = length(cur_t_array);
    p1=mat_struct.p1;
    p1=transpose(p1);
    t1=mat_struct.p2;
    t1=transpose(t1);
    if count~=0
        t_array=[t_array,[cur_t_array;cur_t_array(end)*ones(max_len-cur_len,1)]];
        p=[p,[p1;p1(end)*ones(max_len-cur_len,1)]];
        t=[t,[t1;t1(end)*ones(max_len-cur_len,1)]];
    else
        t_array=[t_array,cur_t_array];
        p=[p,p1];
        t=[t,t1];
        max_len=cur_len;
    end
    A1 = max(p1(1:round(end*0.75)))-min(p1(1:round(end*0.75)));
    A2 = max(p1(round(end*0.75):end))-min(p1(round(end*0.75):end));
    A3 = max(t1)-min(t1);
    rratio = [rratio; A2/A1];
    tratio = [tratio; A3/A1];
    count = count+1;
end

rratio = abs(rratio-rrt);
tratio = abs(tratio-trt);

figure('Position', [0,0,600,300]);
for i=1:4:9
    plot(1e6*t_array(:,i), p(:,i)*1e-3, "LineWidth", 1);
    hold on;
end
fontsize(20,"points");
set(gca, "LineWidth",2);
xlim([10,20]);
ylim([-500,500]);
xlabel("Time [μs]");
ylabel("Pressure [kPa]");
legend("ppw=10", "ppw=6", "ppw=2");
saveas(gcf, fullfile(save_pic_path, ['sandbox/glass_lossy_vert_ppw.png']));

figure('Position', [0,0,600,300]);
for i=1:4:9
    plot(1e6*t_array(:,i), t(:,i)*1e-3, "LineWidth", 1);
    hold on;
end
fontsize(20,"points");
set(gca, "LineWidth",2);
xlim([10,20]);
ylim([-800,800]);
xlabel("Time [μs]");
ylabel("Pressure [kPa]");
legend("ppw=10", "ppw=6", "ppw=2");
saveas(gcf, fullfile(save_pic_path, ['sandbox/glass_lossy_vert_ppw_trans.png']));

figure('Position', [0,0,600,300]);
plot(10:-1:2, (rratio(:,1)),'-o','LineWidth',2);
hold on;
plot(10:-1:2, (tratio(:,1)),'-o','LineWidth',2);
fontsize(20,"points");
set(gca, "LineWidth",2);
xlim([1.8,10.2]);
set(gca, "YScale", "log");
xlabel("Points per Wave");
ylabel("Error");
legend("Reflection","Transmission");
saveas(gcf, fullfile(save_pic_path, ['sandbox/glass_lossy_vert_ppw_conv.png']));

Z2=1200*10;
Z1=1500*1000;
rrt=(Z2-Z1)/(Z2+Z1);
trt=2*Z2/(Z2+Z1);

t_end = config.simulation.t_end;
t_array=[];
rratio=[];
tratio=[];
p=[];
t=[];
count=0;
max_len=0;

for i=15:-1:2
    if exist(fullfile(save_data_path, ['sandbox/air_lossy_vert_refined' num2str(i) '.mat']))
        mat_struct = load(fullfile(save_data_path, ['sandbox/air_lossy_vert_refined' num2str(i) '.mat']));
        cur_t_array=mat_struct.kgrid.t_array;
        cur_t_array=transpose(cur_t_array);
        cur_len = length(cur_t_array);
        p1=mat_struct.p1;
        p1=transpose(p1);
        t1=mat_struct.p2;
        t1=transpose(t1);
        if count~=0
            t_array=[t_array,[cur_t_array;cur_t_array(end)*ones(max_len-cur_len,1)]];
            p=[p,[p1;p1(end)*ones(max_len-cur_len,1)]];
            t=[t,[t1;t1(end)*ones(max_len-cur_len,1)]];
        else
            t_array=[t_array,cur_t_array];
            p=[p,p1];
            t=[t,t1];
            max_len=cur_len;
        end
        A1 = max(p1(1:round(end*0.75)))-min(p1(1:round(end*0.75)));
        A2 = max(p1(round(end*0.75):end))-min(p1(round(end*0.75):end));
        A3 = max(t1)-min(t1);
        rratio = [rratio; -A2/A1];
        tratio = [tratio; A3/A1];
        count = count+1;
    end
end

figure('Position', [0,0,600,300]);
for i=1:5:11
    plot(1e6*t_array(:,i), p(:,i)*1e-3, "LineWidth", 1);
    hold on;
end
fontsize(20,"points");
set(gca, "LineWidth",2);
xlim([10,20]);
ylim([-500,500]);
xlabel("Time [μs]");
ylabel("Pressure [kPa]");
legend("ppw=15", "ppw=7", "ppw=2");
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossy_vert_ppw.png']));

figure('Position', [0,0,600,300]);
for i=1:5:11
    plot(1e6*t_array(:,i), t(:,i)*1e-3, "LineWidth", 1);
    hold on;
end
fontsize(20,"points");
set(gca, "LineWidth",2);
xlim([10,20]);
ylim([-40,40]);
xlabel("Time [μs]");
ylabel("Pressure [kPa]");
legend("ppw=15", "ppw=7", "ppw=2");
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossy_vert_ppw_trans.png']));

figure('Position', [0,0,600,300]);
plot([15,12,10:-1:2], rratio(:,1),'-o','LineWidth',2);
hold on;
plot([15,12,10:-1:2], tratio(:,1),'-o','LineWidth',2);
fontsize(20,"points");
set(gca, "LineWidth",2);
% set(gca, "YScale", "log");
xlim([1.8,15.2]);
xlabel("Points per Wave");
ylabel("Ratio");
legend("Reflection","Transmission");
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossy_vert_ppw_ppw.png']));

rratio = abs(rratio-rrt);
tratio = abs(tratio-trt);

figure('Position', [0,0,600,300]);
plot([15,12,10:-1:2], rratio(:,1),'-o','LineWidth',2);
hold on;
plot([15,12,10:-1:2], tratio(:,1),'-o','LineWidth',2);
fontsize(20,"points");
set(gca, "LineWidth",2);
set(gca, "YScale", "log");
xlim([1.8,15.2]);
xlabel("Points per Wave");
ylabel("Error");
legend("Reflection","Transmission");
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossy_vert_ppw_conv.png']));

Z2=343*1.2;
Z1=1500*1000;
rrt=(Z2-Z1)/(Z2+Z1);
trt=2*Z2/(Z2+Z1);

t_end = config.simulation.t_end;
t_array=[];
rratio=[];
tratio=[];
p=[];
t=[];
count=0;
max_len=0;

for i=5:-0.5:2
    if exist(fullfile(save_data_path, ['sandbox/air_lossy_vert_refined_true' num2str(i) '.mat']))
        mat_struct = load(fullfile(save_data_path, ['sandbox/air_lossy_vert_refined_true' num2str(i) '.mat']));
        cur_t_array=mat_struct.kgrid.t_array;
        cur_t_array=transpose(cur_t_array);
        cur_len = length(cur_t_array);
        p1=mat_struct.p1;
        p1=transpose(p1);
        t1=mat_struct.p2;
        t1=transpose(t1);
        if count~=0
            t_array=[t_array,[cur_t_array;cur_t_array(end)*ones(max_len-cur_len,1)]];
            p=[p,[p1;p1(end)*ones(max_len-cur_len,1)]];
            t=[t,[t1;t1(end)*ones(max_len-cur_len,1)]];
        else
            t_array=[t_array,cur_t_array];
            p=[p,p1];
            t=[t,t1];
            max_len=cur_len;
        end
        A1 = max(p1(1:round(end*0.75)))-min(p1(1:round(end*0.75)));
        A2 = max(p1(round(end*0.75):end))-min(p1(round(end*0.75):end));
        A3 = max(t1)-min(t1);
        rratio = [rratio; -A2/A1];
        tratio = [tratio; A3/A1];
        count = count+1;
    end
end

rratio = (rratio-rrt);
tratio = (tratio-trt);

figure('Position', [0,0,600,300]);
for i=[1,4,5]
    plot(1e6*t_array(:,i), p(:,i)*1e-3, "LineWidth", 1);
    hold on;
end
fontsize(20,"points");
set(gca, "LineWidth",2);
xlim([10,20]);
ylim([-500,500]);
xlabel("Time [μs]");
ylabel("Pressure [kPa]");
legend("ppw=5", "ppw=2.5", "ppw=2");
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossy_vert_ppw_true.png']));

figure('Position', [0,0,600,300]);
for i=[1,4,5]
    plot(1e6*t_array(:,i), t(:,i)*1e-3, "LineWidth", 1);
    hold on;
end
fontsize(20,"points");
set(gca, "LineWidth",2);
xlim([10,20]);
ylim([-40,40]);
xlabel("Time [μs]");
ylabel("Pressure [kPa]");
legend("ppw=5", "ppw=2.5", "ppw=2");
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossy_vert_ppw_trans_true.png']));

figure('Position', [0,0,600,300]);
plot([5,4,3,2.5,2], rratio(:,1)/10,'-o','LineWidth',2);
hold on;
plot([5,4,3,2.5,2], tratio(:,1),'-o','LineWidth',2);
fontsize(20,"points");
set(gca, "LineWidth",2);
set(gca, "YScale", "log");
xlim([1.8,5.2]);
xlabel("Points per Wave");
ylabel("Error");
legend({"Reflection","Transmission"},'Location','northwest');
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossy_vert_ppw_conv_true.png']));

Z2=343*10;
Z1=1500*1000;
rrt=(Z2-Z1)/(Z2+Z1);
trt=2*Z2/(Z2+Z1);

t_end = config.simulation.t_end;
t_array=[];
rratio=[];
tratio=[];
p=[];
t=[];
count=0;
max_len=0;

for i=6:-0.5:2
    if exist(fullfile(save_data_path, ['sandbox/air_lossy_vert_refined_truenew' num2str(i) '.mat']))
        mat_struct = load(fullfile(save_data_path, ['sandbox/air_lossy_vert_refined_truenew' num2str(i) '.mat']));
        cur_t_array=mat_struct.kgrid.t_array;
        cur_t_array=transpose(cur_t_array);
        cur_len = length(cur_t_array);
        p1=mat_struct.p1;
        p1=transpose(p1);
        t1=mat_struct.p2;
        t1=transpose(t1);
        if count~=0
            t_array=[t_array,[cur_t_array;cur_t_array(end)*ones(max_len-cur_len,1)]];
            p=[p,[p1;p1(end)*ones(max_len-cur_len,1)]];
            t=[t,[t1;t1(end)*ones(max_len-cur_len,1)]];
        else
            t_array=[t_array,cur_t_array];
            p=[p,p1];
            t=[t,t1];
            max_len=cur_len;
        end
        A1 = max(p1(1:round(end*0.75)))-min(p1(1:round(end*0.75)));
        A2 = max(p1(round(end*0.75):end))-min(p1(round(end*0.75):end));
        A3 = max(t1)-min(t1);
        rratio = [rratio; -A2/A1];
        tratio = [tratio; A3/A1];
        count = count+1;
    end
end

rratio = abs(rratio-rrt);
tratio = abs(tratio-trt);

figure('Position', [0,0,600,300]);
for i=[1,3,4]
    plot(1e6*t_array(:,i), p(:,i)*1e-3, "LineWidth", 1);
    hold on;
end
fontsize(20,"points");
set(gca, "LineWidth",2);
xlim([10,20]);
ylim([-500,500]);
xlabel("Time [μs]");
ylabel("Pressure [kPa]");
legend("ppw=6", "ppw=2.5", "ppw=2");
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossy_vert_ppw_truenew.png']));

figure('Position', [0,0,600,300]);
for i=[1,3,4]
    plot(1e6*t_array(:,i), t(:,i)*1e-3, "LineWidth", 1);
    hold on;
end
fontsize(20,"points");
set(gca, "LineWidth",2);
xlim([10,20]);
ylim([-40,40]);
xlabel("Time [μs]");
ylabel("Pressure [kPa]");
legend("ppw=6", "ppw=2.5", "ppw=2");
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossy_vert_ppw_trans_truenew.png']));

figure('Position', [0,0,600,300]);
plot([4,3,2.5,2], rratio(:,1),'-o','LineWidth',2);
hold on;
plot([4,3,2.5,2], tratio(:,1),'-o','LineWidth',2);
fontsize(20,"points");
set(gca, "LineWidth",2);
set(gca, "YScale", "log");
xlim([1.8,4.2]);
xlabel("Points per Wave");
ylabel("Error");
legend({"Reflection","Transmission"});
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossy_vert_ppw_conv_truenew.png']));

Z2=0*10;
Z1=1500*1000;
rrt=(Z2-Z1)/(Z2+Z1);
trt=2*Z2/(Z2+Z1);

t_end = config.simulation.t_end;
t_array=[];
rratio=[];
tratio=[];
p=[];
t=[];
count=0;
max_len=0;

for i=12:-1:2
    if exist(fullfile(save_data_path, ['sandbox/air_lossy_vert_refined_void' num2str(i) '.mat']))
        mat_struct = load(fullfile(save_data_path, ['sandbox/air_lossy_vert_refined_void' num2str(i) '.mat']));
        cur_t_array=mat_struct.kgrid.t_array;
        cur_t_array=transpose(cur_t_array);
        cur_len = length(cur_t_array);
        p1=mat_struct.p1;
        p1=transpose(p1);
        t1=mat_struct.p2;
        t1=transpose(t1);
        if count~=0
            t_array=[t_array,[cur_t_array;cur_t_array(end)*ones(max_len-cur_len,1)]];
            p=[p,[p1;p1(end)*ones(max_len-cur_len,1)]];
            t=[t,[t1;t1(end)*ones(max_len-cur_len,1)]];
        else
            t_array=[t_array,cur_t_array];
            p=[p,p1];
            t=[t,t1];
            max_len=cur_len;
        end
        A1 = max(p1(1:round(end*0.75)))-min(p1(1:round(end*0.75)));
        A2 = max(p1(round(end*0.75):end))-min(p1(round(end*0.75):end));
        A3 = max(t1)-min(t1);
        rratio = [rratio; -A2/A1];
        tratio = [tratio; A3/A1];
        count = count+1;
    end
end

rratio = abs(rratio-rrt);
tratio = abs(tratio-trt);

figure('Position', [0,0,600,300]);
for i=[1,4,7]
    plot(1e6*t_array(:,i), p(:,i)*1e-3, "LineWidth", 1);
    hold on;
end
fontsize(20,"points");
set(gca, "LineWidth",2);
xlim([10,20]);
ylim([-500,500]);
xlabel("Time [μs]");
ylabel("Pressure [kPa]");
legend("ppw=8", "ppw=5", "ppw=2");
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossy_vert_ppw_void.png']));

figure('Position', [0,0,600,300]);
plot([12,10:-1:2], rratio(:,1),'-o','LineWidth',2);
fontsize(20,"points");
set(gca, "LineWidth",2);
set(gca, "YScale", "log");
xlim([1.8,12.2]);
xlabel("Points per Wave");
ylabel("Error");
legend({"Reflection"});
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossy_vert_ppw_conv_void.png']));

Z2=0*10;
Z1=1500*1000;
rrt=(Z2-Z1)/(Z2+Z1);
trt=2*Z2/(Z2+Z1);

t_end = config.simulation.t_end;
t_array=[];
rratio=[];
tratio=[];
p=[];
t=[];
count=0;
max_len=0;

for i=8:-1:2
    if exist(fullfile(save_data_path, ['sandbox/air_lossy_vert_refined_voidrev' num2str(i) '.mat']))
        mat_struct = load(fullfile(save_data_path, ['sandbox/air_lossy_vert_refined_voidrev' num2str(i) '.mat']));
        cur_t_array=mat_struct.kgrid.t_array;
        cur_t_array=transpose(cur_t_array);
        cur_len = length(cur_t_array);
        p1=mat_struct.p1;
        p1=transpose(p1);
        t1=mat_struct.p2;
        t1=transpose(t1);
        if count~=0
            t_array=[t_array,[cur_t_array;cur_t_array(end)*ones(max_len-cur_len,1)]];
            p=[p,[p1;p1(end)*ones(max_len-cur_len,1)]];
            t=[t,[t1;t1(end)*ones(max_len-cur_len,1)]];
        else
            t_array=[t_array,cur_t_array];
            p=[p,p1];
            t=[t,t1];
            max_len=cur_len;
        end
        A1 = max(p1(1:round(end*0.75)))-min(p1(1:round(end*0.75)));
        A2 = max(p1(round(end*0.75):end))-min(p1(round(end*0.75):end));
        A3 = max(t1)-min(t1);
        rratio = [rratio; -A2/A1];
        tratio = [tratio; A3/A1];
        count = count+1;
    end
end

rratio = abs(rratio-rrt);
tratio = abs(tratio-trt);

figure('Position', [0,0,600,300]);
for i=[1,4,7]
    plot(1e6*t_array(:,i), p(:,i)*1e-3, "LineWidth", 1);
    hold on;
end
fontsize(20,"points");
set(gca, "LineWidth",2);
xlim([10,20]);
ylim([-500,500]);
xlabel("Time [μs]");
ylabel("Pressure [kPa]");
legend("ppw=8", "ppw=5", "ppw=2");
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossy_vert_ppw_voidrev.png']));

figure('Position', [0,0,600,300]);
plot(8:-1:2, rratio(:,1),'-o','LineWidth',2);
fontsize(20,"points");
set(gca, "LineWidth",2);
set(gca, "YScale", "log");
xlim([1.8,8.2]);
xlabel("Points per Wave");
ylabel("Error");
legend({"Reflection"});
saveas(gcf, fullfile(save_pic_path, ['sandbox/air_lossy_vert_ppw_conv_voidrev.png']));
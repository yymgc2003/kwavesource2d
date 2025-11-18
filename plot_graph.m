location_csv = "/mnt/yamaguchi/kwavesource_gl/location_seed/location1.csv";

location_df = readtable(location_csv);
location = table2array(location_df);

diameter = 13*transpose(location(:,4).^2.* location(:, 5)).^(1/3);

figure;
histogram(diameter, BinWidth=0.5);

saveas(gcf, fullfile("/mnt/yamaguchi/kwavesource_gl/location_seed", 'histogram.png'));

pos_5 = [];
pos_45 = [];
pos_34 = [];
pos_23 = [];
pos_12 = [];
pos_1 = [];
for i = 1:size(diameter, 2)
    if diameter(i) > 5
        pos_5 = [pos_5, (location(i,1)^2+location(i,2)^2)^0.5];
    elseif diameter(i) > 4
        pos_45 = [pos_45, (location(i,1)^2+location(i,2)^2)^0.5];
    elseif diameter(i) > 3
        pos_34 = [pos_34, (location(i,1)^2+location(i,2)^2)^0.5];
    elseif diameter(i) > 2
        pos_23 = [pos_23, (location(i,1)^2+location(i,2)^2)^0.5];
    elseif diameter(i) > 1
        pos_12 = [pos_12, (location(i,1)^2+location(i,2)^2)^0.5];
    else
        pos_1 = [pos_1, (location(i,1)^2+location(i,2)^2)^0.5];
    end
end

subplot(3,2,1);
histogram(pos_5, BinLimits=[0, 1], BinWidth=0.1);
xlabel('r/R');
ylabel('N');
title('d > 5 mm');

subplot(3,2,2);
histogram(pos_45, BinLimits=[0, 1], BinWidth=0.1);
xlabel('r/R');
ylabel('N');
title('4 mm < d < 5 mm');

subplot(3,2,3);
histogram(pos_34, BinLimits=[0, 1], BinWidth=0.1);
set(gca, 'YTick', [0:400:20000]);
xlabel('r/R');
ylabel('N');
title('3 mm < d < 4 mm');

subplot(3,2,4);
histogram(pos_23, BinLimits=[0, 1], BinWidth=0.1);
xlabel('r/R');
ylabel('N');
title('2 mm < d < 3 mm');

subplot(3,2,5);
histogram(pos_12, BinLimits=[0, 1], BinWidth=0.1);
xlabel('r/R');
ylabel('N');
title('1 mm < d < 2 mm');

subplot(3,2,6);
histogram(pos_1, BinLimits=[0, 1], BinWidth=0.1);
xlabel('r/R');
ylabel('N');
title('d < 1 mm');

saveas(gcf, fullfile("/mnt/yamaguchi/kwavesource_gl/location_seed", 'histogram_pos_big.png'));

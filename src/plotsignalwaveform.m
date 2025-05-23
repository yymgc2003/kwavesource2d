function plotsignalwaveform(kgrid, scan_line, save_path, filename)
% 信号波形をプロットし、画像として保存する関数
% kgrid: kWaveGridオブジェクト
% scan_line: プロットするデータ
% save_path: 保存先ディレクトリ
% filename: 保存するファイル名

figure;
plot(kgrid.t_array * 1e6, scan_line * 1e-6, 'b-');
xlabel('Time [\mus]');
ylabel('Pressure [MPa]');
ylim([-10 15]);
title('Signal from transducer transmit');
grid on;
saveas(gcf, fullfile(save_path, filename));
close(gcf);
end 
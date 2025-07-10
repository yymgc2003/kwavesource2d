% Run locationgen from kwavesource directory
% このスクリプトはkwavesourceフォルダで実行してください

% srcフォルダをパスに追加
addpath('src');

% 実行例：100個のサンプルを生成
number_samples = 100;
save_path = pwd;  % 現在のディレクトリに保存

% 関数を実行
generate_spaced_samples(number_samples, save_path);

fprintf('実行完了！\n'); 
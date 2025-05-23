# 2023-07-20 k-Waveシミュレーションデータ変換手順

## 1. ドキュメント作成者

松原貞徳　（smatsubara@fel.t.u-tokyo.ac.jp）

## 2. 改訂履歴

- 1.0:
  - 作成日時: 2023-07-20 10:30:00
  - 更新内容: 初版作成
- 1.1:
  - 作成日時: 2023-07-26 15:45:00
  - 更新内容: PicoScope BatchConvert機能を使用する方法に変更

## 3. このドキュメントの目次

- [1. ドキュメント作成者](#1-ドキュメント作成者)
- [2. 改訂履歴](#2-改訂履歴)
- [3. このドキュメントの目次](#3-このドキュメントの目次)
- [4. このドキュメントの目的・概要](#4-このドキュメントの目的概要)
- [5. 前提条件](#5-前提条件)
- [6. データ変換手順](#6-データ変換手順)
  - [6.1. 環境準備](#61-環境準備)
  - [6.2. 変換コマンドの実行](#62-変換コマンドの実行)
  - [6.3. 変換結果の確認](#63-変換結果の確認)
- [7. トラブルシューティング](#7-トラブルシューティング)
- [8. 参考情報](#8-参考情報)

## 4. このドキュメントの目的・概要

機械学習用データセットを作成するための前処理として、k-Waveシミュレーションで生成された`.psdata`形式のファイルを、Pythonで読み込み可能な`.mat`形式に変換する手順を説明します。本段階では、直接的に機械学習モデルの入力として使える形式のデータセットの作成ではなく、あくまでファイル形式の変換を目的としています。

変換にはPicoScope 7のBatchConvert機能を使用します。この機能により、複数の`.psdata`ファイルを一括で`.mat`形式に変換することが可能です。

## 5. 前提条件

- Windows環境であること
- PicoScope 7 T&M Stableがインストールされていること
  - インストールパス: `C:\Program Files\Pico Technology\PicoScope 7 T&M Stable\PicoScope.exe`
- Scandiumネットワークドライブにアクセス可能であること
  - ソースデータパス (例): `W:\2024\0726`
  - 出力先パス (例): `Z:\database\rawsignal`
- PowerShellが使用可能であること

## 6. データ変換手順

### 6.1. 環境準備

1. Scandiumネットワークドライブへの接続確認
   - エクスプローラーで `W:` ドライブと `Z:` ドライブにアクセスできることを確認
   - アクセスできない場合は、以下のコマンドでネットワークドライブをマウント
   ```powershell
   # 必要に応じてユーザー名とパスワードを指定
   New-PSDrive -Name "W" -PSProvider FileSystem -Root "\\scandium\data" -Persist
   New-PSDrive -Name "Z" -PSProvider FileSystem -Root "\\scandium\database" -Persist
   ```

2. PicoScopeのインストール確認
   - 以下のコマンドでPicoScopeの実行ファイルが存在するか確認
   ```powershell
   Test-Path "C:\Program Files\Pico Technology\PicoScope 7 T&M Stable\PicoScope.exe"
   ```

### 6.2. 変換コマンドの実行

1. 単一フォルダの変換
   - PowerShellを管理者権限で起動
   - 以下のコマンドを実行
   ```powershell
   Start-Process -FilePath "C:\Program Files\Pico Technology\PicoScope 7 T&M Stable\PicoScope.exe" -ArgumentList "BatchConvert", "W:\2024\0726", "Z:\database\rawsignal", ".mat" -Wait -NoNewWindow
   ```

2. 複数フォルダの一括変換
   - 複数のフォルダを処理する場合は、以下のようなスクリプトを作成して実行
   ```powershell
   # batch_convert.ps1
   $sourceFolders = @(
       "W:\2024\0726",
       "W:\2024\0727",
       "W:\2024\0728"
   )
   $outputFolder = "Z:\database\rawsignal"
   $outputFormat = ".mat"
   
   foreach ($folder in $sourceFolders) {
       Write-Host "Converting files in $folder to $outputFormat format..."
       Start-Process -FilePath "C:\Program Files\Pico Technology\PicoScope 7 T&M Stable\PicoScope.exe" -ArgumentList "BatchConvert", $folder, $outputFolder, $outputFormat -Wait -NoNewWindow
       Write-Host "Conversion completed for $folder"
   }
   
   Write-Host "All conversions completed successfully!"
   ```
   
   - スクリプトの実行
   ```powershell
   .\batch_convert.ps1
   ```

3. スケジュール実行
   - 定期的に変換を実行したい場合は、Windowsタスクスケジューラを使用
   ```powershell
   # 毎日午前2時に実行するタスクを作成
   $action = New-ScheduledTaskAction -Execute "PowerShell.exe" -Argument "-File C:\Scripts\batch_convert.ps1"
   $trigger = New-ScheduledTaskTrigger -Daily -At 2am
   Register-ScheduledTask -TaskName "PicoScope Data Conversion" -Action $action -Trigger $trigger -RunLevel Highest
   ```

### 6.3. 変換結果の確認

1. 出力ファイルの確認
   ```powershell
   # 出力ファイル数の確認
   (Get-ChildItem -Path "Z:\database\rawsignal" -Filter "*.mat" -Recurse).Count
   ```

2. MATLABで変換後のファイルを読み込んで確認
   ```matlab
   % load_and_verify.m
   clear all;
   close all;
   
   % 変換後のMATファイルを読み込む
   matFile = 'Z:\database\rawsignal\example_data.mat';
   data = load(matFile);
   
   % データ構造の確認
   disp('データ構造:');
   disp(fieldnames(data));
   
   % データの可視化
   if isfield(data, 'sensor_data')
       figure;
       plot(data.sensor_data);
       title('Sensor Data');
       xlabel('Time');
       ylabel('Amplitude');
       saveas(gcf, 'sensor_data_plot.png');
   end
   ```

## 7. トラブルシューティング

- **エラー：PicoScopeが起動しない**
  - 解決策：PicoScopeのインストールパスが正しいことを確認
  - インストールパスが異なる場合はコマンドの `-FilePath` パラメータを適宜修正

- **エラー：ネットワークドライブにアクセスできない**
  - 解決策：ネットワーク接続とアクセス権限を確認
  - VPNを使用している場合は、VPN接続が安定していることを確認

- **エラー：変換が途中で停止する**
  - 解決策：大量のファイルを変換する場合は、小さなバッチに分割して実行
  - PowerShellスクリプトにリトライロジックを追加
  ```powershell
  $maxRetries = 3
  $retryCount = 0
  $success = $false
  
  while (-not $success -and $retryCount -lt $maxRetries) {
      try {
          Start-Process -FilePath "C:\Program Files\Pico Technology\PicoScope 7 T&M Stable\PicoScope.exe" -ArgumentList "BatchConvert", $folder, $outputFolder, $outputFormat -Wait -NoNewWindow -ErrorAction Stop
          $success = $true
      }
      catch {
          $retryCount++
          Write-Host "Error occurred. Retrying... Attempt $retryCount of $maxRetries"
          Start-Sleep -Seconds 10
      }
  }
  ```

- **エラー：変換後のファイルサイズが0KB**
  - 解決策：ソースファイルが正常であることを確認
  - PicoScopeのバージョンが最新であることを確認

- **エラー：変換元のファイルが破損**
  - 解決策：無視して報告
  - 20240731は破損している。

- **エラー：出力ファイルが破損している**
  - 解決策：ディスク容量に十分な空きがあることを確認
  - アンチウイルスソフトによるファイルアクセスの干渉がないか確認

## 8. 参考情報

- [PicoScopeコマンドラインオプション](https://www.picotech.com/download/manuals/picoscope-7-command-line-manual.pdf)
- [PicoTechnologyサポートページ](https://www.picotech.com/support)
- [PowerShellスクリプト作成ガイド](https://docs.microsoft.com/ja-jp/powershell/scripting/overview)
- [MATLAB `.mat`ファイル形式仕様](https://www.mathworks.com/help/matlab/import_export/mat-file-versions.html) 
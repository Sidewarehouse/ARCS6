% MultiPlotスクリプトのサンプルコード
% 2024/06/25 Yokokura, Yuki
clc;
clear;

% CSVファイルから読み込み
FileName = '../001_空の基本コード/DATA.csv';	% CSVファイル名設定
import MultiPlot.LoadCsvFile;		% LoadCsvFile関数インポート
[t, x1, x2, x3, x4, x5, x6, x7, x8, x9] = LoadCsvFile(FileName);	% 変数値読み込み

% グラフ時間軸オフセット＆スケーリング
%t = t - 0.05;	% 開始時刻のオフセット
%t = t*1e3;		% [s] → [ms] のスケーリング

% グラフ描画
figure(1);
	% MultiPlot全体の設定
	Graph1 = MultiPlot(gcf);			% MultiPlot生成
	Graph1.FigurePosition(55, 55);		% Figure位置の設定(左[px], 下[px])
	Graph1.FigureSize(900, 800);		% Figureサイズの設定(幅[px], 高さ[px])
	Graph1.FigureMargin(130, 70, 20);	% Figure余白の設定(左側[px], 下側[px], 右側[px])
	Graph1.Font('Times New Roman', 18);	% フォントの設定(フォント名, フォントサイズ)
	Graph1.NumOfPlanes(2);				% プロット平面の段数
	%Graph1.XaxisRange(-1, 1, 5);		% X軸範囲の設定(最小値, グリッド間隔, 最大値) ←コメントアウトで自動モード
	Graph1.XaxisLabel('Time [s]');		% X軸ラベル名
	% プロット平面1段目
	Graph1.SelectPlane(1);							% プロット段選択
	Graph1.StairsPlot(t, x1, 'Black', 'Thin');		% 階段プロット
	Graph1.StairsPlot(t, x2, 'Red',   'Normal');	% 線の色は Black/Gray/Red/Green/Blue の5種のみ
	Graph1.StairsPlot(t, x3, 'Green', 'Normal');	% 線の太さは Thin/Normal/Heavy の3種のみ
	Graph1.StairsPlot(t, x4, 'Blue',  'Heavy');		% 重ねられるプロットの数は無制限
	%Graph1.ManualGrid(-2, 0.5, 2);					% Y軸範囲の設定(最小値, グリッド間隔, 最大値) ←コメントアウトで自動モード
	Graph1.Label({'(a)','Var. A','[-]'});			% 縦軸ラベルの設定
	Graph1.Legend({'x_1','x_2','x_3','x_4'}, 'NorthEast', 'Vertical');	% 凡例の設定
	% プロット平面2段目
	Graph1.SelectPlane(2);							% プロット段選択
	Graph1.LinePlot(t, x5, 'Black', 'Thin');		% 線形プロット
	Graph1.LinePlot(t, x6, 'Red',   'Normal');		% 線の色は Black/Gray/Red/Green/Blue の5種のみ
	Graph1.LinePlot(t, x7, 'Green', 'Normal');		% 線の太さは Thin/Normal/Heavy の3種のみ
	Graph1.LinePlot(t, x8, 'Blue',  'Heavy');		% 重ねられるプロットの数は無制限
	Graph1.LinePlot(t, x9, 'Blue',  'Heavy');		% 重ねられるプロットの数は無制限
	%Graph1.ManualGrid(-2, 0.5, 2);					% Y軸範囲の設定(最小値, グリッド間隔, 最大値) ←コメントアウトで自動モード
	Graph1.Label({'(b)','Var. B','[-]'});			% 縦軸ラベルの設定
	Graph1.Legend({'x_5','x_6','x_7','x_8','x_9'}, 'NorthEast', 'Vertical');	% 凡例の設定
	% 画像生成
	Graph1.SavePNGandEPS('DATA.eps');		% PNG画像とEPSファイルを生成
	Graph1.SavePNGandPDF('DATA.pdf');		% PNG画像とPDFファイルを生成


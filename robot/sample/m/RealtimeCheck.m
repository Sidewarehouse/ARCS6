% リアルタイム性の測定データプロット
% 2025/11/25 Yokokura, Yuki
clc;
clear;

% CSVファイルから読み込み
FileName = '../011_リアルタイム性の測定/DATA.csv';	% CSVファイル名設定
import MultiPlot.LoadCsvFile;			% LoadCsvFile関数インポート
[t, Tcmp, Tact] = LoadCsvFile(FileName);% 変数値読み込み

% グラフ時間軸オフセット＆スケーリング
%t = t - 0.05;	% 開始時刻のオフセット
%t = t*1e3;		% [s] → [ms] のスケーリング

% 消費時間と制御周期の時系列描画
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
	Graph1.StairsPlot(t, Tcmp*1e6, 'Black', 'Thin');% 計算消費時間のプロット
	%Graph1.ManualGrid(0, 25, 100);					% Y軸範囲の設定(最小値, グリッド間隔, 最大値) ←コメントアウトで自動モード
	Graph1.Label({'(a)','Consumption Time','Tcmp [\mus]'});			% 縦軸ラベルの設定
	% プロット平面2段目
	Graph1.SelectPlane(2);							% プロット段選択
	Graph1.StairsPlot(t, Tact*1e6, 'Black', 'Thin');% 実際の制御周期時間のプロット
	%Graph1.ManualGrid(0, 50, 200);					% Y軸範囲の設定(最小値, グリッド間隔, 最大値) ←コメントアウトで自動モード
	Graph1.Label({'(b)','Actual Periodic Time','Tact [\mus]'});		% 縦軸ラベルの設定
	% 画像生成
	%Graph1.SavePNGandEPS('DATA.eps');		% PNG画像とEPSファイルを生成
	%Graph1.SavePNGandPDF('DATA.pdf');		% PNG画像とPDFファイルを生成

% 制御周期ヒストグラムの描画
figure(2);
	set(gcf,'Color','white');
	h = histogram(Tact*1e6);
		set(h, 'BinWidth', 0.1);
		set(h, 'FaceColor', 'black');
		set(h, 'EdgeColor', 'black');
	set(gca,'yscale','log');
	set(gca,'FontSize',18, 'FontName','Times New Roman'), ;
	xlabel('Actual Periodic Time Tact [\mus]','FontSize',18);
	ylabel({'(c)','Number of Control Loops [-]'},'FontSize',18);
	axis([90 inf 0.1 inf]);
	grid on;


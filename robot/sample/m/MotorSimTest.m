% モータシミュレータテストコード
% 2026/05/09 Yokokura, Yuki
clc;
clear;

% CSVファイルから読み込み
FileName = '../051_モータシミュレータ/DATA.csv';	% CSVファイル名設定
import MultiPlot.LoadCsvFile;		% LoadCsvFile関数インポート
[t, iq1, wm1, th1, x4, x5, x6, x7, x8, x9] = LoadCsvFile(FileName);	% 変数値読み込み

% グラフ時間軸オフセット＆スケーリング
%t = t - 0.05;	% 開始時刻のオフセット
%t = t*1e3;		% [s] → [ms] のスケーリング

% MATLAB側シミュレーション
Ts = 100e-6;	% [s] 制御周期
Kt = 0.5;		% [Nm/A] トルク定数
Jm = 1.00e-3;	% [kgm^2] モータ側慣性
Dm = 1.00e-3;	% [Nm s/rad] モータ側粘性
s = tf('s');
P = Kt/(Jm*s + Dm);
ts = linspace(0, 10, 1/Ts*10);
tslen = length(ts);
iqs(1:tslen) = 0;
iqs(1/Ts*1:tslen) = 0.3;
wms = lsim(P, iqs, ts);
ths = lsim(1/s, wms, ts);

% グラフ描画
figure(1);
	% MultiPlot全体の設定
	Graph1 = MultiPlot(gcf);			% MultiPlot生成
	Graph1.FigurePosition(55, 55);		% Figure位置の設定(左[px], 下[px])
	Graph1.FigureSize(900, 800);		% Figureサイズの設定(幅[px], 高さ[px])
	Graph1.FigureMargin(160, 70, 20);	% Figure余白の設定(左側[px], 下側[px], 右側[px])
	Graph1.Font('Times New Roman', 18);	% フォントの設定(フォント名, フォントサイズ)
	Graph1.NumOfPlanes(3);				% プロット平面の段数
	Graph1.XaxisRange(0, 1, 10);		% X軸範囲の設定(最小値, グリッド間隔, 最大値) ←コメントアウトで自動モード
	Graph1.XaxisLabel('Time [s]');		% X軸ラベル名
	% プロット平面1段目
	Graph1.SelectPlane(1);							% プロット段選択
	Graph1.StairsPlot(ts, iqs, 'GemBlue', 'Heavy');	% プロット
	Graph1.StairsPlot(t, iq1, 'GemRed', 'Thin');	% プロット
	%Graph1.ManualGrid(-2, 0.5, 2);					% Y軸範囲の設定(最小値, グリッド間隔, 最大値) ←コメントアウトで自動モード
	Graph1.Label({'(a)','Current i_q','[A]'});			% 縦軸ラベルの設定
	Graph1.Legend({'MATLAB', 'ARCS'}, 'SouthEast', 'Vertical');	% 凡例の設定
	% プロット平面2段目
	Graph1.SelectPlane(2);							% プロット段選択
	Graph1.StairsPlot(ts, wms, 'GemBlue', 'Heavy');	% プロット
	Graph1.StairsPlot(t, wm1, 'GemRed', 'Thin');	% プロット
	%Graph1.ManualGrid(-20, 2, 20);					% Y軸範囲の設定(最小値, グリッド間隔, 最大値) ←コメントアウトで自動モード
	Graph1.Label({'(b)','Motor Speed \omega_m','[rad/s]'});			% 縦軸ラベルの設定
	Graph1.Legend({'MATLAB', 'ARCS'}, 'SouthEast', 'Vertical');	% 凡例の設定
	% プロット平面3段目
	Graph1.SelectPlane(3);							% プロット段選択
	Graph1.StairsPlot(ts, ths, 'GemBlue', 'Heavy');	% プロット
	Graph1.StairsPlot(t, th1, 'GemRed', 'Thin');	% プロット
	%Graph1.ManualGrid(-20, 2, 20);					% Y軸範囲の設定(最小値, グリッド間隔, 最大値) ←コメントアウトで自動モード
	Graph1.Label({'(c)','Motor Position \theta_m','[rad]'});		% 縦軸ラベルの設定
	Graph1.Legend({'MATLAB', 'ARCS'}, 'SouthEast', 'Vertical');	% 凡例の設定
	% 画像生成
	%Graph1.SavePNGandEPS('DATA.eps');		% PNG画像とEPSファイルを生成
	%Graph1.SavePNGandPDF('DATA.pdf');		% PNG画像とPDFファイルを生成


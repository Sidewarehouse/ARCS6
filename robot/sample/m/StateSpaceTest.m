% ArcsControl::StateSpaceクラスのテスト用スクリプト
% 2025/01/17 Yokokura, Yuki
clc;
clear;

% CSVファイル名設定
FileName = '../042_任意の状態空間モデル/DATA.csv';

% テスト対象の状態空間モデル
A = [
	 0,  1 ;
	-5, -2
];
b = [
	0 ;
	3
];
C = [
	1, 0
];
D = [
	1
];
System1 = ss(A,b,C,D)

% CSVファイルから変数値読み込み
CsvData  = csvread(FileName);
t = CsvData(:,1);
u = CsvData(:,2);
y = CsvData(:,3);
clear CsvData;
tlen = length(t);

% 真値の計算
Ts = t(tlen)/(tlen - 1);
ti = 0:Ts:t(tlen);
ui = SquareWave(0.5, pi, ti, 1);
yi = lsim(System1, ui, ti);

% グラフ描画
figure(1);
clf;
set(gcf,'PaperPositionMode','manual');
set(gcf,'color',[1 1 1]);
subplot(2,1,1);
	h=stairs(t, u, 'k:');
		set(h,'linewidth',4);
	hold on;
	h=stairs(t, y, 'r');
		set(h,'linewidth',3);
	h=stairs(ti, ui, 'g:');
		set(h,'linewidth',2);
	h=stairs(ti, yi, 'k');
		set(h,'linewidth',1);
	hold off;
	xlabel('Time [s]','FontSize',12);
	ylabel('Input and Output [-]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	%axis([0 10 -inf inf]);
	legend('ARCS Input','ARCS Output','MATLAB Input','MATLAB Output','Location','SouthEast','Orientation','Vertical');
subplot(2,1,2);
	h=stairs(t, ti, 'k');
		set(h,'linewidth',1);
	xlabel('Ideal Time [s]','FontSize',12);
	ylabel('Actual Time [s]','FontSize',12);
	grid on;


% EPSファイル生成(ローカルで実行のこと)
% print(gcf,'-depsc2','-tiff',strcat(FileName,'.eps'));

% 方形波生成関数
function y = SquareWave(freq, phase, time, starttime)
	tlen = length(time);
	r(1:tlen) = 0;
	y(1:tlen) = 0;
	
	% 開始時間だけシフトされた正弦波を生成
	for i = 1:tlen
		if time(i) < starttime
			r(i) = 0;			% 開始時刻より前のときはゼロ出力
		else
			r(i) = sin(2*pi*freq*(time(i) - starttime) + phase);	% 開始時刻以降は方形波出力
			% 方形波に変換
			if 0 <= r(i)
				y(i) =  1;
			else
				y(i) = -1;
			end
		end
	end
end



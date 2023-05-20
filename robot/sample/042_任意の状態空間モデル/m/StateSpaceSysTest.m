% StateSpaceSystemクラスのテスト用スクリプト
% 2022/02/20 Yokokura, Yuki
clc;
clear;

% CSVファイル名設定
FileName = '../DATA.csv';

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
Ts = t(tlen)/(tlen-1);
ti = 0:Ts:t(tlen);
yi = lsim(System1, u, ti);

% グラフ描画
figure(1);
clf;
set(gcf,'PaperPositionMode','manual');
set(gcf,'color',[1 1 1]);
h=stairs(t, u, 'k:');
	set(h,'linewidth',2);
hold on;
h=stairs(t, y, 'r');
	set(h,'linewidth',3);
h=stairs(t, yi, 'k');
	set(h,'linewidth',1);
hold off;
xlabel('Time [s]','FontSize',12);
ylabel('Input and Output [-]','FontSize',12);
set(gca,'FontSize',12);
grid on;
%axis([0 10 -inf inf]);
legend('Input','Output','Ideal','Location','SouthEast','Orientation','Vertical');

% EPSファイル生成(ローカルで実行のこと)
% print(gcf,'-depsc2','-tiff',strcat(FileName,'.eps'));

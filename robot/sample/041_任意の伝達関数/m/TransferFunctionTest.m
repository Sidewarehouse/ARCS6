% TransferFunctioクラスのテスト用スクリプト
% 2022/11/11 Yokokura, Yuki
clc;
clear;

% CSVファイル名設定
FileName = '../DATA.csv';

% テスト対象の伝達関数
s = tf('s');
G1 = (5*s^2 + 6*s + 7)/(2*s^3 + 3*s^2 + 4*s + 5)
G2 = (4*s^3 + 5*s^2 + 6*s + 7)/(2*s^3 + 3*s^2 + 4*s + 5)
w = 30;
z = 0.2;
G3 = w^2/(s^2 + 2*z*w*s + w^2)

% CSVファイルから変数値読み込み
CsvData  = csvread(FileName);
t = CsvData(:,1);
u = CsvData(:,2);
y1 = CsvData(:,3);
y2 = CsvData(:,4);
y3 = CsvData(:,5);
clear CsvData;
tlen = length(t);

% 真値の計算
Ts = t(tlen)/(tlen-1);
ti = 0:Ts:t(tlen);
y1i = lsim(G1, u, ti);
y2i = lsim(G2, u, ti);
y3i = lsim(G3, u, ti);

% グラフ描画
figure(1);
	clf;
	set(gcf,'PaperPositionMode','manual');
	set(gcf,'color',[1 1 1]);
subplot(4,1,1);
	h=plot(t, u, 'k');
		set(h,'linewidth',2);
	xlabel('Time [s]','FontSize',12);
	ylabel('Input u [-]','FontSize',12);
	set(gca,'FontSize',12);
subplot(4,1,2);
	hold on;
	h=plot(t, y1i, 'r');
		set(h,'linewidth',3);
	h=plot(t, y1, 'k');
		set(h,'linewidth',1);
	hold off;
	xlabel('Time [s]','FontSize',12);
	ylabel('Output y1 [-]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	%axis([0 10 -inf inf]);
	legend('MATLAB','ARCS6','Location','SouthEast','Orientation','Vertical');
subplot(4,1,3);
	hold on;
	h=plot(t, y2i, 'r');
		set(h,'linewidth',3);
	h=plot(t, y2, 'k');
		set(h,'linewidth',1);
	hold off;
	xlabel('Time [s]','FontSize',12);
	ylabel('Output y2 [-]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	%axis([0 10 -inf inf]);
	legend('MATLAB','ARCS6','Location','SouthEast','Orientation','Vertical');
subplot(4,1,4);
	hold on;
	h=plot(t, y3i, 'r');
		set(h,'linewidth',3);
	h=plot(t, y3, 'k');
		set(h,'linewidth',1);
	hold off;
	xlabel('Time [s]','FontSize',12);
	ylabel('Output y3 [-]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	%axis([0 10 -inf inf]);
	legend('MATLAB','ARCS6','Location','SouthEast','Orientation','Vertical');

% EPSファイル生成(ローカルで実行のこと)
% print(gcf,'-depsc2','-tiff',strcat(FileName,'.eps'));

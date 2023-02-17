% ARCSが出力するDATA.csvをMATLABに読み込むスクリプトの一例
% 2020/04/03 Yokokura, Yuki
clc;
clear;

% CSVファイル名設定
FileName = '../DATA.csv';

% CSVファイルから変数値読み込み
CsvData  = csvread(FileName);
t = CsvData(:,1);
A1 = CsvData(:,2);
A2 = CsvData(:,3);
B1 = CsvData(:,4);
B2 = CsvData(:,5);
C1 = CsvData(:,6);
C2 = CsvData(:,7);
D1 = CsvData(:,8);
D2 = CsvData(:,9);
clear CsvData;
tlen = length(t);

% ------- データを間引く場合
%{
RedRate = 20;	% 間引く要素数
t = t(1:RedRate:tlen);
A = A(1:RedRate:tlen);
B = B(1:RedRate:tlen);
C = C(1:RedRate:tlen);
tlen = length(t);
%}

% グラフ描画
figure(1);
clf;
set(gcf,'PaperPositionMode','manual');
set(gcf,'color',[1 1 1]);
%set(gcf,'Position',[100 100 800 900]);
subplot(4,1,1);
	h=plot(t, A1, 'k');
		set(h,'linewidth',4);
	hold on;
	h=plot(t, A2, 'r');
		set(h,'linewidth',2);
	hold off;
	xlabel({'Time [s]','(a)'},'FontSize',12);
	ylabel(' Variable A [-]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	%axis([0 10 -inf inf]);
	%set(gca,'YTickMode','manual');
	%set(gca,'XTick', 0:0.1:0.5);
	%set(gca,'YTick', 0:0.1:0.5);
	legend('A1','A2','Location','SouthEast','Orientation','Vertical');
	%legend boxoff;
subplot(4,1,2);
	h=plot(t, B1, 'k');
		set(h,'linewidth',4);
	hold on;
	h=plot(t, B2, 'r');
		set(h,'linewidth',2);
	hold off;
	xlabel({'Time [s]','(b)'},'FontSize',12);
	ylabel(' Variable B [-]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	%axis([0 10 -inf inf]);
	%set(gca,'YTickMode','manual');
	%set(gca,'XTick', 0:0.1:0.5);
	%set(gca,'YTick', 0:0.1:0.5);
	legend('B1','B2','Location','SouthEast','Orientation','Vertical');
	%legend boxoff;
subplot(4,1,3);
	h=plot(t, C1, 'k');
		set(h,'linewidth',4);
	hold on;
	h=plot(t, C2, 'r');
		set(h,'linewidth',2);
	hold off;
	xlabel({'Time [s]','(c)'},'FontSize',12);
	ylabel(' Variable C [-]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	%axis([0 10 -inf inf]);
	%set(gca,'YTickMode','manual');
	%set(gca,'XTick', 0:0.1:0.5);
	%set(gca,'YTick', 0:0.1:0.5);
	legend('C1','C2','Location','SouthEast','Orientation','Vertical');
	%legend boxoff;
subplot(4,1,4);
	h=plot(t, D1, 'k');
		set(h,'linewidth',4);
	hold on;
	h=plot(t, D2, 'r');
		set(h,'linewidth',2);
	hold off;
	xlabel({'Time [s]','(d)'},'FontSize',12);
	ylabel(' Variable D [-]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	%axis([0 10 -inf inf]);
	%set(gca,'YTickMode','manual');
	%set(gca,'XTick', 0:0.1:0.5);
	%set(gca,'YTick', 0:0.1:0.5);
	legend('D1','D2','Location','SouthEast','Orientation','Vertical');
	%legend boxoff;

% EPSファイル生成(ローカルで実行のこと)
% print(gcf,'-depsc2','-tiff',strcat(FileName,'.eps'));

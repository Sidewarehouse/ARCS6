% ARCSが出力するDATA.csvをMATLABに読み込むスクリプトの一例
% 2022/03/08 Yokokura, Yuki
clc;
clear;

% CSVファイル名設定
FileName = '../DATA.csv';

% CSVファイルから変数値読み込み
CsvData  = csvread(FileName);
t = CsvData(:,1);
ps = CsvData(:,6);
fG = CsvData(:,7);
pG = CsvData(:,9);
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

psref(1:tlen) = 0;

% グラフ描画
figure(1);
clf;
set(gcf,'PaperPositionMode','manual');
set(gcf,'color',[1 1 1]);
%set(gcf,'Position',[100 100 800 900]);
subplot(3,1,1);
	h=plot(t, psref, 'k');
		set(h,'linewidth',4);
	hold on;
	h=plot(t, ps, 'r');
		set(h,'linewidth',2);
	hold off;
	xlabel({'Time [s]','(a)'},'FontSize',12);
	ylabel('Momentum [Nms]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	%axis([0 10 -inf inf]);
	%set(gca,'YTickMode','manual');
	%set(gca,'XTick', 0:0.1:0.5);
	%set(gca,'YTick', 0:0.1:0.5);
	legend('Reference','Response','Location','SouthEast','Orientation','Vertical');
	%legend boxoff;
subplot(3,1,2);
	h=plot(t, fG, 'k');
		set(h,'linewidth',4);
	xlabel({'Time [s]','(b)'},'FontSize',12);
	ylabel('Force [N]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	%axis([0 10 -inf inf]);
	%set(gca,'YTickMode','manual');
	%set(gca,'XTick', 0:0.1:0.5);
	%set(gca,'YTick', 0:0.1:0.5);
	%legend boxoff;
subplot(3,1,3);
	h=plot(t, pG, 'k');
		set(h,'linewidth',4);
	xlabel({'Time [s]','(c)'},'FontSize',12);
	ylabel('Position [m]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	%axis([0 10 -inf inf]);
	%set(gca,'YTickMode','manual');
	%set(gca,'XTick', 0:0.1:0.5);
	%set(gca,'YTick', 0:0.1:0.5);
	%legend boxoff;

% EPSファイル生成(ローカルで実行のこと)
% print(gcf,'-depsc2','-tiff',strcat(FileName,'.eps'));


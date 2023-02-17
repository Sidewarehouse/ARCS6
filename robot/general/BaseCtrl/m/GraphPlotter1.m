% ARCSが出力するDATA.csvをMATLABに読み込むスクリプトの一例
% 2020/04/06 Yokokura, Yuki
clc;
clear;

% CSVファイル名設定
FileName = '../DATA.csv';

% CSVファイルから変数値読み込み
CsvData  = csvread(FileName);
t = CsvData(:,1);
A = CsvData(:,2);
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
h=plot(t, A, 'k');
	set(h,'linewidth',2);
xlabel('Time [s]','FontSize',12);
ylabel(' Variable A [-]','FontSize',12);
set(gca,'FontSize',12);
grid on;
%axis([0 10 -inf inf]);
%set(gca,'YTickMode','manual');
%set(gca,'XTick', 0:0.1:0.5);
%set(gca,'YTick', 0:0.1:0.5);
%legend('Line 1','Line 2','Location','SouthEast','Orientation','Vertical');
%legend boxoff;

% EPSファイル生成(ローカルで実行のこと)
% print(gcf,'-depsc2','-tiff',strcat(FileName,'.eps'));

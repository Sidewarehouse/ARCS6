% ARCSが出力するDATA.csvをMATLABに読み込むスクリプトの一例
% 重力項同定用
% 2022/02/17 Yokokura, Yuki
clc;
clear;

% CSVファイル名設定
FileName = '../DATA.csv';

% CSVファイルから変数値読み込み
CsvData  = csvread(FileName);
t = CsvData(:,1);
thmpro = CsvData(:,2);
thm = CsvData(:,3);
taus = CsvData(:,5);
clear CsvData;
tlen = length(t);

% グラフ描画
figure(1);
clf;
set(gcf,'PaperPositionMode','manual');
set(gcf,'color',[1 1 1]);
subplot(2,1,1);
	h=plot(t, thmpro, 'k');
		set(h,'linewidth',4);
	hold on;
	h=plot(t, thm, 'r');
		set(h,'linewidth',2);
	hold off;
	xlabel({'Time [s]','(a)'},'FontSize',12);
	ylabel('Motor-side position [rad]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	%axis([0 10 -inf inf]);
	%set(gca,'YTickMode','manual');
	%set(gca,'XTick', 0:0.1:0.5);
	%set(gca,'YTick', 0:0.1:0.5);
	legend('Profile','Actual','Location','SouthEast','Orientation','Vertical');
	%legend boxoff;
subplot(2,1,2);
	h=plot(t, taus, 'k');
		set(h,'linewidth',4);
	xlabel({'Time [s]','(b)'},'FontSize',12);
	ylabel('Torsion torque [Nm]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	%axis([0 10 -inf inf]);
	%set(gca,'YTickMode','manual');
	%set(gca,'XTick', 0:0.1:0.5);
	%set(gca,'YTick', 0:0.1:0.5);
	%legend boxoff;

% EPSファイル生成(ローカルで実行のこと)
% print(gcf,'-depsc2','-tiff',strcat(FileName,'.eps'));

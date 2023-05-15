% ARCS���o�͂���DATA.csv��MATLAB�ɓǂݍ��ރX�N���v�g�̈��
% 2020/04/10 Yuki YOKOKURA
clc;
clear;

% CSV�t�@�C�����ݒ�
FileName = '../DATA.csv';

% CSV�t�@�C������ϐ��l�ǂݍ���
CsvData  = csvread(FileName);
t = CsvData(:,1);
Var1 = CsvData(:,2);
Var2 = CsvData(:,3);
Var3 = CsvData(:,4);
Var4 = CsvData(:,5);
Var5 = CsvData(:,6);
clear CsvData;
tlen = length(t);

% ------- �f�[�^���Ԉ����ꍇ
%{
RedRate = 20;	% �Ԉ����v�f��
t = t(1:RedRate:tlen);
A = A(1:RedRate:tlen);
B = B(1:RedRate:tlen);
C = C(1:RedRate:tlen);
tlen = length(t);
%}

% �O���t�`��
figure(1);
clf;
set(gcf,'PaperPositionMode','manual');
set(gcf,'color',[1 1 1]);
%set(gcf,'Position',[100 100 800 900]);
subplot(5,1,1);
	h=plot(t, Var1, 'k');
		set(h,'linewidth',2);
	xlabel({'Time [s]','(a)'},'FontSize',12);
	ylabel('Variable 1 [-]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	%axis([0 10 -inf inf]);
	%set(gca,'YTickMode','manual');
	%set(gca,'XTick', 0:0.1:0.5);
	%set(gca,'YTick', 0:0.1:0.5);
	%legend('Line 1','Line 2','Location','SouthEast','Orientation','Vertical');
	%legend boxoff;
subplot(5,1,2);
	h=plot(t, Var2, 'k');
		set(h,'linewidth',2);
	xlabel({'Time [s]','(b)'},'FontSize',12);
	ylabel('Variable 2 [-]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	%axis([0 10 -inf inf]);
	%set(gca,'YTickMode','manual');
	%set(gca,'XTick', 0:0.1:0.5);
	%set(gca,'YTick', 0:0.1:0.5);
	%legend('Line 1','Line 2','Location','SouthEast','Orientation','Vertical');
	%legend boxoff;
subplot(5,1,3);
	h=plot(t, Var3, 'k');
		set(h,'linewidth',2);
	xlabel({'Time [s]','(c)'},'FontSize',12);
	ylabel('Variable 3 [-]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	%axis([0 10 -inf inf]);
	%set(gca,'YTickMode','manual');
	%set(gca,'XTick', 0:0.1:0.5);
	%set(gca,'YTick', 0:0.1:0.5);
	%legend('Line 1','Line 2','Location','SouthEast','Orientation','Vertical');
	%legend boxoff;
subplot(5,1,4);
	h=plot(t, Var4, 'k');
		set(h,'linewidth',2);
	xlabel({'Time [s]','(d)'},'FontSize',12);
	ylabel('Variable 4 [-]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	%axis([0 10 -inf inf]);
	%set(gca,'YTickMode','manual');
	%set(gca,'XTick', 0:0.1:0.5);
	%set(gca,'YTick', 0:0.1:0.5);
	%legend('Line 1','Line 2','Location','SouthEast','Orientation','Vertical');
	%legend boxoff;
subplot(5,1,5);
	h=plot(t, Var5, 'k');
		set(h,'linewidth',2);
	xlabel({'Time [s]','(e)'},'FontSize',12);
	ylabel('Variable 5 [-]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	%axis([0 10 -inf inf]);
	%set(gca,'YTickMode','manual');
	%set(gca,'XTick', 0:0.1:0.5);
	%set(gca,'YTick', 0:0.1:0.5);
	%legend('Line 1','Line 2','Location','SouthEast','Orientation','Vertical');
	%legend boxoff;

% EPS�t�@�C������(���[�J���Ŏ��s�̂���)
% print(gcf,'-depsc2','-tiff',strcat(FileName,'.eps'));

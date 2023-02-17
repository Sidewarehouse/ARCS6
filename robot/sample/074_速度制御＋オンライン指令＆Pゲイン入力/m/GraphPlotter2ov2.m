% ARCS���o�͂���DATA.csv��MATLAB�ɓǂݍ��ރX�N���v�g�̈��
% 2020/06/11 Yokokura, Yuki
clc;
clear;

% CSV�t�@�C�����ݒ�
FileName = '../DATA.csv';

% CSV�t�@�C������ϐ��l�ǂݍ���
CsvData  = csvread(FileName);
t = CsvData(:,1);
A1 = CsvData(:,2);
A2 = CsvData(:,3);
B1 = CsvData(:,4);
B2 = CsvData(:,5);
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
subplot(2,1,1);
	h=plot(t, A1, 'k');
		set(h,'linewidth',4);
	hold on;
	h=plot(t, A2, 'r');
		set(h,'linewidth',2);
	hold off;
	xlabel({'Time [s]','(a)'},'FontSize',12);
	ylabel('Variable A [-]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	%axis([0 10 -inf inf]);
	%set(gca,'YTickMode','manual');
	%set(gca,'XTick', 0:0.1:0.5);
	%set(gca,'YTick', 0:0.1:0.5);
	legend('A1','A2','Location','SouthEast','Orientation','Vertical');
	%legend boxoff;
subplot(2,1,2);
	h=plot(t, B1, 'k');
		set(h,'linewidth',4);
	hold on;
	h=plot(t, B2, 'r');
		set(h,'linewidth',2);
	hold off;
	xlabel({'Time [s]','(b)'},'FontSize',12);
	ylabel('Variable B [-]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	%axis([0 10 -inf inf]);
	%set(gca,'YTickMode','manual');
	%set(gca,'XTick', 0:0.1:0.5);
	%set(gca,'YTick', 0:0.1:0.5);
	legend('B1','B2','Location','SouthEast','Orientation','Vertical');
	%legend boxoff;

% EPS�t�@�C������(���[�J���Ŏ��s�̂���)
% print(gcf,'-depsc2','-tiff',strcat(FileName,'.eps'));

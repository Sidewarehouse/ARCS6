% ARCS���o�͂���DATA.csv��MATLAB�ɓǂݍ��ރX�N���v�g�̈��
% ���A���^�C�����`�F�b�N�p�̃O���t�����X�N���v�g
% 2020/04/13 Yokokura, Yuki
clc;
clear;

% CSV�t�@�C�����ݒ�
FileName = '../DATA.csv';

% CSV�t�@�C������ϐ��l�ǂݍ���
CsvData  = csvread(FileName);
t = CsvData(:,1);
Tcmp = CsvData(:,2);	% [s] �����
Tact = CsvData(:,3);	% [s] ���ۂ̌v�������������
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
set(gcf,'renderer','Painters');
set(gcf,'PaperPositionMode','auto');
set(gcf,'color',[1 1 1]);
%set(gcf,'Position',[100 100 800 900]);
subplot(2,2,1);
	h=plot(t, Tcmp*1e6, 'k+');
		set(h,'linewidth',2);
	xlabel({'Time t [s]','(a)'},'FontSize',12);
	ylabel('Consumption Time Tcmp [\mus]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	%axis([0 10 -inf inf]);
	%set(gca,'YTickMode','manual');
	%set(gca,'XTick', 0:0.1:0.5);
	%set(gca,'YTick', 0:0.1:0.5);
	%legend('Consumption Tcmp','Actual Periodic Tact','Location','SouthEast','Orientation','Vertical');
	%legend boxoff;
subplot(2,2,3);
	h=plot(t, Tact*1e6, 'r+');
		set(h,'linewidth',2);
	xlabel({'Time t [s]','(b)'},'FontSize',12);
	ylabel('Actual Periodic Time Tact [\mus]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	%axis([0 inf 0 inf]);
	%set(gca,'YTickMode','manual');
	%set(gca,'XTick', 0:0.1:0.5);
	%set(gca,'YTick', 0:0.1:0.5);
	%legend('Line 1','Line 2','Location','SouthEast','Orientation','Vertical');
	%legend boxoff;
subplot(2,2,[2,4]);
	h = histogram(Tact*1e6);
		set(h, 'BinWidth', 0.5);
		set(h, 'FaceColor', 'red');
		set(h, 'EdgeColor', 'red');
	set(gca,'yscale','log');
	set(gca,'FontSize',12);
	xlabel({'Actual Periodic Time Tact [\mus]','(c)'},'FontSize',12);
	ylabel('Number of Control Loops [-]','FontSize',12);
	axis([90 inf 0.1 inf]);
	grid on;
	title('Histogram of Actual Periodic Time');
	
% EPS�t�@�C������
print(gcf,'-depsc2','-tiff','-r300','-painters',strcat(FileName,'.eps'));


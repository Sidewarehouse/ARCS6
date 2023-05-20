% ARCS���o�͂���DATA.csv��MATLAB�ɓǂݍ��ރX�N���v�g�̈��
% ��/�ʒu�n�C�u���b�h�����
% 2022/03/09 Yokokura, Yuki
clc;
clear;

% CSV�t�@�C�����ݒ�
FileName = '../DATA.csv';

% CSV�t�@�C������ϐ��l�ǂݍ���
CsvData  = csvread(FileName);
t = CsvData(:,1);
px = CsvData(:,2);
py = CsvData(:,3);
pz = CsvData(:,4);
fx = CsvData(:,5);
fy = CsvData(:,6);
fz = CsvData(:,7);
pxref = CsvData(:,8);
pyref = CsvData(:,9);
fzref = CsvData(:,10);
clear CsvData;
tlen = length(t);

% �f�[�^���Ԉ����ꍇ�̗�
%{
RedRate = 20;	% �Ԉ����v�f��
t = t(1:RedRate:tlen);
A = A(1:RedRate:tlen);
tlen = length(t);
%}

% �O���t�`��
figure(1);
	Width  = 1000;					% [px] ��
	Height = 1000;					% [px] ����
	FontSize = 14;					% [pt] �����T�C�Y
	FontName = 'Times New Roman';	% �t�H���g��
	AxisColor = [0,0,0];			% ���̐F
	BackColor = [1,1,1];			% �w�i�F
	clf;
	set(gcf,'PaperPositionMode','auto','Position',[100 10 Width Height],'Color',BackColor);
subplot(3,1,1);
	h = stairs(t, pxref*1e3, 'k');
		set(h,'linewidth',3);
	hold on;
	h = stairs(t, px*1e3, 'r');
		set(h,'linewidth',2);
	hold off;
	xlabel({'Time [s]','(a)'},'FontSize',FontSize,'FontName',FontName);
	ylabel({'Position','X-axis [mm]'},'FontSize',FontSize,'FontName',FontName);
	set(gca,'FontSize',FontSize,'FontName',FontName,'Color',BackColor,'XColor',AxisColor,'YColor',AxisColor);
	grid on;
	%axis([0 3 -1e-3 1e-3]);
	%set(gca,'XTickMode','manual','YTickMode','manual','XTick', -50:50:300,'YTick', -1:1:6);
	legend('Reference','Response','Location','NorthWest','Orientation','Vertical','TextColor',AxisColor);
subplot(3,1,2);
	h = stairs(t, pyref*1e3, 'k');
		set(h,'linewidth',3);
	hold on;
	h = stairs(t, py*1e3, 'r');
		set(h,'linewidth',2);
	hold off;
	xlabel({'Time [s]','(b)'},'FontSize',FontSize,'FontName',FontName);
	ylabel({'Position','Y-axis [mm]'},'FontSize',FontSize,'FontName',FontName);
	set(gca,'FontSize',FontSize,'FontName',FontName,'Color',BackColor,'XColor',AxisColor,'YColor',AxisColor);
	grid on;
	%axis([0 3 -0.25 0.05]);
	%set(gca,'XTickMode','manual','YTickMode','manual','XTick', -50:50:300,'YTick', -1:1:6);
	legend('Reference','Response','Location','SouthWest','Orientation','Vertical','TextColor',AxisColor);
subplot(3,1,3);
	h = stairs(t, fzref, 'k');
		set(h,'linewidth',3);
	hold on;
	h = stairs(t, fz, 'r');
		set(h,'linewidth',2);
	hold off;
	xlabel({'Time [s]','(c)'},'FontSize',FontSize,'FontName',FontName);
	ylabel({'Position','Z-axis [N]'},'FontSize',FontSize,'FontName',FontName);
	set(gca,'FontSize',FontSize,'FontName',FontName,'Color',BackColor,'XColor',AxisColor,'YColor',AxisColor);
	grid on;
	%axis([0 3 -0.2 1.4]);
	%set(gca,'XTickMode','manual','YTickMode','manual','XTick', -50:50:300,'YTick', -1:1:6);
	legend('Reference','Response','Location','NorthEast','Orientation','Vertical','TextColor',AxisColor);

% PNG,EPS�t�@�C�������̗�
%saveas(gcf, strcat(FileName,'.png'));
%print(gcf,'-depsc2','-tiff','-painters',strcat(FileName,'.eps'));


% ARCS���o�͂���DATA.csv��MATLAB�ɓǂݍ��ރX�N���v�g�̈��
% 2020/04/06 Yokokura, Yuki
clc;
clear;

% CSV�t�@�C�����ݒ�
FileName = '../DATA.csv';

% CSV�t�@�C������ϐ��l�ǂݍ���
CsvData  = csvread(FileName);
t = CsvData(:,1);
A = CsvData(:,2);
B = CsvData(:,3);
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
h=plot(t, A, 'k');
	set(h,'linewidth',4);
hold on;
h=plot(t, B, 'r');
	set(h,'linewidth',2);
hold off;
xlabel('Time [s]','FontSize',12);
ylabel(' Variable A [-]','FontSize',12);
set(gca,'FontSize',12);
grid on;
%axis([0 10 -inf inf]);
%set(gca,'YTickMode','manual');
%set(gca,'XTick', 0:0.1:0.5);
%set(gca,'YTick', 0:0.1:0.5);
legend('A','B','Location','SouthEast','Orientation','Vertical');
%legend boxoff;

% EPS�t�@�C������(���[�J���Ŏ��s�̂���)
% print(gcf,'-depsc2','-tiff',strcat(FileName,'.eps'));

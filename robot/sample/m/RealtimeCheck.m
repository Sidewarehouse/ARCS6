% ���A���^�C�����̑���f�[�^�v���b�g
% 2024/06/25 Yokokura, Yuki
clc;
clear;

% CSV�t�@�C������ǂݍ���
FileName = '../011_���A���^�C�����̑���/DATA.csv';	% CSV�t�@�C�����ݒ�
import MultiPlot.LoadCsvFile;			% LoadCsvFile�֐��C���|�[�g
[t, Tcmp, Tact] = LoadCsvFile(FileName);% �ϐ��l�ǂݍ���

% �O���t���Ԏ��I�t�Z�b�g���X�P�[�����O
%t = t - 0.05;	% �J�n�����̃I�t�Z�b�g
%t = t*1e3;		% [s] �� [ms] �̃X�P�[�����O

% �O���t�`��
figure(1);
	% MultiPlot�S�̂̐ݒ�
	Graph1 = MultiPlot(gcf);			% MultiPlot����
	Graph1.FigurePosition(55, 55);		% Figure�ʒu�̐ݒ�(��[px], ��[px])
	Graph1.FigureSize(900, 800);		% Figure�T�C�Y�̐ݒ�(��[px], ����[px])
	Graph1.FigureMargin(130, 70, 20);	% Figure�]���̐ݒ�(����[px], ����[px], �E��[px])
	Graph1.Font('Times New Roman', 18);	% �t�H���g�̐ݒ�(�t�H���g��, �t�H���g�T�C�Y)
	Graph1.NumOfPlanes(3);				% �v���b�g���ʂ̒i��
	%Graph1.XaxisRange(-1, 1, 5);		% X���͈͂̐ݒ�(�ŏ��l, �O���b�h�Ԋu, �ő�l) ���R�����g�A�E�g�Ŏ������[�h
	Graph1.XaxisLabel('Time [s]');		% X�����x����
	% �v���b�g����1�i��
	Graph1.SelectPlane(1);							% �v���b�g�i�I��
	Graph1.StairsPlot(t, Tcmp*1e6, 'Black', 'Thin');% �v�Z����Ԃ̃v���b�g
	Graph1.ManualGrid(0, 25, 100);					% Y���͈͂̐ݒ�(�ŏ��l, �O���b�h�Ԋu, �ő�l) ���R�����g�A�E�g�Ŏ������[�h
	Graph1.Label({'(a)','Consumption Time','Tcmp [\mus]'});			% �c�����x���̐ݒ�
	% �v���b�g����2�i��
	Graph1.SelectPlane(2);							% �v���b�g�i�I��
	Graph1.StairsPlot(t, Tact*1e6, 'Black', 'Thin');% ���ۂ̐���������Ԃ̃v���b�g
	Graph1.ManualGrid(0, 50, 200);					% Y���͈͂̐ݒ�(�ŏ��l, �O���b�h�Ԋu, �ő�l) ���R�����g�A�E�g�Ŏ������[�h
	Graph1.Label({'(b)','Actual Periodic Time','Tact [\mus]'});		% �c�����x���̐ݒ�
	% �v���b�g����3�i��
	Graph1.SelectPlane(3);							% �v���b�g�i�I��
	hold on;
		h = histogram(Tact*1e6);
			%set(h, 'BinWidth', 0.5);
			set(h, 'FaceColor', 'black');
			set(h, 'EdgeColor', 'black');
		set(gca,'yscale','log');
		set(gca,'FontSize',18);
		xlabel('Actual Periodic Time Tact [\mus]','FontSize',18);
		ylabel({'(c)','Number of','Control Loops [-]'},'FontSize',18);
		axis([90 inf 0.1 inf]);
		grid on;
	hold off;
	% �摜����
	%Graph1.SavePNGandEPS('DATA.eps');		% PNG�摜��EPS�t�@�C���𐶐�
	%Graph1.SavePNGandPDF('DATA.pdf');		% PNG�摜��PDF�t�@�C���𐶐�


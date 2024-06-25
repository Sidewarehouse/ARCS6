% MultiPlot�X�N���v�g�̃T���v���R�[�h
% 2024/06/25 Yokokura, Yuki
clc;
clear;

% CSV�t�@�C������ǂݍ���
FileName = '../BaseCtrl/DATA.csv';	% CSV�t�@�C�����ݒ�
import MultiPlot.LoadCsvFile;		% LoadCsvFile�֐��C���|�[�g
[t, x1, x2, x3, x4, x5, x6, x7, x8, x9] = LoadCsvFile(FileName);	% �ϐ��l�ǂݍ���

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
	Graph1.NumOfPlanes(2);				% �v���b�g���ʂ̒i��
	%Graph1.XaxisRange(-1, 1, 5);		% X���͈͂̐ݒ�(�ŏ��l, �O���b�h�Ԋu, �ő�l) ���R�����g�A�E�g�Ŏ������[�h
	Graph1.XaxisLabel('Time [s]');		% X�����x����
	% �v���b�g����1�i��
	Graph1.SelectPlane(1);							% �v���b�g�i�I��
	Graph1.StairsPlot(t, x1, 'Black', 'Thin');		% �K�i�v���b�g
	Graph1.StairsPlot(t, x2, 'Red',   'Normal');	% ���̐F�� Black/Gray/Red/Green/Blue ��5��̂�
	Graph1.StairsPlot(t, x3, 'Green', 'Normal');	% ���̑����� Thin/Normal/Heavy ��3��̂�
	Graph1.StairsPlot(t, x4, 'Blue',  'Heavy');		% �d�˂���v���b�g�̐��͖�����
	%Graph1.ManualGrid(-2, 0.5, 2);					% Y���͈͂̐ݒ�(�ŏ��l, �O���b�h�Ԋu, �ő�l) ���R�����g�A�E�g�Ŏ������[�h
	Graph1.Label({'(a)','Var. A','[-]'});			% �c�����x���̐ݒ�
	Graph1.Legend({'x_1','x_2','x_3','x_4'}, 'NorthEast', 'Vertical');	% �}��̐ݒ�
	% �v���b�g����2�i��
	Graph1.SelectPlane(2);							% �v���b�g�i�I��
	Graph1.LinePlot(t, x5, 'Black', 'Thin');		% ���`�v���b�g
	Graph1.LinePlot(t, x6, 'Red',   'Normal');		% ���̐F�� Black/Gray/Red/Green/Blue ��5��̂�
	Graph1.LinePlot(t, x7, 'Green', 'Normal');		% ���̑����� Thin/Normal/Heavy ��3��̂�
	Graph1.LinePlot(t, x8, 'Blue',  'Heavy');		% �d�˂���v���b�g�̐��͖�����
	Graph1.LinePlot(t, x9, 'Blue',  'Heavy');		% �d�˂���v���b�g�̐��͖�����
	%Graph1.ManualGrid(-2, 0.5, 2);					% Y���͈͂̐ݒ�(�ŏ��l, �O���b�h�Ԋu, �ő�l) ���R�����g�A�E�g�Ŏ������[�h
	Graph1.Label({'(b)','Var. B','[-]'});			% �c�����x���̐ݒ�
	Graph1.Legend({'x_5','x_6','x_7','x_8','x_9'}, 'NorthEast', 'Vertical');	% �}��̐ݒ�
	% �摜����
	Graph1.SavePNGandEPS('DATA.eps');		% PNG�摜��EPS�t�@�C���𐶐�
	Graph1.SavePNGandPDF('DATA.pdf');		% PNG�摜��PDF�t�@�C���𐶐�


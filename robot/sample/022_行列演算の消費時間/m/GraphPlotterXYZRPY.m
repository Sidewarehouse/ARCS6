% ARCS���o�͂���DATA.csv��MATLAB�ɓǂݍ��ރX�N���v�g�̈��
% 6�����{�b�g�p
% 2020/03/01 Yuki YOKOKURA
clc;
clear;

% CSV�t�@�C�����ݒ�
FileName = '../DATA.csv';

% CSV�t�@�C������ϐ��l�ǂݍ���
CsvData  = csvread(FileName);
t = CsvData(:,1);   % [s] ����
Px = CsvData(:,2);  % [m] X�ʒu
Py = CsvData(:,3);  % [m] Y�ʒu
Pz = CsvData(:,4);  % [m] Z�ʒu
Ar = CsvData(:,5);  % [rad] ���[���p
Ap = CsvData(:,6);  % [rad] �s�b�`�p
Ay = CsvData(:,7);  % [rad] ���[�p
clear CsvData;
tlen = length(t);

% �p���x�N�g���p�f�[�^���Ԉ���
RedRate = 200;	% �Ԉ����v�f��
Pxa = Px(1:RedRate:tlen);
Pya = Py(1:RedRate:tlen);
Pza = Pz(1:RedRate:tlen);
Ara = Ar(1:RedRate:tlen);
Apa = Ap(1:RedRate:tlen);
Aya = Ay(1:RedRate:tlen);

% �I�C���[�p����x�N�g���ɕϊ�
L = 0.3;    % �x�N�g���̒���
Vx = L*cos(Aya).*cos(-Apa);
Vy = L*sin(Aya).*cos(-Apa);
Vz = L*sin(-Apa);

% �O���t�`��
figure(1);
clf;
set(gcf,'PaperPositionMode','manual');
set(gcf,'color',[1 1 1]);
%set(gcf,'Position',[100 100 800 900]);
subplot(2,2,1);
	h=plot(t, Px, 'r', t, Py, 'k', t, Pz, 'b');
		set(h,'linewidth',2);
	xlabel({'Time [s]','(a)'},'FontSize',12);
	ylabel('Position [m]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	%axis([0 10 -inf inf]);
	%set(gca,'YTickMode','manual');
	%set(gca,'XTick', 0:0.1:0.5);
	%set(gca,'YTick', 0:0.1:0.5);
    legend('X','Y','Z','Location','SouthEast','Orientation','Vertical');
	%legend boxoff;
subplot(2,2,3);
	h=plot(t, Ar*180/pi, 'r', t, Ap*180/pi, 'k', t, Ay*180/pi, 'b');
		set(h,'linewidth',2);
	xlabel({'Time [s]','(b)'},'FontSize',12);
	ylabel('Attitude [deg]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	%axis([0 10 -inf inf]);
	%set(gca,'YTickMode','manual');
	%set(gca,'XTick', 0:0.1:0.5);
	%set(gca,'YTick', 0:0.1:0.5);
	legend('Roll','Pitch','Yaw','Location','SouthEast','Orientation','Vertical');
	%legend boxoff;
subplot(2,2,[2,4]);
    h = plot3(Px, Py, Pz, 'k');
		set(h,'linewidth',2);
    hold on;
    h = quiver3(Pxa, Pya, Pza, Vx, Vy, Vz, 0, 'r');
		set(h,'linewidth',2);
    h = plot3(0, 0, 0, 'ko');
		set(h,'linewidth',2);
    hold off;
    xlabel('Position X [m]','FontSize',12);
    ylabel('Position Y [m]','FontSize',12);
    zlabel('Position Z [m]','FontSize',12);
	set(gca,'FontSize',12);
	axis([-0.2 1.2 -1.2 1.2 0 2]);
	set(gca,'XTickMode','manual');
	set(gca,'YTickMode','manual');
	set(gca,'ZTickMode','manual');
	set(gca,'XTick', -0.2:0.2:1.2);
	set(gca,'YTick', -1.2:0.2:1.2);
	set(gca,'ZTick',  0.0:0.2:2.0);
    grid on;

% EPS�t�@�C������(���[�J���Ŏ��s�̂���)
% print(gcf,'-depsc2','-tiff',strcat(FileName,'.eps'));

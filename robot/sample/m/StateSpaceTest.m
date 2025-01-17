% ArcsControl::StateSpace�N���X�̃e�X�g�p�X�N���v�g
% 2025/01/17 Yokokura, Yuki
clc;
clear;

% CSV�t�@�C�����ݒ�
FileName = '../042_�C�ӂ̏�ԋ�ԃ��f��/DATA.csv';

% �e�X�g�Ώۂ̏�ԋ�ԃ��f��
A = [
	 0,  1 ;
	-5, -2
];
b = [
	0 ;
	3
];
C = [
	1, 0
];
D = [
	1
];
System1 = ss(A,b,C,D)

% CSV�t�@�C������ϐ��l�ǂݍ���
CsvData  = csvread(FileName);
t = CsvData(:,1);
u = CsvData(:,2);
y = CsvData(:,3);
clear CsvData;
tlen = length(t);

% �^�l�̌v�Z
Ts = t(tlen)/(tlen - 1);
ti = 0:Ts:t(tlen);
ui = SquareWave(0.5, pi, ti, 1);
yi = lsim(System1, ui, ti);

% �O���t�`��
figure(1);
clf;
set(gcf,'PaperPositionMode','manual');
set(gcf,'color',[1 1 1]);
subplot(2,1,1);
	h=stairs(t, u, 'k:');
		set(h,'linewidth',4);
	hold on;
	h=stairs(t, y, 'r');
		set(h,'linewidth',3);
	h=stairs(ti, ui, 'g:');
		set(h,'linewidth',2);
	h=stairs(ti, yi, 'k');
		set(h,'linewidth',1);
	hold off;
	xlabel('Time [s]','FontSize',12);
	ylabel('Input and Output [-]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	%axis([0 10 -inf inf]);
	legend('ARCS Input','ARCS Output','MATLAB Input','MATLAB Output','Location','SouthEast','Orientation','Vertical');
subplot(2,1,2);
	h=stairs(t, ti, 'k');
		set(h,'linewidth',1);
	xlabel('Ideal Time [s]','FontSize',12);
	ylabel('Actual Time [s]','FontSize',12);
	grid on;


% EPS�t�@�C������(���[�J���Ŏ��s�̂���)
% print(gcf,'-depsc2','-tiff',strcat(FileName,'.eps'));

% ���`�g�����֐�
function y = SquareWave(freq, phase, time, starttime)
	tlen = length(time);
	r(1:tlen) = 0;
	y(1:tlen) = 0;
	
	% �J�n���Ԃ����V�t�g���ꂽ�����g�𐶐�
	for i = 1:tlen
		if time(i) < starttime
			r(i) = 0;			% �J�n�������O�̂Ƃ��̓[���o��
		else
			r(i) = sin(2*pi*freq*(time(i) - starttime) + phase);	% �J�n�����ȍ~�͕��`�g�o��
			% ���`�g�ɕϊ�
			if 0 <= r(i)
				y(i) =  1;
			else
				y(i) = -1;
			end
		end
	end
end



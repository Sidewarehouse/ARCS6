% TwoInertiaSimulator�N���X�̃e�X�g�p�X�N���v�g
% 2023/10/27 Yokokura, Yuki
clc;
clear;

% CSV�t�@�C�����ݒ�
FileName = '../DATA.csv';

% 2�����n�p�����[�^
Jl = 0.1;	% [kgm^2]		���ב�����
Dl = 1e-12;	% [Nm/(rad/s)]	���ב��S��
Ds = 1;		% [Nm/(rad/s)]	�˂���S��
Ks = 10000;	% [Nm/rad]		�˂��ꍄ��
Jm = 1e-5;	% [kgm^2]		���[�^������
Dm = 1e-4;	% [Nm/(rad/s)]	���[�^���S��
Rg = 100;	% [-]			������
Kt = 0.1;	% [Nm/A]		�g���N�萔

% CSV�t�@�C������ϐ��l�ǂݍ���
CsvData  = csvread(FileName);
t    = CsvData(:,1);
iq   = CsvData(:,2);
tdis = CsvData(:,3);
wl   = CsvData(:,4);
ths  = CsvData(:,5);
wm   = CsvData(:,6);
taus = CsvData(:,7);
clear CsvData;
tlen = length(t);

% �e�X�g�Ώۂ̏�ԋ�ԃ��f��
% 2�������U�n�̏�ԋ�ԃ��f���ł̕\��
% xp = [wl, ths, wm]^T
Ap = [
	-(Dl + Ds)/Jl,	Ks/Jl,			Ds/(Jl*Rg)			;
	-1,				0,				1.0/Rg				;
	Ds/(Jm*Rg),		-Ks/(Jm*Rg),	-(Ds/Rg^2 + Dm)/Jm
];
bp = [
	0		, -1/Jl ;
	0		,  0    ;
	Kt/Jm	,  0
]
Cp = eye(3);
dp = zeros(3,2);
System1 = ss(Ap,bp,Cp,dp)

% �^�l�����̌v�Z
Ts = t(tlen)/(tlen - 1);
ti = 0:Ts:t(tlen);
iqi   = 0.5*SquareWave(0.1, 0, ti, 0.5);
tdisi =   2*SquareWave(0.2, 0, ti, 2.0);
ui = [iqi ; tdisi];
yi = lsim(System1, ui, ti);
wli  = yi(:,1);
thsi = yi(:,2);
wmi  = yi(:,3);
tausi = -Ds*wli + Ks*thsi + Ds/Rg*wmi;

% �O���t�`��
figure(1);
clf;
set(gcf,'PaperPositionMode','manual');
set(gcf,'color',[1 1 1]);
subplot(6,1,1);
	h=stairs(ti, iqi, 'k');
		set(h,'linewidth',3);
	hold on;
	h=stairs(t, iq, 'r');
		set(h,'linewidth',1);
	hold off;
	xlabel('Time [s]','FontSize',12);
	ylabel('q-axis Current [A]','FontSize',12);
	set(gca,'FontSize',12);
	legend('MATLAB','ARCS','Location','SouthEast','Orientation','Vertical');
	grid on;
subplot(6,1,2);
	h=stairs(ti, tdisi, 'k');
		set(h,'linewidth',3);
	hold on;
	h=stairs(t, tdis, 'r');
		set(h,'linewidth',1);
	hold off;
	xlabel('Time [s]','FontSize',12);
	ylabel({'Disturbance','Torque [Nm]'},'FontSize',12);
	set(gca,'FontSize',12);
	grid on;
subplot(6,1,3);
	h=stairs(ti, wli, 'k');
		set(h,'linewidth',3);
	hold on;
	h=stairs(t, wl, 'r');
		set(h,'linewidth',1);
	hold off;
	xlabel('Time [s]','FontSize',12);
	ylabel({'Load-side','Velocity [rad/s]'},'FontSize',12);
	set(gca,'FontSize',12);
	grid on;
subplot(6,1,4);
	h=stairs(ti, thsi, 'k');
		set(h,'linewidth',3);
	hold on;
	h=stairs(t, ths, 'r');
		set(h,'linewidth',1);
	hold off;
	xlabel('Time [s]','FontSize',12);
	ylabel({'Torsion Angle [rad]'},'FontSize',12);
	set(gca,'FontSize',12);
	grid on;
subplot(6,1,5);
	h=stairs(ti, wmi, 'k');
		set(h,'linewidth',3);
	hold on;
	h=stairs(t, wm, 'r');
		set(h,'linewidth',1);
	hold off;
	xlabel('Time [s]','FontSize',12);
	ylabel({'Motor-side','Velocity [rad/s]'},'FontSize',12);
	set(gca,'FontSize',12);
	grid on;
subplot(6,1,6);
	h=stairs(ti, tausi, 'k');
		set(h,'linewidth',3);
	hold on;
	h=stairs(t, taus, 'r');
		set(h,'linewidth',1);
	hold off;
	xlabel('Time [s]','FontSize',12);
	ylabel({'Torsion','Torque [Nm]'},'FontSize',12);
	set(gca,'FontSize',12);
	grid on;

figure(2);
set(gcf,'PaperPositionMode','manual');
set(gcf,'color',[1 1 1]);
	h=stairs(t, ti, 'k');
		set(h,'linewidth',1);
	xlabel('Ideal Time [s]','FontSize',12);
	ylabel('Actual Time [s]','FontSize',12);
	set(gca,'FontSize',12);
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



%
% ��������̍ő卽��������fa����A���S���Y��
%						2019/09/18 Yuki YOKOKURA
clc;
clear;

% CSV�t�@�C�����ݒ�
FileName = '../DATA.csv';

% �p�����[�^�ݒ�
Ts   = 250e-6;	% [s] �T���v�����O����
P    = 4;		% [-] �ɑΐ�
Ra   = 0.671;	% [��] �d�@�q��R
La   = 1.40e-3;	% [H]  �d�@�q�C���_�N�^���X

% CSV�t�@�C������ϐ��l�ǂݍ���
CsvData  = csvread(FileName);
t        = CsvData(:,1);
theta_rm = CsvData(:,2);
Ia       = CsvData(:,3);
Ib       = CsvData(:,4);
Va       = CsvData(:,5);
Vb       = CsvData(:,6);
clear CsvData;
tlen = length(t);

% �d�C�p���x�̌v�Z
g = 6280;	% [rad/s] ����ш�
s = tf('s');
Gpd = s*g/(s + g);
tsim = (0:Ts:(tlen - 1)*Ts).';
omega_re = lsim(Gpd, theta_rm - theta_rm(1), tsim)*P;


% ���f�[�^�̕\��
figure(1);
clf;
subplot(3,1,1);
	h = plot(t, theta_rm);
	set(h,'LineWidth',2);
	xlabel('Time [s]');
	ylabel('Position theta\_rm [rad]');
	grid on;
subplot(3,1,2);
	h = plot(t, omega_re);
	set(h,'LineWidth',2);
	xlabel('Time [s]');
	ylabel('Electric Speed omega\_re [rad/s]');
	grid on;
subplot(3,1,3);
	h = plot(t, Ia, t, Ib);
	set(h,'LineWidth',2);
	xlabel('Time [s]');
	ylabel('Current Ia, Ib [A]');
	grid on;

% �U�N�d������ƐU���v�Z
L = g/(s + g);
Pn = 1/(La*s + Ra);
Left  = L;
Right = minreal(L/Pn);
Loa = lsim(Left, Va, tsim);
Roa = lsim(Right, Ia, tsim);
Ea = Loa - Roa;
Lob = lsim(Left, Vb, tsim);
Rob = lsim(Right, Ib, tsim);
Eb = Lob - Rob;
Emax = sqrt(Ea.^2 + Eb.^2);

% �U�N�d���g�`�ƐU���g�`�̕\��
figure(2);
subplot(2,1,1);
	h = plot(tsim, Ea, tsim, Eb);
	set(h,'LineWidth',2);
	xlabel('Time [s]');
	ylabel('Back EMF Ea, Eb [V]');
	grid on;
subplot(2,1,2);
	h = plot(tsim, Emax);
	set(h,'LineWidth',2);
	xlabel('Time [s]');
	ylabel('Back EMF Amplitude [V]');
	grid on;

% �ŏ����@�ɂ��ő卽�������̓���
x = abs(omega_re);
y = Emax;
A(1:length(x),1) = x;
%A(1:length(x),2) = 1; % �d���I�t�Z�b�g�����肷��Ƃ��̂݃R�����g�A�E�g
u = A\y;
Phifa = u(1)	% [V/(rad/s)] �ő卽������
%b = u(2)		% [V]		  �d���I�t�Z�b�g
Wrefit  = linspace(0, max(abs(omega_re)));
Emaxfit = Phifa*Wrefit;
figure(3);
	h = plot(abs(omega_re), Emax, 'k+', Wrefit, Emaxfit, 'r');
	set(h,'LineWidth',2);
	xlabel('Elec. Anglular Velocity [rad/s]');
	ylabel('Back EMF Amplitude [V]');
    title(strcat('Max. Flux Linkage Phifa = ', sprintf(' %e', Phifa),' [V/(rad/s)]'));
	grid on;


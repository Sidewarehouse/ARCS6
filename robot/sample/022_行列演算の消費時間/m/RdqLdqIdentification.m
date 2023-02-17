%
% dq���̒�R�ƃC���_�N�^���X�̓���A���S���Y��
%						2019/09/18 Yuki YOKOKURA
clc;
clear;

% CSV�t�@�C�����ݒ�
FileName = '../DATA.csv';

% �p�����[�^�ݒ�
Ts       = 250e-6;	% [s] �T���v�����O����
Tfcs_sta = 5.5;		% [s] �v���Ɏg���f�[�^�̊J�n����
Tfcs_end = 6.9;		% [s] �v���Ɏg���f�[�^�̏I������
Tstep    = 6.0;		% [s] �X�e�b�v���͊J�n����
Tstd_sta = 6.5;	% [s] ����ԊJ�n����
Tstd_end = 6.9;	% [s] ����ԏI������

% CSV�t�@�C������ϐ��l�ǂݍ���
CsvData  = csvread(FileName);
t   = CsvData(:,1);
Vdq = CsvData(:,2);
Idq = CsvData(:,3);
clear CsvData;
tlen = length(t);

% ���f�[�^�̕\��
figure(1);
clf;
subplot(2,1,1);
	h = plot(t, Vdq);
	set(h,'LineWidth',2);
	xlabel('Time [s]');
	ylabel('Voltage Vd or Vq [V]');
	grid on;
subplot(2,1,2);
	h = plot(t, Idq);
	set(h,'LineWidth',2);
	xlabel('Time [s]');
	ylabel('Current Id or Iq [A]');
	grid on;

% �f�[�^�؂�o���C���d���l�̌v����63.2���̌v�Z
Iastd = mean(Idq((Tstd_sta/Ts):(Tstd_end/Ts)));
Ia632 = Iastd*(1 - exp(-1));
t = t((Tfcs_sta/Ts):(Tfcs_end/Ts)) - Tstep;
Ia = Idq((Tfcs_sta/Ts):(Tfcs_end/Ts));

% 63.2���̏ꏊ�̒T��
ep = 0.01;	% �T���̋��e��
idx = round(mean(find( (Ia632 - ep) < Ia & Ia < (Ia632 + ep) ))); % ������63.2�����炢�ɓ���v�f�ԍ����擾
if(isnan(idx) == true)
	warning('eq������������63.2���t�߂̃f�[�^��������܂���Bep��傫�����ĉ������B');
else
	tau = t(idx);   % ���萔���擾
end;

% ��R�ƃC���_�N�^���X�̌v�Z
Rdq = max(Vdq)/Iastd;
Ldq = Rdq*tau;

% �v���f�[�^�\��
figure(2);
	h = plot(t, Ia, [t(1) t(length(t))], [Iastd Iastd], [t(1) t(length(t))], [Ia632 Ia632], [0 0], [0 Iastd], [tau tau], [0 Iastd]);
	set(h,'LineWidth',2);
	xlabel('Time [s]');
	ylabel('Current Id or Iq [A]');
    title(strcat('Time constant \tau = ',sprintf(' %e',tau),' [s]  Resistance Rd or Rq = ',sprintf(' %e', Rdq),' [Ohm]  Inductance Ld or Lq = ',sprintf(' %e', Ldq),' [H]'));
	grid on;

% EPS�t�@�C������
print(gcf,'-depsc2','-tiff','RdqLdqIdentResults.eps');

% ����l�ł̃V�~�����[�V����
s = tf('s');
P = 1/(Ldq*s + Rdq);
Iasim = lsim(P, Vdq((Tfcs_sta/Ts):(Tfcs_end/Ts)), [0:Ts:((Tfcs_end/Ts) - (Tfcs_sta/Ts))*Ts]);
figure(3);
	h = plot(t, Ia,'k', t, Iasim,'r');
	set(h,'LineWidth',2);
	xlabel('Time [s]');
	ylabel('Current Id or Iq [A]');
	legend('Measured','Identified');
	grid on;


% ARCS-FRA�p�̎��g����������X�N���v�g V2021A
% Yokokura, Yuki & Muto Hirotaka & Shunsuke Suzuki 2021/12/22
%
% �Q�l�����F ���� �O, �t�R �u�Y:�u�𗬃C���s�[�_���X����@�̑I�� -FRA��FFT-�v�h�H�Z�p, vol. 35, no. 5, Jan. 1986.
%
% �g�����F
% 1. FRAgenerator�N���X��p���Ď��@���N�C����f�[�^��DATA.csv�ɕۑ�
%              ��
% 2. ���̃X�N���v�g����FRA�̑���p�����[�^��FRAgenerator�̐ݒ�ƍ��킹��
%              ��
% 3. ���̃X�N���v�g�����s�����FRA�̌v�Z���s�����g��������������
clc;
clear all;

%% FRA�̑���p�����[�^�̐ݒ�(���@�̃p�����[�^�ɍ��킹��)
FileName = '../DATA.csv';	% CSV�t�@�C����
Ts = 100e-6;% [s]  �T���v�����O����
Fsta =   3;	% [Hz] �J�n���g��
Fend = 100;	% [Hz] �I�����g��
Fstep = 1;	% [Hz] ���g���X�e�b�v
Ni = 10;	% [-]  �ϕ�����
Au = 0.1;	% [-]  �U��
Bu = -0.45;	% [-]  �o�C�A�X
Tsta = 0;	% [s] FRA�J�n����

%% ���g���x�N�g������
freq = Fsta:Fstep:Fend;
flen = length(freq);

%% ���@�f�[�^�ǂݍ���
data = csvread(FileName);
t = data(:,1);     % [s]  ����
f = data(:,2);     % [Hz] ���g��
iqref = data(:,3); % [A]  �d���w�ߒl
wmres = data(:,4); % [rad]���x����
clear data;

%% �o�C�A�X�␳
iqref = iqref - Bu;

%% �J�n�������ȑO�̃f�[�^���폜����O����
t = t - Tsta;			% �J�n�������������ɃV�t�g
t(1:Tsta/Ts) = [];		% �J�n�������ȑO�͍폜
f(1:Tsta/Ts) = [];		% �J�n�������ȑO�͍폜
iqref(1:Tsta/Ts) = [];	% �J�n�������ȑO�͍폜
wmres(1:Tsta/Ts) = [];	% �J�n�������ȑO�͍폜
tlen = length(t);

%% ���f�[�^�̕\��
figure(1); set(gcf,'color',[1 1 1]);
subplot(3,1,1); plot(t,f);      %���g���̕ω�
	xlabel('Time [s]','FontSize',14,'FontName','Times New Roman')
	ylabel('Frequency [rad/s]','FontSize',14,'FontName','Times New Roman')
	set(gca,'FontSize',14);
	grid on;
subplot(3,1,2); plot(t,iqref);  %���͐M��
	xlabel('Time [s]','FontSize',14,'FontName','Times New Roman')
	ylabel('Current [A]','FontSize',14,'FontName','Times New Roman')
	set(gca,'FontSize',14);
	grid on;
subplot(3,1,3); plot(t,wmres);  %�o�͐M��
	xlabel('Time [s]','FontSize',14,'FontName','Times New Roman')
	ylabel('Velocity [rad/s]','FontSize',14,'FontName','Times New Roman')
	set(gca,'FontSize',14);
	grid on;

%% FRA�̌v�Z
y = wmres;
tini = 0;
Ar(1:flen) = 0;
Ai(1:flen) = 0;
j = 1;
for i=1:tlen-1,
	Ar(j) = Ar(j) + y(i)*cos(2*pi*f(i)*(t(i) - tini))*Ts;
	Ai(j) = Ai(j) + y(i)*sin(2*pi*f(i)*(t(i) - tini))*Ts;
	if(f(i) ~= f(i+1))
		if(j < flen)
			Ar(j) = 2*f(i)/Ni*Ar(j);
			Ai(j) = 2*f(i)/Ni*Ai(j);
			tini = t(i+1);
			j = j + 1;
		else
			break;
		end;
	end
end
%figure(127); plot(freq,Ar,'x-', freq,Ai,'x-');	% �m�F�p

%% �Ō�̃f�[�^�������������Ȃ�̂ŏ����Ă���
Ar(flen) = [];
Ai(flen) = [];
freq(flen) = [];
flen = length(freq);

%% ���g�������̌v�Z
Ay = sqrt(Ar.^2 + Ai.^2);          % �o�͐U���v�Z
G = 20*log10(Ay./Au);              % �Q�C�������v�Z
P = unwrap(-atan2(Ai,Ar))*180/pi;  % �ʑ������v�Z

%% ���g�������̕`��
figure(2); set(gcf,'color',[1 1 1]);
subplot(2,1,1); h=semilogx(freq,G,'x-');
	set(h,'linewidth',2);
	xlabel('Frequency [Hz]','FontName','Times New Roman','FontSize',14)
	ylabel('Gain [dB]','FontName','Times New Roman','FontSize',14)
	set(gca,'FontSize',14);
	xlim([Fsta,Fend])
	grid on;
subplot(2,1,2); h=semilogx(freq,P,'x-');
	set(h,'linewidth',2);
	xlabel('Frequency [Hz]','FontName','Times New Roman','FontSize',14)
	ylabel('Phase [deg]','FontName','Times New Roman','FontSize',14)
	set(gca,'FontSize',14);
	xlim([Fsta,Fend])
	grid on;

% EPS�t�@�C������(���[�J���Ŏ��s�̂���)
% print(gcf,'-depsc2','-tiff','FRA.eps');

% �f�[�^�ۑ�
save(strcat(FileName, '_FRA.mat'), 'Ts', 'Fsta', 'Fend', 'Fstep', 'Ni', 'Au', 'Bu', 'Tsta', 'freq', 'G', 'P');


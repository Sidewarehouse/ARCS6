% ARCS-Control�̃e�X�g�p�X�N���v�g
% Yokokura, Yuki 2024/07/23
clc;
clear;
format short;
%format longE;

% �v�����g��ԋ�ԃ��f��(�K��)
Ap = [-1, -2, -9 ;
      0, -1, -7 ;
     -2,  6, -1 ]
bp = [1 ; 0 ; 0]
cp = [1, 0, 0]
dp = [0];
sys = ss(Ap, bp, cp ,dp);

fprintf('\n');
disp '�� �A�����A�v�m�t�������̉��̃e�X�g'
Qp = eye(3);
Xp = lyap(Ap, Qp)
O = norm(Ap*Xp + Xp*Ap' + Qp)	% �[���ɂȂ邩�m�F

fprintf('\n');
disp '�� ���U���A�v�m�t�������̉��̃e�X�g'
Xpd = dlyap(Ap, Qp)

fprintf('\n');
disp '�� ����/�ϑ��O���~�A���̃e�X�g'
Wc = gram(sys,'c')
Wo = gram(sys,'o')

fprintf('\n');
disp '�� ���t�����̃e�X�g'
[sys2,g,T,Ti] = balreal(sys);
Ah = sys2.a
bh = sys2.b
ch = sys2.c
Wc2 = gram(sys2,'c')
Wo2 = gram(sys2,'o')

fprintf('\n');
disp '�� ���U���̃e�X�g'
Ts = 100e-6;% [s]      �T���v�����O����
Jl = 1;		% [kgm^2]  ���ב�����
Dl = 0.3;	% [Nms/rad]���ב��S��
Ds = 1;		% [Nms/rad]�˂���S��
Ks = 500;	% [Nm/rad] 2�����Ԃ̂΂˒萔
Jm = 1e-4;	% [kgm^2]  ���[�^����
Dm = 0.1;	% [Nms/rad]���[�^�S��
Rg = 100;	%          ������
Kt = 0.04;	% [Nm/A]   �g���N�萔
Atc = [
		-(Dl + Ds)/Jl,	Ks/Jl,			Ds/(Jl*Rg)			  ;
		-1,				0,				1.0/Rg				  ;
		Ds/(Jm*Rg),		-Ks/(Jm*Rg),	-(Ds/(Rg*Rg) + Dm)/Jm ];
Btc = [
		0		, -1.0/Jl;
		0		, 0	     ;
		Kt/Jm	, 0		 ];
Ct = eye(3);
Ct2 = [ 0, 0, 1];
Dt = zeros(3,2);
sys3 = ss(Atc, Btc, Ct, Dt)
sys3d = c2d(sys3, Ts);	% ���̂܂ܗ��U������ꍇ
Ad = sys3d.a
bd = sys3d.b
Ob = obsv(Atc, Ct2)

%{

sys3 = balreal(sys3);	% ���t�����Ă��痣�U������ꍇ
sys3d = c2d(sys3, Ts);
Ad = sys3d.a
bd = sys3d.b
%}

fprintf('\n');
disp '�� ���U�n��ԋ�ԃ��f��'
load('../040_���䗝�_/ArcsControlTest.mat');
A = [-0.2516 -0.1684;2.784 0.3549];
B = [0;3];
C = [0 1];
D = 0;
Ts = 1;
sys1 = ss(A,B,C,D,Ts);
figure(1);
	step(sys1);
	hold on;
	plot(k1, y1, 'ro');
	plot(k1, y1n, 'go');
	hold off;
	ylabel('Output y1 [*]');
	xlabel('Time Step k1 [-]');
	grid on;
	legend('MATLAB','DiscStateSpace::GetOutput','DiscStateSpace::GetNextOutput');
	title('MATLAB �� ArcsControl �̔�r (DiscStateSpace�N���X)');
	
fprintf('\n');
disp '�� �A���n��ԋ�ԃ��f��'
figure(2);
	t_ = linspace(0,1,1000);
	y_ = step(sys3, t_);
subplot(3,1,1)
	plot(t_, y_(:,1), 'b');
	hold on;
	stairs(t2, y2(1,:), 'r');
	stairs(t2, y2n(1,:), 'g');
	hold off;
	ylabel('Output y2(1) [*]');
	xlabel('Time [s]');
	grid on;
	title('MATLAB �� ArcsControl �̔�r (StateSpace�N���X)');
subplot(3,1,2)
	plot(t_, y_(:,2), 'b');
	hold on;
	stairs(t2, y2(2,:), 'r');
	stairs(t2, y2n(2,:), 'g');
	hold off;
	ylabel('Output y2(2) [*]');
	xlabel('Time [s]');
	grid on;
subplot(3,1,3)
	plot(t_, y_(:,3), 'b');
	hold on;
	stairs(t2, y2(3,:), 'r');
	stairs(t2, y2n(3,:), 'g');
	hold off;
	ylabel('Output y2(3) [*]');
	xlabel('Time [s]');
	grid on;
	legend('MATLAB','StateSpace::GetOutput','StateSpace::GetNextOutput','Location','SouthEast');

fprintf('\n');
disp '�� �A���n�`�B�֐�'
figure(3);
	sys3 = tf([9],[1, 1.5, 9]);
	t_ = linspace(0,10,1000);
	y_ = step(sys3, t_);
	plot(t_, y_);
	hold on;
	stairs(t3, y3, 'r');
	stairs(t3, y3n, 'g');
	hold off;
	ylabel('Output y3 [*]');
	xlabel('Time [s]');
	grid on;
	legend('MATLAB','TransFunc::GetOutput','TransFunc::GetNextOutput','Location','SouthEast');
	title('MATLAB �� ArcsControl �̔�r (TransFunc�N���X)');


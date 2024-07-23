% ARCS-Control�̃e�X�g�p�X�N���v�g
% Yokokura, Yuki 2024/07/23
clc;
clear;
format short;
%format long;

% �v�����g��ԋ�ԃ��f��(�K��)
Ap = [-1, -2, -9 ;
      0, -1, -7 ;
     -2,  6, -1 ]
bp = [1 ; 0 ; 0]
cp = [1, 0, 0]
dp = [0];
sys = ss(Ap, bp, cp ,dp);

disp '�� �A�����A�v�m�t�������̉��̃e�X�g'
Qp = eye(3);
Xp = lyap(Ap, Qp)
O = norm(Ap*Xp + Xp*Ap' + Qp)	% �[���ɂȂ邩�m�F

disp '�� ���U���A�v�m�t�������̉��̃e�X�g'
Xpd = dlyap(Ap, Qp)


disp '�� ����/�ϑ��O���~�A���̃e�X�g'
Wc = gram(sys,'c')
Wo = gram(sys,'o')

disp '�� ���t�����̃e�X�g'
[sys2,g,T,Ti] = balreal(sys);
Ah = sys2.a
bh = sys2.b
ch = sys2.c
Wc2 = gram(sys2,'c')
Wo2 = gram(sys2,'o')

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
Ac = [
		-(Dl + Ds)/Jl,	Ks/Jl,			Ds/(Jl*Rg)			  ;
		-1,				0,				1.0/Rg				  ;
		Ds/(Jm*Rg),		-Ks/(Jm*Rg),	-(Ds/(Rg*Rg) + Dm)/Jm ];
Bc = [
		0		, -1.0/Jl;
		0		, 0	     ;
		Kt/Jm	, 0		 ];
C = [1, 0, 0];
D = [0, 0];
sys3 = ss(Ac, Bc, C, D);

sys3d = c2d(sys3, Ts);	% ���̂܂ܗ��U������ꍇ
Ad = sys3d.a
bd = sys3d.b

sys3 = balreal(sys3);	% ���t�����Ă��痣�U������ꍇ
sys3d = c2d(sys3, Ts);
Ad = sys3d.a
bd = sys3d.b




% ARCS-Control�̃e�X�g�p�X�N���v�g
% Yokokura, Yuki 2022/07/27
clc;
clear;
format long;

% �v�����g��ԋ�ԃ��f��(�K��)
A = [-1, -2, -9 ;
      0, -1, -7 ;
     -2,  6, -1 ]
b = [1 ; 0 ; 0]
c = [1, 0, 0]
d = [0];
sys = ss(A,b,c,d);


disp '�� �A�����A�v�m�t�������̉��̃e�X�g'
Q = eye(3);
X = lyap(A,Q)
O = A*X + X*A.' + Q	% �[���ɂȂ邩�m�F


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
Ts = 100e-6; % [s]      �T���v�����O����
Jm = 1e-4;   % [kgm^2]  ���[�^����
Jl = 1;      % [kgm^2]  ���ב�����
Dm = 0.1;    % [Nms/rad]���[�^�S��
Dl = 0.3;    % [Nms/rad]���ב��S��
Ks = 500;    % [Nm/rad] 2�����Ԃ̂΂˒萔
Rg = 50;     %          ������
Kt = 0.04;   % [Nm/A]   �g���N�萔
Ac = [
		-Dm/Jm,      0, -Ks/(Rg*Jm) ;
			 0, -Dl/Jl,       Ks/Jl ;
		1.0/Rg,     -1,           0
];
Bc = [
		Kt/Jm,       0 ;
			0, -1.0/Jl ;
			0,       0
];
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



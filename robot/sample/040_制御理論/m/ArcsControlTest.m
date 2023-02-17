% ARCS-Controlのテスト用スクリプト
% Yokokura, Yuki 2022/07/27
clc;
clear;
format long;

% プラント状態空間モデル(適当)
A = [-1, -2, -9 ;
      0, -1, -7 ;
     -2,  6, -1 ]
b = [1 ; 0 ; 0]
c = [1, 0, 0]
d = [0];
sys = ss(A,b,c,d);


disp '◆ 連続リアプノフ方程式の解のテスト'
Q = eye(3);
X = lyap(A,Q)
O = A*X + X*A.' + Q	% ゼロになるか確認


disp '◆ 可制御/可観測グラミアンのテスト'
Wc = gram(sys,'c')
Wo = gram(sys,'o')


disp '◆ 平衡実現のテスト'
[sys2,g,T,Ti] = balreal(sys);
Ah = sys2.a
bh = sys2.b
ch = sys2.c
Wc2 = gram(sys2,'c')
Wo2 = gram(sys2,'o')


disp '◆ 離散化のテスト'
Ts = 100e-6; % [s]      サンプリング時間
Jm = 1e-4;   % [kgm^2]  モータ慣性
Jl = 1;      % [kgm^2]  負荷側慣性
Dm = 0.1;    % [Nms/rad]モータ粘性
Dl = 0.3;    % [Nms/rad]負荷側粘性
Ks = 500;    % [Nm/rad] 2慣性間のばね定数
Rg = 50;     %          減速比
Kt = 0.04;   % [Nm/A]   トルク定数
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

sys3d = c2d(sys3, Ts);	% そのまま離散化する場合
Ad = sys3d.a
bd = sys3d.b

sys3 = balreal(sys3);	% 平衡化してから離散化する場合
sys3d = c2d(sys3, Ts);
Ad = sys3d.a
bd = sys3d.b



% ARCS-Controlのテスト用スクリプト
% Yokokura, Yuki 2024/07/23
clc;
clear;
%format short;
format longE;

% プラント状態空間モデル(適当)
Ap = [-1, -2, -9 ;
      0, -1, -7 ;
     -2,  6, -1 ]
bp = [1 ; 0 ; 0]
cp = [1, 0, 0]
dp = [0];
sys = ss(Ap, bp, cp ,dp);

disp '◆ 連続リアプノフ方程式の解のテスト'
Qp = eye(3);
Xp = lyap(Ap, Qp)
O = norm(Ap*Xp + Xp*Ap' + Qp)	% ゼロになるか確認

disp '◆ 離散リアプノフ方程式の解のテスト'
Xpd = dlyap(Ap, Qp)


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
Ts = 100e-6;% [s]      サンプリング時間
Jl = 1;		% [kgm^2]  負荷側慣性
Dl = 0.3;	% [Nms/rad]負荷側粘性
Ds = 1;		% [Nms/rad]ねじれ粘性
Ks = 500;	% [Nm/rad] 2慣性間のばね定数
Jm = 1e-4;	% [kgm^2]  モータ慣性
Dm = 0.1;	% [Nms/rad]モータ粘性
Rg = 100;	%          減速比
Kt = 0.04;	% [Nm/A]   トルク定数
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
sys3d = c2d(sys3, Ts);	% そのまま離散化する場合
Ad = sys3d.a
bd = sys3d.b
Ob = obsv(Atc, Ct2)

%{

sys3 = balreal(sys3);	% 平衡化してから離散化する場合
sys3d = c2d(sys3, Ts);
Ad = sys3d.a
bd = sys3d.b
%}



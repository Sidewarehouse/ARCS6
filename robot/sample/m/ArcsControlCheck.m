% ARCS-Controlのテスト用スクリプト
% Yokokura, Yuki 2024/07/23
clc;
clear;
format short;
%format longE;

% プラント状態空間モデル(適当)
Ap = [-1, -2, -9 ;
      0, -1, -7 ;
     -2,  6, -1 ]
bp = [1 ; 0 ; 0]
cp = [1, 0, 0]
dp = [0];
sys = ss(Ap, bp, cp ,dp);

fprintf('\n');
disp '◆ 連続リアプノフ方程式の解のテスト'
Qp = eye(3);
Xp = lyap(Ap, Qp)
O = norm(Ap*Xp + Xp*Ap' + Qp)	% ゼロになるか確認

fprintf('\n');
disp '◆ 離散リアプノフ方程式の解のテスト'
Xpd = dlyap(Ap, Qp)

fprintf('\n');
disp '◆ 可制御/可観測グラミアンのテスト'
Wc = gram(sys,'c')
Wo = gram(sys,'o')

fprintf('\n');
disp '◆ 可制御性行列/可観測性行列のテスト'
Uc = ctrb(Ap, bp)
Uo = obsv(Ap, cp)

fprintf('\n');
disp '◆ 平衡実現のテスト'
[sys2,g,T,Ti] = balreal(sys);
Ah = sys2.a
bh = sys2.b
ch = sys2.c
Wc2 = gram(sys2,'c')
Wo2 = gram(sys2,'o')

fprintf('\n');
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

fprintf('\n');
disp '◆ 離散系状態空間モデル'
load('../040_制御理論/ArcsControlTest.mat');
A = [-0.2516 -0.1684;2.784 0.3549];
B = [0;3];
C = [0 1];
D = 0;
Ts = 1;
sys1 = ss(A,B,C,D,Ts);
figure(1);
	step(sys1, k1);
	hold on;
	plot(k1, y1, 'ro');
	plot(k1, y1n, 'gx');
	hold off;
	ylabel('Output y1 [*]');
	xlabel('Time Step k1 [-]');
	grid on;
	legend('MATLAB','DiscStateSpace::GetOutput','DiscStateSpace::GetNextOutput');
	title('MATLAB と ArcsControl の比較 (DiscStateSpaceクラス)');

fprintf('\n');
disp '◆ 離散系伝達関数'
figure(5);
	sys5 = tf([2, 0],[4, 0, 3, -1], Ts);
	step(sys5, k5);
	hold on;
	stairs(k5, y5, 'ro');
	stairs(k5, y5n, 'gx');
	hold off;
	ylabel('Output y5 [*]');
	xlabel('Time Step k5 [-]');
	grid on;
	legend('MATLAB','DiscTransFunc::GetOutput','DiscTransFunc::GetNextOutput','Location','SouthEast');
	title('MATLAB と ArcsControl の比較 (DiscTransFuncクラス)');
figure(6);
	sys6 = tf([2, 1, 0],[4, 0, 3], Ts);
	step(sys6, k6);
	hold on;
	stairs(k6, y6, 'ro');
	hold off;
	ylabel('Output y6 [*]');
	xlabel('Time Step k6 [-]');
	grid on;
	legend('MATLAB','DiscTransFunc::GetOutput','Location','SouthEast');
	title('MATLAB と ArcsControl の比較 (DiscTransFuncクラス)');

fprintf('\n');
disp '◆ 連続系状態空間モデル'
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
	title('MATLAB と ArcsControl の比較 (StateSpaceクラス)');
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
figure(21);
	A2b = [
		   0,    1,    0 ;
		   0,    0,    1 ;
		-2.5,   -2, -1.5
	];
	b2b = [
		0 ;
		0 ;
		1
	];
	c2b = [ -1.5, -1, -0.5 ];
	d2b = [ 2 ];
	sys2b = ss(A2b, b2b, c2b, d2b);
	[sys2bb, g, Tl, Tr]  = balreal(sys2b)
	t_ = linspace(0,10,100);
	y_ = step(sys2b, t_);
	plot(t_, y_);
	hold on;
	stairs(t2b, y2b, 'r');
	hold off;
	ylabel('Output y2b [*]');
	xlabel('Time [s]');
	grid on;
	legend('MATLAB','StateSpace::GetOutput','Location','SouthEast');
	%figure(100); pzmap(sys2b); grid on;
figure(22);
	A2c = [
		-3.4182, 0 ;
		 1     , 0 
	];
	B2c = [
		4.4545e3, -9.0909e3 ;
		0       ,  0
	];
	C2c = [
		1, 0 ;
		0, 1
	];
	D2c = zeros(2,2);
	sys2c = ss(A2c, B2c, C2c, D2c);
	t_ = linspace(0,10,100);
	y_ = step(sys2c, t_);
subplot(2,1,1)
	plot(t_, y_(:,1));
	hold on;
	stairs(t2c, y2c(1,:), 'r');
	hold off;
	ylabel('Output y2c y1 [*]');
	grid on;
subplot(2,1,2)
	plot(t_, y_(:,2));
	hold on;
	stairs(t2c, y2c(2,:), 'r');
	hold off;
	ylabel('Output y2c y2 [*]');
	xlabel('Time [s]');
	grid on;
	legend('MATLAB','StateSpace::GetOutput','Location','SouthEast');
	
return;

fprintf('\n');
disp '◆ 連続系伝達関数'
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
	title('MATLAB と ArcsControl の比較 (TransFuncクラス)');
figure(4);
	sys4 = tf([4, 5, 6, 7],[2, 3, 4, 5]);
	t_ = linspace(0,10,1000);
	y_ = step(sys4, t_);
	plot(t_, y_);
	hold on;
	stairs(t4, y4, 'r');
	hold off;
	ylabel('Output y4 [*]');
	xlabel('Time [s]');
	grid on;
	legend('MATLAB','TransFunc::GetOutput','Location','SouthEast');
	title('MATLAB と ArcsControl の比較 (TransFuncクラス)');


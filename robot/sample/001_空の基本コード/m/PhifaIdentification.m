%
% αβ軸上の最大鎖交磁束Φfa同定アルゴリズム
%						2019/09/18 Yuki YOKOKURA
clc;
clear;

% CSVファイル名設定
FileName = '../DATA.csv';

% パラメータ設定
Ts   = 250e-6;	% [s] サンプリング時間
P    = 4;		% [-] 極対数
Ra   = 0.671;	% [Ω] 電機子抵抗
La   = 1.40e-3;	% [H]  電機子インダクタンス

% CSVファイルから変数値読み込み
CsvData  = csvread(FileName);
t        = CsvData(:,1);
theta_rm = CsvData(:,2);
Ia       = CsvData(:,3);
Ib       = CsvData(:,4);
Va       = CsvData(:,5);
Vb       = CsvData(:,6);
clear CsvData;
tlen = length(t);

% 電気角速度の計算
g = 6280;	% [rad/s] 推定帯域
s = tf('s');
Gpd = s*g/(s + g);
tsim = (0:Ts:(tlen - 1)*Ts).';
omega_re = lsim(Gpd, theta_rm - theta_rm(1), tsim)*P;


% 生データの表示
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

% 誘起電圧推定と振幅計算
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

% 誘起電圧波形と振幅波形の表示
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

% 最小二乗法による最大鎖交磁束の同定
x = abs(omega_re);
y = Emax;
A(1:length(x),1) = x;
%A(1:length(x),2) = 1; % 電圧オフセットも同定するときのみコメントアウト
u = A\y;
Phifa = u(1)	% [V/(rad/s)] 最大鎖交磁束
%b = u(2)		% [V]		  電圧オフセット
Wrefit  = linspace(0, max(abs(omega_re)));
Emaxfit = Phifa*Wrefit;
figure(3);
	h = plot(abs(omega_re), Emax, 'k+', Wrefit, Emaxfit, 'r');
	set(h,'LineWidth',2);
	xlabel('Elec. Anglular Velocity [rad/s]');
	ylabel('Back EMF Amplitude [V]');
    title(strcat('Max. Flux Linkage Phifa = ', sprintf(' %e', Phifa),' [V/(rad/s)]'));
	grid on;


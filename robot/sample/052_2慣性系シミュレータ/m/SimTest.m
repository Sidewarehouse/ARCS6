% TwoInertiaSimulatorクラスのテスト用スクリプト
% 2023/10/27 Yokokura, Yuki
clc;
clear;

% CSVファイル名設定
FileName = '../DATA.csv';

% 2慣性系パラメータ
Jl = 0.1;	% [kgm^2]		負荷側慣性
Dl = 1e-12;	% [Nm/(rad/s)]	負荷側粘性
Ds = 1;		% [Nm/(rad/s)]	ねじれ粘性
Ks = 10000;	% [Nm/rad]		ねじれ剛性
Jm = 1e-5;	% [kgm^2]		モータ側慣性
Dm = 1e-4;	% [Nm/(rad/s)]	モータ側粘性
Rg = 100;	% [-]			減速比
Kt = 0.1;	% [Nm/A]		トルク定数

% CSVファイルから変数値読み込み
CsvData  = csvread(FileName);
t    = CsvData(:,1);
iq   = CsvData(:,2);
tdis = CsvData(:,3);
wl   = CsvData(:,4);
ths  = CsvData(:,5);
wm   = CsvData(:,6);
taus = CsvData(:,7);
clear CsvData;
tlen = length(t);

% テスト対象の状態空間モデル
% 2慣性共振系の状態空間モデルでの表現
% xp = [wl, ths, wm]^T
Ap = [
	-(Dl + Ds)/Jl,	Ks/Jl,			Ds/(Jl*Rg)			;
	-1,				0,				1.0/Rg				;
	Ds/(Jm*Rg),		-Ks/(Jm*Rg),	-(Ds/Rg^2 + Dm)/Jm
];
bp = [
	0		, -1/Jl ;
	0		,  0    ;
	Kt/Jm	,  0
]
Cp = eye(3);
dp = zeros(3,2);
System1 = ss(Ap,bp,Cp,dp)

% 真値応答の計算
Ts = t(tlen)/(tlen - 1);
ti = 0:Ts:t(tlen);
iqi   = 0.5*SquareWave(0.1, 0, ti, 0.5);
tdisi =   2*SquareWave(0.2, 0, ti, 2.0);
ui = [iqi ; tdisi];
yi = lsim(System1, ui, ti);
wli  = yi(:,1);
thsi = yi(:,2);
wmi  = yi(:,3);
tausi = -Ds*wli + Ks*thsi + Ds/Rg*wmi;

% グラフ描画
figure(1);
clf;
set(gcf,'PaperPositionMode','manual');
set(gcf,'color',[1 1 1]);
subplot(6,1,1);
	h=stairs(ti, iqi, 'k');
		set(h,'linewidth',3);
	hold on;
	h=stairs(t, iq, 'r');
		set(h,'linewidth',1);
	hold off;
	xlabel('Time [s]','FontSize',12);
	ylabel('q-axis Current [A]','FontSize',12);
	set(gca,'FontSize',12);
	legend('MATLAB','ARCS','Location','SouthEast','Orientation','Vertical');
	grid on;
subplot(6,1,2);
	h=stairs(ti, tdisi, 'k');
		set(h,'linewidth',3);
	hold on;
	h=stairs(t, tdis, 'r');
		set(h,'linewidth',1);
	hold off;
	xlabel('Time [s]','FontSize',12);
	ylabel({'Disturbance','Torque [Nm]'},'FontSize',12);
	set(gca,'FontSize',12);
	grid on;
subplot(6,1,3);
	h=stairs(ti, wli, 'k');
		set(h,'linewidth',3);
	hold on;
	h=stairs(t, wl, 'r');
		set(h,'linewidth',1);
	hold off;
	xlabel('Time [s]','FontSize',12);
	ylabel({'Load-side','Velocity [rad/s]'},'FontSize',12);
	set(gca,'FontSize',12);
	grid on;
subplot(6,1,4);
	h=stairs(ti, thsi, 'k');
		set(h,'linewidth',3);
	hold on;
	h=stairs(t, ths, 'r');
		set(h,'linewidth',1);
	hold off;
	xlabel('Time [s]','FontSize',12);
	ylabel({'Torsion Angle [rad]'},'FontSize',12);
	set(gca,'FontSize',12);
	grid on;
subplot(6,1,5);
	h=stairs(ti, wmi, 'k');
		set(h,'linewidth',3);
	hold on;
	h=stairs(t, wm, 'r');
		set(h,'linewidth',1);
	hold off;
	xlabel('Time [s]','FontSize',12);
	ylabel({'Motor-side','Velocity [rad/s]'},'FontSize',12);
	set(gca,'FontSize',12);
	grid on;
subplot(6,1,6);
	h=stairs(ti, tausi, 'k');
		set(h,'linewidth',3);
	hold on;
	h=stairs(t, taus, 'r');
		set(h,'linewidth',1);
	hold off;
	xlabel('Time [s]','FontSize',12);
	ylabel({'Torsion','Torque [Nm]'},'FontSize',12);
	set(gca,'FontSize',12);
	grid on;

figure(2);
set(gcf,'PaperPositionMode','manual');
set(gcf,'color',[1 1 1]);
	h=stairs(t, ti, 'k');
		set(h,'linewidth',1);
	xlabel('Ideal Time [s]','FontSize',12);
	ylabel('Actual Time [s]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;


% EPSファイル生成(ローカルで実行のこと)
% print(gcf,'-depsc2','-tiff',strcat(FileName,'.eps'));

% 方形波生成関数
function y = SquareWave(freq, phase, time, starttime)
	tlen = length(time);
	r(1:tlen) = 0;
	y(1:tlen) = 0;
	
	% 開始時間だけシフトされた正弦波を生成
	for i = 1:tlen
		if time(i) < starttime
			r(i) = 0;			% 開始時刻より前のときはゼロ出力
		else
			r(i) = sin(2*pi*freq*(time(i) - starttime) + phase);	% 開始時刻以降は方形波出力
			% 方形波に変換
			if 0 <= r(i)
				y(i) =  1;
			else
				y(i) = -1;
			end
		end
	end
end



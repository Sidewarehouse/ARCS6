% ARCS-FRA用の周波数特性測定スクリプト V2021A
% Yokokura, Yuki & Muto Hirotaka & Shunsuke Suzuki 2021/12/22
%
% 参考文献： 水流 徹, 春山 志郎:「交流インピーダンス測定法の選択 -FRAとFFT-」防食技術, vol. 35, no. 5, Jan. 1986.
%
% 使い方：
% 1. FRAgeneratorクラスを用いて実機を励起，測定データをDATA.csvに保存
%              ↓
% 2. このスクリプト内のFRAの測定パラメータをFRAgeneratorの設定と合わせる
%              ↓
% 3. このスクリプトを実行するとFRAの計算が行われ周波数特性が得られる
clc;
clear all;

%% FRAの測定パラメータの設定(実機のパラメータに合わせる)
FileName = '../DATA.csv';	% CSVファイル名
Ts = 100e-6;% [s]  サンプリング時間
Fsta =   3;	% [Hz] 開始周波数
Fend = 100;	% [Hz] 終了周波数
Fstep = 1;	% [Hz] 周波数ステップ
Ni = 10;	% [-]  積分周期
Au = 0.1;	% [-]  振幅
Bu = -0.45;	% [-]  バイアス
Tsta = 0;	% [s] FRA開始時刻

%% 周波数ベクトル生成
freq = Fsta:Fstep:Fend;
flen = length(freq);

%% 実機データ読み込み
data = csvread(FileName);
t = data(:,1);     % [s]  時間
f = data(:,2);     % [Hz] 周波数
iqref = data(:,3); % [A]  電流指令値
wmres = data(:,4); % [rad]速度応答
clear data;

%% バイアス補正
iqref = iqref - Bu;

%% 開始時刻より以前のデータを削除する前処理
t = t - Tsta;			% 開始時刻分だけ左にシフト
t(1:Tsta/Ts) = [];		% 開始時刻より以前は削除
f(1:Tsta/Ts) = [];		% 開始時刻より以前は削除
iqref(1:Tsta/Ts) = [];	% 開始時刻より以前は削除
wmres(1:Tsta/Ts) = [];	% 開始時刻より以前は削除
tlen = length(t);

%% 生データの表示
figure(1); set(gcf,'color',[1 1 1]);
subplot(3,1,1); plot(t,f);      %周波数の変化
	xlabel('Time [s]','FontSize',14,'FontName','Times New Roman')
	ylabel('Frequency [rad/s]','FontSize',14,'FontName','Times New Roman')
	set(gca,'FontSize',14);
	grid on;
subplot(3,1,2); plot(t,iqref);  %入力信号
	xlabel('Time [s]','FontSize',14,'FontName','Times New Roman')
	ylabel('Current [A]','FontSize',14,'FontName','Times New Roman')
	set(gca,'FontSize',14);
	grid on;
subplot(3,1,3); plot(t,wmres);  %出力信号
	xlabel('Time [s]','FontSize',14,'FontName','Times New Roman')
	ylabel('Velocity [rad/s]','FontSize',14,'FontName','Times New Roman')
	set(gca,'FontSize',14);
	grid on;

%% FRAの計算
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
%figure(127); plot(freq,Ar,'x-', freq,Ai,'x-');	% 確認用

%% 最後のデータだけおかしくなるので消しておく
Ar(flen) = [];
Ai(flen) = [];
freq(flen) = [];
flen = length(freq);

%% 周波数特性の計算
Ay = sqrt(Ar.^2 + Ai.^2);          % 出力振幅計算
G = 20*log10(Ay./Au);              % ゲイン特性計算
P = unwrap(-atan2(Ai,Ar))*180/pi;  % 位相特性計算

%% 周波数特性の描画
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

% EPSファイル生成(ローカルで実行のこと)
% print(gcf,'-depsc2','-tiff','FRA.eps');

% データ保存
save(strcat(FileName, '_FRA.mat'), 'Ts', 'Fsta', 'Fend', 'Fstep', 'Ni', 'Au', 'Bu', 'Tsta', 'freq', 'G', 'P');


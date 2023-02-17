%
% dq軸の抵抗とインダクタンスの同定アルゴリズム
%						2019/09/18 Yuki YOKOKURA
clc;
clear;

% CSVファイル名設定
FileName = '../DATA.csv';

% パラメータ設定
Ts       = 250e-6;	% [s] サンプリング時間
Tfcs_sta = 5.5;		% [s] 計測に使うデータの開始時刻
Tfcs_end = 6.9;		% [s] 計測に使うデータの終了時刻
Tstep    = 6.0;		% [s] ステップ入力開始時刻
Tstd_sta = 6.5;	% [s] 定常状態開始時刻
Tstd_end = 6.9;	% [s] 定常状態終了時刻

% CSVファイルから変数値読み込み
CsvData  = csvread(FileName);
t   = CsvData(:,1);
Vdq = CsvData(:,2);
Idq = CsvData(:,3);
clear CsvData;
tlen = length(t);

% 生データの表示
figure(1);
clf;
subplot(2,1,1);
	h = plot(t, Vdq);
	set(h,'LineWidth',2);
	xlabel('Time [s]');
	ylabel('Voltage Vd or Vq [V]');
	grid on;
subplot(2,1,2);
	h = plot(t, Idq);
	set(h,'LineWidth',2);
	xlabel('Time [s]');
	ylabel('Current Id or Iq [A]');
	grid on;

% データ切り出し，定常電流値の計測と63.2％の計算
Iastd = mean(Idq((Tstd_sta/Ts):(Tstd_end/Ts)));
Ia632 = Iastd*(1 - exp(-1));
t = t((Tfcs_sta/Ts):(Tfcs_end/Ts)) - Tstep;
Ia = Idq((Tfcs_sta/Ts):(Tfcs_end/Ts));

% 63.2％の場所の探索
ep = 0.01;	% 探索の許容幅
idx = round(mean(find( (Ia632 - ep) < Ia & Ia < (Ia632 + ep) ))); % 応答が63.2％ぐらいに入る要素番号を取得
if(isnan(idx) == true)
	warning('eqが小さすぎて63.2％付近のデータが見つかりません。epを大きくして下さい。');
else
	tau = t(idx);   % 時定数を取得
end;

% 抵抗とインダクタンスの計算
Rdq = max(Vdq)/Iastd;
Ldq = Rdq*tau;

% 計測データ表示
figure(2);
	h = plot(t, Ia, [t(1) t(length(t))], [Iastd Iastd], [t(1) t(length(t))], [Ia632 Ia632], [0 0], [0 Iastd], [tau tau], [0 Iastd]);
	set(h,'LineWidth',2);
	xlabel('Time [s]');
	ylabel('Current Id or Iq [A]');
    title(strcat('Time constant \tau = ',sprintf(' %e',tau),' [s]  Resistance Rd or Rq = ',sprintf(' %e', Rdq),' [Ohm]  Inductance Ld or Lq = ',sprintf(' %e', Ldq),' [H]'));
	grid on;

% EPSファイル生成
print(gcf,'-depsc2','-tiff','RdqLdqIdentResults.eps');

% 同定値でのシミュレーション
s = tf('s');
P = 1/(Ldq*s + Rdq);
Iasim = lsim(P, Vdq((Tfcs_sta/Ts):(Tfcs_end/Ts)), [0:Ts:((Tfcs_end/Ts) - (Tfcs_sta/Ts))*Ts]);
figure(3);
	h = plot(t, Ia,'k', t, Iasim,'r');
	set(h,'LineWidth',2);
	xlabel('Time [s]');
	ylabel('Current Id or Iq [A]');
	legend('Measured','Identified');
	grid on;


% クーロン摩擦と粘性摩擦係数の同定スクリプト
% 2021/12/10 Yokokura, Yuki
clc;
clear;

% CSVファイル名設定
FileName = '../FricIdentExp2L-1.csv';

% パラメータ設定
Ktn  = 0.31;	% [Nm/A] トルク定数
Nave = 3000;	% [-] 定常偏差検出のための移動平均点数
wzero = 1e-1;   % [rad/s] 零とみなす速度

% CSVファイルから変数値読み込み
CsvData  = csvread(FileName);
t = CsvData(:,1);
iqref = CsvData(:,2);
wm = CsvData(:,3);
clear CsvData;
tlen = length(t);

% 定常速度を検出するための移動メディアンフィルタ
wm_filt = movmedian(wm, [Nave 0]);

% 定常状態における各値の検出
j = 0;
for i = 2:tlen
	if iqref(i - 1) ~= iqref(i)			% 電流指令値が変化したら，
		j = j + 1;
		t_stead(j)  = t(i - 1);			% 定常状態の時刻
		iq_stead(j) = iqref(i - 1);		% 定常電流
		wm_stead(j) = wm_filt(i - 1);	% 定常速度
	end
end
stead_len = length(wm_stead);

% 速度vsトルク特性への変換
Ktiq_stead = Ktn*iq_stead;			% [Nm] 印加トルクの計算
Ktiq_stead_abs = abs(Ktiq_stead);	% [Nm] 正回転と負回転の特性を同一視する
wm_stead_abs   = abs(wm_stead);		% [rad/s] 正回転と負回転の特性を同一視する

% 速度が零の場合はフィッティング元データから除外する
% (最小二乗法の同定精度を上げるための前処理)
Ktiq_lsq = Ktiq_stead_abs;
wm_lsq = wm_stead_abs;
for i = stead_len:-1:1,
   if(-wzero < wm_lsq(i) && wm_lsq(i) < wzero)	% 速度が零なら，
       Ktiq_lsq(i) = [];	% データ削除
       wm_lsq(i) = [];		% データ削除
   end
end
lsq_len = length(wm_lsq);

% 最小二乗法による摩擦同定
A(1:lsq_len,1) = wm_lsq.';
A(1:lsq_len,2) = 1;
b = Ktiq_lsq.';
x = A\b;						% バックスラッシュ演算子で A*x = b をxについて解く
Dall  = x(1);					% [Nm/(rad/s)] 粘性摩擦係数(1次関数の傾き)
Fclmb = x(2);					% [Nm] クーロン摩擦力(1次関数の切片)
Ktiq_fit = wm_lsq*Dall + Fclmb;	% [Nm] 摩擦トルクフィッティング直線の計算
fprintf('Dall  = %3.2e;	// [Nm/(rad/s)] 全粘性摩擦係数\n', Dall);
fprintf('Fclmb = %3.2e;	// [Nm] クーロン摩擦力\n', Fclmb);

% グラフ描画
figure(1);
clf;
set(gcf,'PaperPositionMode','auto');
set(gcf,'color',[1 1 1]);
subplot(2,2,1);
	h=plot(t, iqref, 'k');
		set(h,'linewidth',2);
	xlabel({'Time [s]','(a)'},'FontSize',12);
	ylabel('q-axis current [A]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
subplot(2,2,3);
	h=plot(t, wm, 'k', t, wm_filt, 'r', t_stead, wm_stead, 'go');
		set(h,'linewidth',2);
	xlabel({'Time [s]','(b)'},'FontSize',12);
	ylabel('Velocity [rad/s]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	legend('Raw','Moving Median','Steady State','Location','NorthWest','Orientation','Vertical');
subplot(2,2,[2,4]);
	h=plot(wm_stead_abs, Ktiq_stead_abs, 'go', wm_lsq, Ktiq_lsq, 'k+', wm_lsq, Ktiq_fit, 'r');
		set(h,'linewidth',2);
	xlabel({'Velocity [rad/s]','(c)'},'FontSize',12);
	ylabel('Motor Torque [Nm]','FontSize',12);
	set(gca,'FontSize',12);
	grid on;
	legend('Steady State','Non-zero Speed','Fitted by Least Sq.','Location','SouthEast','Orientation','Vertical');
	title( strcat('Dall = ', sprintf('%4.3e',Dall), ' [Nm/(rad/s)]', '   Fclmb = ', sprintf('%4.3e',Fclmb), ' [Nm]') );
saveas(gcf, strcat(FileName,'.png'));

% EPSファイル生成(ローカルで実行のこと)
% print(gcf,'-depsc2','-tiff',strcat(FileName,'.eps'));

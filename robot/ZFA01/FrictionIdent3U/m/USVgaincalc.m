% USVユニサーボ電流制御ゲイン計算スクリプト
% Yokokura, Yuki - 2021/12/09
clc;
clear;

% プラントパラメータ設定
Rn  = 0.375;	% [Ω]  抵抗
Ln  = 2.11e-3;	% [H]   インダクタンス

% 制御器パラメータ設定
T  = 1/10e3;	% [s] 制御周期
Wc = 2*pi*300;	% [rad/s] 電流制御帯域
Zc = 1;			% [-] 制動係数

% 連続系の極から離散系の極への変換
a = -Zc*Wc + Wc*sqrt(Zc^2 - 1)	% 連続系S平面上の極1
b = -Zc*Wc - Wc*sqrt(Zc^2 - 1)	% 連続系S平面上の極2
alpha = exp(a*T)	% 離散系Z平面上の極1
beta  = exp(b*T)	% 離散系Z平面上の極2

% 離散系I-P電流制御器のゲイン計算
Ki = Rn/T*(1 - alpha - beta + alpha*beta)/(1 - exp(-Rn/Ln*T))
Kp = Rn*(exp(-Rn/Ln*T) - alpha*beta)/(1 - exp(-Rn/Ln*T))


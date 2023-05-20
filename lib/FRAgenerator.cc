//! @file FRAgenerator.hh
//! @brief FRA用信号生成器
//!
//! Frequency Response Analysis のための信号生成器
//!
//! @date 2022/03/31
//! @author Yokokura, Yuki & Muto Hirotaka
//
// Copyright (C) 2011-2022 Yokokura, Yuki & Muto Hirotaka
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include <cmath>
#include "FRAgenerator.hh"
#include "ARCSeventlog.hh"

using namespace ARCS;

//! @brief コンストラクタ
FRAgenerator::FRAgenerator(
	const double FreqMin,		// [Hz]開始周波数
	const double FreqMax,	 	// [Hz]終了周波数
	const double FreqStep,		// [Hz]周波数ステップ
	const double NumIntg,		// [-] 積分周期 (1周波数につき何回sin波を入力するか)
	const double Ampl,			// [-] 振幅
	const double Bias,			// [-] バイアス
	const double TimeSta		// [s] FRA開始時刻
)
	: Fmin(FreqMin), Fmax(FreqMax), Fstep(FreqStep), Ni(NumIntg), Au(Ampl), Bu(Bias), Tsta(TimeSta),
	  isEnd(false), f(Fmin), tini(0), tauclmb(0)
{
	PassedLog();
}

//! @brief デストラクタ
FRAgenerator::~FRAgenerator(){
	PassedLog();
}

//! @brief FRA信号出力関数(引数で返す版)
//! @param[in] Time 	[s]  時刻
//! @param[out]	Freq	[Hz] 周波数
//! @param[out]	Output	[*]  FRA信号出力
void FRAgenerator::GetSignal(const double Time, double& Freq, double& Output){
	const double t = Time - Tsta;			// [s] FRA開始時刻補正後の時刻
	
	if(Tsta <= Time && isEnd == false){		// FRA開始時刻を過ぎて，且つ最大周波数以下であれば，
		Output = Au*cos(2.0*M_PI*f*(t - tini)) + Bu;	 // [A] FRA信号を生成
		
		if (Ni/f <= t - tini){ 	// Ni回必要な時間になるまでf[Hz]を加える。Ni回終わったら....
			if(f <= Fmax){		// 最大周波数以下なら次の周波数へ移行
				tini = t;		// [s] 次の周波数の初期時刻
				f = f + Fstep;	// [Hz]次の周波数 (double型の累積加算は…まあご愛嬌)
			}else{				
				isEnd = true;	// 最大周波数に達したら終了
			}
		}
	}else{
		Output = Bu;	// FRA非動作時はバイアス値を出力
	}
	
	Freq = f;		// [Hz] 周波数を出力
}

//! @brief FRA信号出力関数(周波数を引数で返し，出力を戻り値で返す版)
//! @param[in] Time 	[s]  時刻
//! @param[out]	Freq	[Hz] 周波数
//! @return	Output	[*]  FRA信号出力
double FRAgenerator::GetSignal(const double Time, double& Freq){
	double Output;
	GetSignal(Time, Freq, Output);
	return Output;
}

//! @brief FRA信号出力関数(タプルで返す版)
//! @param[in] Time [s] 時刻
//! @return {f [Hz] 周波数, outsig [*] FRA信号} のタプル
std::tuple<double, double> FRAgenerator::GetSignal(const double Time){
	double Freq, Output;
	GetSignal(Time, Freq, Output);
	return {Freq, Output};	// 周波数とFRA信号をタプルで返す
}

//! @brief FRA信号出力関数(クーロン摩擦補償付き)
//! @param[in]	Time 		[s]  時刻
//! @param[in]	Velocity	[rad/s] モータ速度
//! @param[out]	Freq		[Hz] 周波数
//! @param[out] ObsrvOut	[*] 観測用FRA信号出力(クーロン摩擦トルク補償前のFRA信号)
//! @param[out]	DriveOut	[*] 駆動用FRA信号出力(クーロン摩擦トルク補償後のFRA信号)
void FRAgenerator::GetSignalWithFricComp(const double Time, const double Velocity, double& Freq, double& ObsrvOut, double& DriveOut){
	double Output;
	GetSignal(Time, Freq, Output);
	ObsrvOut = Output;
	if(0 <= Velocity){
		// モータ速度が正なら，
		DriveOut = Output + tauclmb;
	}else{
		// モータ速度が負なら，
		DriveOut = Output - tauclmb;
	}
}

//! @brief クーロン摩擦トルクを設定する関数
//! @param[in]	ClmbFric	[Nm] クーロン摩擦トルク
void FRAgenerator::SetCoulombFriction(const double ClmbFric){
	tauclmb = ClmbFric;
}

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

#ifndef FRAGENERATOR
#define FRAGENERATOR

#include <tuple>

namespace ARCS {	// ARCS名前空間
//! @brief FRA用信号生成器
class FRAgenerator {
	public:
		FRAgenerator(
			const double FreqMin,		// [Hz]開始周波数
			const double FreqMax,	 	// [Hz]終了周波数
			const double FreqStep,		// [Hz]周波数ステップ
			const double NumIntg,		// [-] 積分周期 (1周波数につき何回sin波を入力するか)
			const double Ampl,			// [-] 振幅
			const double Bias,			// [-] バイアス
			const double TimeSta		// [s] FRA開始時刻
		);				//!< コンストラクタ
		~FRAgenerator();//!< デストラクタ
		void GetSignal(const double Time, double& Freq, double& Output);//!< FRA信号出力関数(引数で返す版)
		double GetSignal(const double Time, double& Freq);				//!< FRA信号出力関数(周波数を引数で返し，出力を戻り値で返す版)
		std::tuple<double, double> GetSignal(const double Time);		//!< FRA信号出力関数(タプルで返す版)
		void GetSignalWithFricComp(const double Time, const double Velocity, double& Freq, double& ObsrvOut, double& DriveOut);	//!< FRA信号出力関数(クーロン摩擦補償付き)
		void SetCoulombFriction(const double ClmbFric);					//!< クーロン摩擦トルクを設定する関数
		
	private:
		FRAgenerator(const FRAgenerator&) = delete;						//!< コピーコンストラクタ使用禁止
		const FRAgenerator& operator=(const FRAgenerator&) = delete;	//!< 代入演算子使用禁止
		const double Fmin;	//!< [Hz]開始周波数
		const double Fmax;	//!< [Hz]終了周波数
		const double Fstep;	//!< [Hz]周波数ステップ
		const double Ni;	//!< [-] 積分周期  (1周波数につき何回sin波を入力するか)
		const double Au;	//!< [-] 振幅
		const double Bu;	//!< [-] バイアス
		const double Tsta;	//!< [s] FRA開始時刻
		bool isEnd;			//!< FRA終了フラグ
		double f;			//!< [Hz]現在の周波数	
		double tini;		//!< [s] 各周波数ごとの初期時刻
		double tauclmb;		//!< [Nm] クーロン摩擦トルク
};
}

#endif


//! @file TwoStepWave.hh
//! @brief 2段ステップ信号生成器
//!
//! 2段ステップの信号を出力する。
//!
//! @date 2019/02/26
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef TWOSTEPWAVE
#define TWOSTEPWAVE

namespace ARCS {	// ARCS名前空間
	//! @brief 2段ステップ信号生成器
	class TwoStepWave {
		public:
			TwoStepWave(
				double Stp1Hgt, double Stp2Hgt,
				double Stp1StaT, double Stp2StaT, double Stp2EndT, double Stp1EndT
			);	// コンストラクタ  ステップ1の高さ, ステップ2の高さ，ステップ1の開始時刻, ステップ2の開始時刻, ステップ2の終了時刻, ステップ1の終了時刻
			~TwoStepWave();				// デストラクタ
			double GetSignal(double t) const;	// 2段ステップの信号を出力する t : 時刻
			
		private:
			TwoStepWave(const TwoStepWave&);					// コピーコンストラクタ使用禁止
			const TwoStepWave& operator=(const TwoStepWave&);	// 代入演算子使用禁止
			const double Height1;	// ステップ1の高さ
			const double Height2;	// ステップ1の高さ
			const double T1sta;	// ステップ1の開始時刻
			const double T2sta;	// ステップ2の開始時刻
			const double T2end;	// ステップ2の終了時刻
			const double T1end;	// ステップ1の終了時刻
	};
}

#endif




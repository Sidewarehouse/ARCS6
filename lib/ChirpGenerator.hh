//! @file ChirpGenerator.hh
//! @brief チャープ信号生成器クラス
//!
//! チャープ信号生成器
//!
//! @date 2019/02/22
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef CHIRPGENERATOR
#define CHIRPGENERATOR

namespace ARCS {	// ARCS名前空間
	//! @brief チャープ信号生成器クラス
	class ChirpGenerator {
		public:
			ChirpGenerator(double freq_min, double freq_max, double time_max);	//!< コンストラクタ
			~ChirpGenerator();	//!< デストラクタ
			double GetSignal(double time) const;	//!< チャープ信号の取得
			double GetFrequency(double time) const;	//!< 現在の周波数 [Hz] の取得
			
		private:
			ChirpGenerator(const ChirpGenerator&) = delete;					//!< コピーコンストラクタ使用禁止
			const ChirpGenerator& operator=(const ChirpGenerator&) = delete;//!< 代入演算子使用禁止
			double fmin;	//!< 開始周波数
			double fmax;	//!< 終了周波数
			double tmax;	//!< 終了時間
			double A;		//!< 係数
	};
}

#endif




//! @file SshapeGenerator.hh
//! @brief S字軌道生成器
//!
//! 入力された位置データからS字軌道を生成するクラス
//!
//! @date 2021/08/02
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2021 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef SSHAPEGENERATOR
#define SSHAPEGENERATOR

#include <cassert>
#include "MovingAverage.hh"

// ARCS組込み用マクロ
#ifdef ARCS_IN
	// ARCSに組み込まれる場合
	#include "ARCSassert.hh"
	#include "ARCSeventlog.hh"
#else
	// ARCSに組み込まれない場合
	#define arcs_assert(a) (assert(a))
	#define PassedLog()
	#define EventLog(a)
	#define EventLogVar(a)
#endif

namespace ARCS {	// ARCS名前空間
//! @brief S字軌道生成器
//! @tparam	N	移動平均の点数
template <size_t N>
class SshapeGenerator {
	public:
		//! @brief コンストラクタ
		SshapeGenerator()
			: RampStage(), SshapeStage()
		{
			PassedLog();
		}

		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		SshapeGenerator(SshapeGenerator&& r)
			// :
		{
			
		}

		//! @brief デストラクタ
		~SshapeGenerator(){
			PassedLog();
		}
		
		//! @brief S字軌道の最新値を取得する関数
		//! @param[in]	input	入力値
		//! @return	S字軌道の値
		double GetShapedSignal(const double input){
			return SshapeStage.GetSignal(RampStage.GetSignal(input));		// 移動平均を2回掛けてS字にする
		}
		
		//! @brief S字軌道の初期値を設定する関数
		//! @param[in]	init	初期値
		void SetInitialValue(const double init){
			RampStage.Fill(init);	// 移動平均のメモリを初期値で埋める
			SshapeStage.Fill(init);	// 移動平均のメモリを初期値で埋める
		}
		
	private:
		SshapeGenerator(const SshapeGenerator&) = delete;					//!< コピーコンストラクタ使用禁止
		const SshapeGenerator& operator=(const SshapeGenerator&) = delete;	//!< 代入演算子使用禁止
		MovingAverage<N> RampStage;		//!< ランプ軌道生成段
		MovingAverage<N> SshapeStage;	//!< S字軌道生成段
		
};
}

#endif


//! @file SshapeGenerators.hh
//! @brief S字軌道生成器(ベクトル版)
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

#ifndef SSHAPEGENERATORS
#define SSHAPEGENERATORS

#include <cassert>
#include <array>
#include "ArcsMatrix.hh"
#include "SshapeGenerator.hh"

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
//! @tparam M	縦ベクトルの高さ
template <size_t N, size_t M>
class SshapeGenerators {
	public:
		//! @brief コンストラクタ
		SshapeGenerators()
			: SshapeGeneratorVec()
		{
			PassedLog();
		}

		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		SshapeGenerators(SshapeGenerators&& r)
			// :
		{
			
		}

		//! @brief デストラクタ
		~SshapeGenerators(){
			PassedLog();
		}
		
		//! @brief S字軌道の最新値を取得する関数(引数で返す版)
		//! @param[in]	input	入力ベクトル
		//! @param[out]	output	S字軌道の最新の出力ベクトル
		void GetShapedSignal(const ArcsMat<M,1>& input, ArcsMat<M,1>& output){
			for(size_t i = 1; i <= M; ++i){
				output[i] = SshapeGeneratorVec[i-1].GetShapedSignal(input[i]);
			}
		}
		
		//! @brief S字軌道の最新値を取得する関数(ベクトルで返す版)
		//! @param[in]	input	入力ベクトル
		//! @return	output	S字軌道の最新の出力ベクトル
		ArcsMat<M,1> GetShapedSignal(const ArcsMat<M,1>& input){
			ArcsMat<M,1> ret;
			GetShapedSignal(input, ret);
			return ret;
		}
		
		//! @brief S字軌道の初期値を設定する関数
		//! @param[in]	init	初期値ベクトル
		void SetInitialValue(const ArcsMat<M,1>& init){
			for(size_t i = 1; i <= M; ++i) SshapeGeneratorVec[i-1].SetInitialValue(init[i]);
		}
		
	private:
		SshapeGenerators(const SshapeGenerators&) = delete;					//!< コピーコンストラクタ使用禁止
		const SshapeGenerators& operator=(const SshapeGenerators&) = delete;//!< 代入演算子使用禁止
		std::array<SshapeGenerator<N>, M> SshapeGeneratorVec;				//!< S字軌道生成器の配列
		
};
}

#endif


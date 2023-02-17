//! @file Shuffle.hh
//! @brief シャッフルクラス
//!
//! メルセンヌ・ツイスタによって行列のランダムシャッフルを行うクラス
//!
//! @date 2020/04/19
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef SHUFFLE
#define SHUFFLE

#include <cassert>
#include <cmath>
#include "Matrix.hh"
#include "RandomGenerator.hh"

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
//! @brief シャッフルクラス
class Shuffle {
	public:
		//! @brief コンストラクタ
		Shuffle()
			: Rand(0,1)	// 0～1の範囲の乱数
		{
			
		}

		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		Shuffle(Shuffle&& r)
			: Rand(std::move(r.Rand))
		{
			
		}

		//! @brief デストラクタ
		~Shuffle(){
			
		}
		
		//! @brief 行列の行をランダムに入れ替える
		//! @param[in,out]	U	入出力行列
		template <size_t N, size_t M>
		void ShuffleMatrixRow(Matrix<N,M>& U){
			// フィッシャー-イェーツのシャッフルアルゴリズム
			size_t j = 0;
			for(size_t i = M - 1; i != 0; --i){
				j = std::floor(Rand.GetDoubleRandom()*i);
				swaprow(U, j+1, i+1);
			}
		}

		//! @brief 行列1と行列2の行を同時にランダムに入れ替える
		//! @param[in,out]	U1	入出力行列1
		//! @param[in,out]	U2	入出力行列2
		template <size_t N, size_t M, size_t L>
		void ShuffleMatrixRow(Matrix<N,M>& U1, Matrix<L,M>& U2){
			// フィッシャー-イェーツのシャッフルアルゴリズム
			size_t j = 0;
			for(size_t i = M - 1; i != 0; --i){
				j = std::floor(Rand.GetDoubleRandom()*i);
				swaprow(U1, j+1, i+1);
				swaprow(U2, j+1, i+1);
			}
		}

	private:
		Shuffle(const Shuffle&) = delete;					//!< コピーコンストラクタ使用禁止
		const Shuffle& operator=(const Shuffle&) = delete;	//!< 代入演算子使用禁止
		RandomGenerator Rand;	//!< メルセンヌ・ツイスタ乱数生成器
};
}

#endif


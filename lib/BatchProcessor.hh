//! @file BatchProcessor.hh
//! @brief バッチ処理器クラス
//!
//! ニューラルネットワークデータセット用のバッチ処理を行うクラス
//!
//! @date 2020/05/12
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef BATCHPROCESSOR
#define BATCHPROCESSOR

#include <cassert>
#include "Matrix.hh"

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
//! @brief バッチ処理器クラス
class BatchProcessor {
	public:
		//! @brief i番目のミニバッチデータを取得する関数
		//! @tparam	N	データセットのチャネル数(データセット行列の幅)
		//! @tparam	D	データセットのデータ数(データセット行列の高さ)
		//! @tparam M	ミニバッチ数(ミニバッチ行列の幅)
		//! @param[in]	i	ミニバッチ番号
		//! @return	切り出したミニバッチデータ行列
		template <size_t N, size_t D, size_t M>
		static Matrix<M,N> GetMiniBatchData(const Matrix<N,D>& U, const size_t i){
			Matrix<N,M> Buff;
			size_t k = 1;
			for(size_t j = (i - 1)*M; j < (i - 1)*M + M; j++){
				setrow(Buff, getrow(U, j+1), k);	// ミニバッチサイズごとに切り出す
				++k;
			}
			return tp(Buff);
		}
		
	private:
		BatchProcessor() = delete;										//!< コンストラクタ使用禁止
		BatchProcessor(BatchProcessor&& r) = delete;					//!< ムーブコンストラクタ使用禁止
		BatchProcessor(const BatchProcessor&) = delete;					//!< コピーコンストラクタ使用禁止
		const BatchProcessor& operator=(const BatchProcessor&) = delete;//!< 代入演算子使用禁止
		~BatchProcessor() = delete;										//!< デストラクタ使用禁止
};
}

#endif


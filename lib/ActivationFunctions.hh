//! @file ActivationFunctions.hh
//! @brief 活性化関数
//!
//! ニューラルネット用の活性化関数
//!
//! @date 2021/08/08
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2021 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef ACTIVATIONFUNCTIONS
#define ACTIVATIONFUNCTIONS

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
//! @brief 活性化関数のタイプの定義
enum class ActvFunc {
	STEP,		//!< ステップ関数
	IDENTITY,	//!< 恒等関数
	SIGMOID,	//!< シグモイド関数
	ReLU,		//!< ReLU関数
	TANH,		//!< tanh関数
	SOFTMAX		//!< ソフトマックス関数
};

//! @brief 活性化関数
class ActivationFunctions {
	public:
		//! @brief ステップ活性化関数
		//! @param[in]	u	入力データ
		//! @return	出力
		static double Step(double u){
			double ret;
			if(0 <= u){
				ret = 1;
			}else{
				ret = 0;
			}
			return ret;
		}
		
		//! @brief 恒等活性化関数
		//! @param[in]	u	入力データ
		//! @return	出力
		static double Identity(double u){
			return u;
		}
		
		//! @brief シグモイド関数活性化関数
		//! @param[in]	u	入力データ
		//! @return	出力
		static double Sigmoid(double u){
			return 1.0/(1.0 + exp(-u));
		}
		
		//! @brief シグモイド関数活性化関数の微分
		//! @param[in]	u	入力データ
		//! @return	出力
		static double DerivativeSigmoid(double u){
			return (1.0 - Sigmoid(u))*Sigmoid(u);
		}
		
		//! @brief ReLU活性化関数
		//! @param[in]	u	入力データ
		//! @return	出力
		static double ReLU(double u){
			double ret;
			if(0 <= u){
				ret = u;
			}else{
				ret = 0;
			}
			return ret;
		}
		
		//! @brief ReLU活性化関数の微分
		//! @param[in]	u	入力データ
		//! @return	出力
		static double DerivativeReLU(double u){
			return Step(u);
		}
		
		//! @brief Softmax関数
		//! @tparam	M	縦ベクトルの長さ
		//! @param[in]	u	入力縦ベクトル
		//! @return	出力
		template <size_t M>
		static Matrix<1,M> Softmax(const Matrix<1,M>& u){
			const Matrix<1,1> SigExp = sumcolumn(expe(u));	// 要素ごとにexpを通して縦方向に加算
			return expe(u)/(SigExp.GetElement(1,1));		// 要素ごとに加算結果で除算して全体に対しての割合を返す
		}
		
		//! @brief 活性化関数
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		//! @tparam T	活性化関数の種類
		//! @tparam	N	行列の幅
		//! @tparam M	行列の高さ
		template <ActvFunc T, size_t N, size_t M>
		static void f(const Matrix<N,M>& U, Matrix<N,M>& Y){
			// 活性化関数の種類の設定に従って呼ぶ関数をコンパイル時に変える
			if constexpr(T == ActvFunc::STEP){
				// ステップ活性化関数の場合
				// 行列の要素ごとに活性化関数を掛ける
				for(size_t n = 1; n <= N; ++n){
					for(size_t m = 1; m <= M; ++m){
						Y.SetElement(n, m, Step(U.GetElement(n, m)));
					}
				}
			}
			if constexpr(T == ActvFunc::IDENTITY){
				// 恒等活性化関数の場合
				// 行列の要素ごとに活性化関数を掛ける
				for(size_t n = 1; n <= N; ++n){
					for(size_t m = 1; m <= M; ++m){
						Y.SetElement(n, m, Identity(U.GetElement(n, m)));
					}
				}
			}
			if constexpr(T == ActvFunc::SIGMOID){
				// シグモイド活性化関数の場合
				// 行列の要素ごとに活性化関数を掛ける
				for(size_t n = 1; n <= N; ++n){
					for(size_t m = 1; m <= M; ++m){
						Y.SetElement(n, m, Sigmoid(U.GetElement(n, m)));
					}
				}
			}
			if constexpr(T == ActvFunc::TANH){
				// tanh活性化関数の場合
				// 行列の要素ごとにtanh関数を掛ける
				Y = tanhe(U);
			}
			if constexpr(T == ActvFunc::ReLU){
				// ReLU活性化関数の場合
				// 行列の要素ごとに活性化関数を掛ける
				for(size_t n = 1; n <= N; ++n){
					for(size_t m = 1; m <= M; ++m){
						Y.SetElement(n, m, ReLU(U.GetElement(n, m)));
					}
				}
			}
			if constexpr(T == ActvFunc::SOFTMAX){
				// ソフトマックス関数の場合
				// 行列の縦ベクトルごとに活性化関数を掛ける
				for(size_t n = 1; n <= N; ++n){
					setcolumn(Y, Softmax(getcolumn(U, n)), n);
				}
			}
		}

		//! @brief 活性化関数の微分
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		//! @tparam T	活性化関数の種類
		//! @tparam	N	行列の幅
		//! @tparam M	行列の高さ
		template <ActvFunc T, size_t N, size_t M>
		static void fp(const Matrix<N,M>& U, Matrix<N,M>& Y){
			// 活性化関数の種類の設定に従って呼ぶ関数をコンパイル時に変える
			if constexpr(T == ActvFunc::STEP){
				// ステップ活性化関数の場合
				for(size_t n = 1; n <= N; ++n){
					for(size_t m = 1; m <= M; ++m){
						Y.SetElement(n, m, 0);
					}
				}
			}
			if constexpr(T == ActvFunc::IDENTITY){
				// 恒等活性化関数の場合
				for(size_t n = 1; n <= N; ++n){
					for(size_t m = 1; m <= M; ++m){
						Y.SetElement(n, m, 1);
					}
				}
			}
			if constexpr(T == ActvFunc::SIGMOID){
				// シグモイド活性化関数の場合
				for(size_t n = 1; n <= N; ++n){
					for(size_t m = 1; m <= M; ++m){
						Y.SetElement(n, m, DerivativeSigmoid(U.GetElement(n, m)));
					}
				}
			}
			if constexpr(T == ActvFunc::TANH){
				// tanh活性化関数の場合
				Y = Matrix<N,M>::ones() - (U & U);
			}
			if constexpr(T == ActvFunc::ReLU){
				// ReLU活性化関数の場合
				for(size_t n = 1; n <= N; ++n){
					for(size_t m = 1; m <= M; ++m){
						Y.SetElement(n, m, DerivativeReLU(U.GetElement(n, m)));
					}
				}
			}
			if constexpr(T == ActvFunc::SOFTMAX){
				// ソフトマックス関数は内部層では非対応
				arcs_assert(false && "Softmax is not supported in internal layers.");
			}
		}
		
	private:
		ActivationFunctions() = delete;							//!< コンストラクタ使用禁止
		ActivationFunctions(ActivationFunctions&& r) = delete;	//!< ムーブコンストラクタ使用禁止
		~ActivationFunctions() = delete;						//!< デストラクタ使用禁止
		ActivationFunctions(const ActivationFunctions&) = delete;					//!< コピーコンストラクタ使用禁止
		const ActivationFunctions& operator=(const ActivationFunctions&) = delete;	//!< 代入演算子使用禁止
};
}

#endif


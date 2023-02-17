//! @file SimplePerceptron.hh
//! @brief 単純パーセプトロンクラス
//!
//! 1層の単純なパーセプトロン
//!
//! @date 2020/03/05
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2020 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef SIMPLEPERCEPTRON
#define SIMPLEPERCEPTRON

#include <cassert>
#include <array>
#include "Matrix.hh"
#include "ActivationFunctions.hh"

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
	//! @brief 単純パーセプトロンクラス
	//! @tparam N 入力の数, M 学習データ数
	template <unsigned int N, unsigned int M>
	class SimplePerceptron {
		public:
			//! @brief コンストラクタ
			SimplePerceptron()
				: w(), eta(0.5)
			{
				
			}
			
			//! @brief ムーブコンストラクタ
			//! @param[in]	r	右辺値
			SimplePerceptron(SimplePerceptron&& r)
				: w(r.w), eta(r.eta)
			{
				
			}
			
			//! @brief デストラクタ
			~SimplePerceptron(){
				
			}
			
			//! @brief 学習をする関数
			//! @param[in]	x	入力データベクトルの配列
			//! @param[in]	t	教師データベクトル
			//! @param[in]	Gain	学習ゲイン
			//! @param[in]	Epoch	エポック数
			void Train(const std::array<Matrix<1,N+1>, M>& x, const Matrix<1,M>& t, double Gain, unsigned int Epoch){
			 	double Error = 0;	// 誤差
				eta = Gain;			// 学習ゲインをセット
				// エポック数だけ回す
				for(unsigned int i = 0; i < Epoch; ++i){
					printf("Epoch: %u, ", i);
					// データ数だけ回す
					for(unsigned int j = 0; j < M; ++j){
						Error = WeightCalculation(x[j], t.GetElement(1,j+1));	// 重み修正計算
						printf("%f, ", Error);
					}
					printf("\n");
				}
			}
			
			//! @brief 前進計算
			//! @param[in]	x	入力データ
			double ForwardCalculation(const Matrix<1,N+1>& x){
				Matrix<1,1> u = tp(w)*x;	// 重み乗算加算
				return ActivationFunctions::Step(u.GetElement(1,1));	// 活性化関数
			}
			
		private:
			SimplePerceptron(const SimplePerceptron&) = delete;					//!< コピーコンストラクタ使用禁止
			const SimplePerceptron& operator=(const SimplePerceptron&) = delete;//!< 代入演算子使用禁止
			Matrix<1,N+1> w;	//!< 重みベクトル
			double eta;			//!< 学習ゲイン
			
			//! @brief 重み計算
			//! @param[in]	x	入力データ
			//! @param[in]	t	教師データ
			//! @return	誤差
			double WeightCalculation(const Matrix<1,N+1>& x, double t){
				double z = ForwardCalculation(x);	// 前進計算
				double K = (t - z)*eta;				// 教師データとの差分に学習ゲインを乗算
				w += K*x;	// 重み修正計算
				return K;	// 誤差を返す
			}
	};
}

#endif


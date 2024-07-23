//! @file RandomGenerator.hh
//! @brief 乱数生成器
//!
//! メルセンヌ・ツイスタによる指定範囲の一様乱数とガウシアン(正規分布)乱数を生成をするクラス
//!
//! @date 2024/07/23
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#ifndef RANDOMGENERATOR
#define RANDOMGENERATOR

#include <random>
#include "ArcsMatrix.hh"
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
//! @brief 乱数生成器
//! @tparam	T	乱数生成の型 (デフォルト値 = double型)
template<typename T = double>
class RandomGenerator {
	public:
		//! @brief コンストラクタ
		//! @param[in]	MinOrMean	乱数の最小値, ガウシアンの場合は平均値
		//! @param[in]	MaxOrStdDev	乱数の最大値, ガウシアンの場合は標準偏差
		RandomGenerator(const double MinOrMean, const double MaxOrStdDev)
			: RandomDevice(), MersenneTwister(RandomDevice()),
				RandomInt(MinOrMean, MaxOrStdDev),
				RandomDouble(MinOrMean, MaxOrStdDev),
				GaussianRandom(MinOrMean, MaxOrStdDev)
		{
			PassedLog();
			static_assert(std::is_same_v<T,int> || std::is_same_v<T,double>, "RandomGen: Type Error");	// int型か、double型のみ対応
		}

		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		RandomGenerator(RandomGenerator&& r)
			: RandomDevice(), MersenneTwister(r.MersenneTwister),
				RandomInt(r.RandomInt), RandomDouble(r.RandomDouble), GaussianRandom(r.GaussianRandom)
		{

		}

		//! @brief デストラクタ
		~RandomGenerator(){
			PassedLog();
		}
		
		//! @brief 一様乱数を返す関数
		//! @return	一様乱数
		T GetRandom(void){
			// 型によって呼び出す関数を変える
			if constexpr(std::is_same_v<T,int>){
				// int型の場合
				return RandomInt(MersenneTwister);
			}else if constexpr(std::is_same_v<T,double>){
				// double型の場合
				return RandomDouble(MersenneTwister);
			}else{
				arcs_assert(false);	// ここには来ない
			}
		}

		//! @brief 正規分布(ガウシアン)乱数を返す関数
		//! @return	一様乱数
		T GetGaussRandom(void){
			static_assert(std::is_same_v<T,double>, "RandomGen: Type Error");	// double型のみ対応
			return GaussianRandom(MersenneTwister);
		}
		
		//! @brief 乱数シードのリセット
		void ResetSeed(void){
			RandomDevice.entropy();
			MersenneTwister.seed(RandomDevice());
		}

		//! @brief 乱数行列を生成する関数
		//! @param[out]	Y	乱数で埋められた行列
		template <size_t M, size_t N>
		void GetRandomMatrix(ArcsMat<M,N,T>& Y){
			for(size_t n = 1; n <= N; ++n){
				for(size_t m = 1; m <= M; ++m){
					Y(m,n) = GetRandom();	// 行列の各要素に乱数値を書き込む
				}
			}
		}

		//! @brief ガウシアン乱数行列を生成する関数
		//! @param[out]	Y	ガウシアン乱数で埋められた行列
		template <size_t M, size_t N>
		void GetGaussRandomMatrix(ArcsMat<M,N,T>& Y){
			for(size_t n = 1; n <= N; ++n){
				for(size_t m = 1; m <= M; ++m){
					Y(m,n) = GetGaussRandom();	// 行列の各要素に乱数値を書き込む
				}
			}
		}

		//! @brief 乱数行列を生成する関数
		//! @param[out]	Y	乱数行列
		template <size_t N, size_t M>
		void GetRandomMatrix(Matrix<N,M>& Y){
			for(size_t n = 1; n <= N; ++n){
				for(size_t m = 1; m <= M; ++m){
					Y.SetElement(n, m, GetRandom());	// 行列の各要素に乱数値を書き込む
				}
			}
		}

		//! @brief ガウシアン乱数行列を生成する関数
		//! @param[out]	Y	乱数行列
		template <size_t N, size_t M>
		void GetGaussianRandomMatrix(Matrix<N,M>& Y){
			for(size_t n = 1; n <= N; ++n){
				for(size_t m = 1; m <= M; ++m){
					Y.SetElement(n, m, GetGaussRandom());// 行列の各要素に乱数値を書き込む
				}
			}
		}
		
	private:
		RandomGenerator(const RandomGenerator&) = delete;					//!< コピーコンストラクタ使用禁止
		const RandomGenerator& operator=(const RandomGenerator&) = delete;	//!< 代入演算子使用禁止
		std::random_device RandomDevice;				//!< 非決定的乱数生成器
		std::mt19937 MersenneTwister;					//!< メルセンヌ・ツイスタ(32bit版)
		std::uniform_int_distribution<> RandomInt;		//!< 整数用一様乱数
		std::uniform_real_distribution<> RandomDouble;	//!< 浮動小数点用一様乱数
		std::normal_distribution<> GaussianRandom;		//!< 正規分布(ガウシアン)乱数
};
}

#endif


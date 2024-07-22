//! @file RandomGenerator.hh
//! @brief 乱数生成器
//!
//! メルセンヌ・ツイスタによる指定範囲の一様乱数とガウシアン(正規分布)乱数を生成をするクラス
//!
//! @date 2024/07/22
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef RANDOMGENERATOR
#define RANDOMGENERATOR

#include <random>
#include "ArcsMatrix.hh"
#include "Matrix.hh"

namespace ARCS {	// ARCS名前空間
//! @brief 乱数生成器
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
			
		}
		
		//! @brief 一様乱数を返す関数(整数版)
		//! @return	整数乱数
		int GetIntegerRandom(void){
			return RandomInt(MersenneTwister);
		}

		//! @brief 一様乱数を返す関数(浮動小数点版)
		//! @return 浮動小数点乱数
		double GetDoubleRandom(void){
			return RandomDouble(MersenneTwister);
		}
		
		//! @brief 正規分布(ガウシアン)乱数を返す関数
		//! @return 浮動小数点乱数
		double GetGaussianRandom(void){
			return GaussianRandom(MersenneTwister);
		}

		//! @brief 乱数シードのリセット
		void ResetSeed(void){
			RandomDevice.entropy();
			MersenneTwister.seed(RandomDevice());
		}

		//! @brief 乱数行列を生成する関数
		//! @param[out]	Y	乱数行列
		template <size_t M, size_t N, typename T = double>
		void GetRandomMatrix(ArcsMat<M,N,T>& Y){
			for(size_t n = 1; n <= N; ++n){
				for(size_t m = 1; m <= M; ++m){
					Y(m,n) = GetDoubleRandom();	// 行列の各要素に乱数値を書き込む
				}
			}
		}

		//! @brief 乱数行列を生成する関数
		//! @param[out]	Y	乱数行列
		template <size_t N, size_t M>
		void GetRandomMatrix(Matrix<N,M>& Y){
			for(size_t n = 1; n <= N; ++n){
				for(size_t m = 1; m <= M; ++m){
					Y.SetElement(n, m, GetDoubleRandom());	// 行列の各要素に乱数値を書き込む
				}
			}
		}

		//! @brief ガウシアン乱数行列を生成する関数
		//! @param[out]	Y	乱数行列
		template <size_t N, size_t M>
		void GetGaussianRandomMatrix(Matrix<N,M>& Y){
			for(size_t n = 1; n <= N; ++n){
				for(size_t m = 1; m <= M; ++m){
					Y.SetElement(n, m, GetGaussianRandom());// 行列の各要素に乱数値を書き込む
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


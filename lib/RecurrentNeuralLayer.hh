//! @file RecurrentNeuralLayer.hh
//! @brief 再帰ニューラルレイヤクラス
//!
//! 再帰ニューラルネット1層分のレイヤ
//!
//! @date 2021/08/19
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2021 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef RECURRENTNEURALLAYER
#define RECURRENTNEURALLAYER

#include <cassert>
#include <array>
#include <cmath>
#include <string>
#include <array>
#include "Matrix.hh"
#include "NeuralNetParamDef.hh"
#include "ActivationFunctions.hh"
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
//! @brief 再帰ニューラルレイヤクラス
//! @tparam N	入力チャネル数
//! @tparam P	パーセプトロン数
//! @tparam T	時系列方向のデータ数
//! @tparam W	入力ウィンドウ幅
//! @tparam M	ミニバッチ数
//! @tparam	AF	活性化関数の種類
//! @tparam	IT	初期化のタイプ
//! @tparam	GD	勾配降下法のタイプ
//! @tparam	DD	ドロップアウト動作フラグ
template <
	size_t N,
	size_t P,
	size_t T,
	size_t W,
	size_t M,
	ActvFunc AF,
	NnInitTypes IT = NnInitTypes::XAVIER,
	NnDescentTypes GD = NnDescentTypes::MOMENTUM,
	NnDropout DD = NnDropout::DISABLE
>
class RecurrentNeuralLayer {
	public:
		// レイヤーの入出力変数
		std::array<Matrix<1,P>, W + 2> z;	//!< 活性化関数通過後の状態ベクトルの時系列配列(範囲 t = 1 … T, t = 0 と T + 1 の分も確保)
		std::array<Matrix<1,N>, W + 2> dLdz;//!< 勾配ベクトルの時系列配列(範囲 t = 1 … T, t = 0 と T + 1 の分も確保)
		Matrix<1,N> dLde;					//!< 勾配ベクトルの時系列配列(出力層用)
		
		//! @brief コンストラクタ
		RecurrentNeuralLayer()
		 : z(), dLdz(), dLde(), u(), d(), y(), e(), Wl(), Wt(), b(), dWl(), dWt(), db(), fpu(), 
		   DropRand(0, 1), DropMask(),
		   epsilon(0.001), alpha(0.00001)
		{
			
		}
		
		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		RecurrentNeuralLayer(RecurrentNeuralLayer&& r)
		{
			
		}
		
		//! @brief デストラクタ
		~RecurrentNeuralLayer(){
			
		}
		
		//! @brief 順伝播計算(ベクトル入出力版)
		//! @param[in]	zprev	前の層からの入力ベクトル
		//! @param[in]	t		時刻インデックス
		void PropagateForward(const std::array<Matrix<1,N>, W + 2>& zprev, const size_t t){
			u.at(t) = Wl*zprev.at(t) + Wt*z.at(t - 1) + b;		// 重み乗算加算とバイアス加算
			ActivationFunctions::f<AF,1,P>(u.at(t), z.at(t));	// 活性化関数
		}
		
		//! @brief 逆伝播計算(入力層と内部層用，ベクトル入出力版)
		//! @param[in]	dLdznext	後ろの層からの勾配ベクトル
		//! @param[in]	t			時刻インデックス
		void PropagateBackward(const std::array<Matrix<1,P>, W + 2>& dLdznext, const size_t t){
			ActivationFunctions::fp<AF,1,P>(u.at(t), fpu);			// 活性化関数の微分を通す計算
			d.at(t) = fpu & (dLdznext.at(t) + tp(Wt)*d.at(t + 1));	// 誤差ベクトルの計算 ←ここのdLdznextのtの関数はどうするの？
			
			dLdz.at(t) = tp(Wl)*d.at(t);		// 前の層に渡すための勾配計算
		}
		
		//! @brief 逆伝播計算(-)
		//! @param[in]	dLdznext	後ろの層からの勾配ベクトル
		//! @param[in]	t			時刻インデックス
		void PropagateBackward(const Matrix<1,P>& dLdenext, const size_t t){
			ActivationFunctions::fp<AF,1,P>(z.at(t), fpu);		// 活性化関数の微分を通す計算
			d.at(t) = fpu & (dLdenext + tp(Wt)*d.at(t + 1));	// 誤差ベクトルの計算 ←ここのdLdznextのtの関数はどうするの？
			
			dLdz.at(t) = tp(Wl)*d.at(t);		// 前の層に渡すための勾配計算
		}
		
		//! @brief 順伝播計算(出力層用，ベクトル入出力版)
		void PropagateForwardForOutput(const std::array<Matrix<1,N>, W + 2>& zprev){
			ActivationFunctions::f<AF,1,P>(Wl*zprev.at(W) + b, y);	// 重み乗算加算とバイアス加算と活性化関数
		}
		
		//! @brief 出力の重み誤差ベクトル計算(出力層用，ベクトル入出力訓練版)
		//! @param[in]	r	目標値ベクトル
		//! @param[in]	t	時刻インデックス
		void PropagateBackwardForOutput(const Matrix<1,P>& r){
			e = y - r;				// 誤差ベクトルの計算
			dLde = tp(Wl)*e;		// 前の層に渡すための勾配計算
		}
		
		//! @brief 重み行列とバイアスベクトルの更新
		//! @param[in]	zprev	前の層からの入力ベクトル
		void UpdateWeightAndBias(const std::array<Matrix<1,N>, W + 2>& zprev){
			// 重み・バイアス勾配の初期化
			dWl.FillAllZero();
			dWt.FillAllZero();
			db.FillAllZero();
			
			// 重み・バイアス勾配の計算
			for(size_t t = 1; t <= W; ++t){
				dWl += d.at(t)*tp(zprev.at(t));	// 階層方向の重み勾配
				dWt += d.at(t)*tp(z.at(t - 1));	// 時刻方向の重み勾配
				db  += d.at(t);					// バイアス勾配
			}
			
			GetUpdatedValue(dWl, dWt, db);	// 確率的勾配降下法
		}
		
		//! @brief 重み行列とバイアスベクトルの更新(出力層用)
		//! @param[in]	zprev	前の層からの入力ベクトル
		void UpdateWeightAndBiasForOutput(const std::array<Matrix<1,N>, W + 2>& zprev){
			// 重み・バイアス勾配の初期化
			dWl.FillAllZero();
			db.FillAllZero();
			
			// 重み・バイアス勾配の計算
			for(size_t t = 1; t <= W; ++t){
				dWl += e*tp(zprev.at(t));	// 階層方向の重み勾配
				db  += e;					// バイアス勾配
			}
			
			GetUpdatedValue(dWl, dWt, db);	// 確率的勾配降下法
		}
		
		void ClearStateVars(void){
			for(size_t t = 0; t < W + 2; ++t){
				z.at(t).FillAllZero();
				dLdz.at(t).FillAllZero();
				u.at(t).FillAllZero();
				d.at(t).FillAllZero();
			}
		}
		
		//! @brief 重み行列の初期化
		//! @param[in]	Nprev	前層のユニット数
		void InitWeight(const size_t Nprev){
			if constexpr(IT == NnInitTypes::XAVIER){
				// Xavierの初期化
				InitWeightByRandom( 1.0/sqrt((double)Nprev) );
			}
			if constexpr(IT == NnInitTypes::HE){
				// Heの初期化
				InitWeightByRandom( sqrt(2.0/(double)Nprev) );
			}
		}
		
		//! @brief 重み行列とバイアスベクトルの表示
		void DispWeightAndBias(void){
			PrintMat(Wl);
			PrintMat(Wt);
			PrintMat(b);
		}
		
		void DispError(void){
			PrintMatrix(e, "% 16.8f");
		}
		
		//! @brief ドロップアウトマスクの生成
		void GenerateDropMask(void){
			if constexpr(DD == NnDropout::ENABLE){
				// ドロップアウトするときのみ下記を計算
				DropRand.GetRandomMatrix(DropMask);	// 乱数生成
				
				// ドロップアウト率よりも大きければマスクを0にする
				for(size_t i = 0; i < P; ++i){
					if(DropRate < DropMask.GetElement(1,i+1)){
						DropMask.SetElement(1,i+1, 0);
					}else{
						DropMask.SetElement(1,i+1, 1);
					}
				}
			}
		}
		
	private:
		static constexpr Matrix<1,M> l = Matrix<1,M>::ones();	//!< 1ベクトル
		std::array<Matrix<1,P>, W + 2> u;	//!< 状態ベクトルの時系列配列(範囲 t = 1 … T, t = 0 と T + 1 の分も確保)
		std::array<Matrix<1,P>, W + 2> d;	//!< 誤差ベクトルの時系列配列(範囲 t = 1 … T, t = 0 と T + 1 の分も確保)
		Matrix<1,P> y, e;
		Matrix<N,P> Wl;			//!< 階層方向の重み行列
		Matrix<P,P> Wt;			//!< 時刻方向の重み行列
		Matrix<1,P> b;			//!< バイアスベクトル
		Matrix<N,P> dWl;		//!< 階層方向の重み更新差分値
		Matrix<P,P> dWt;		//!< 時刻方向の重み更新差分値
		Matrix<1,P> db;			//!< バイアス更新差分値
		Matrix<1,P> fpu;		//!< 活性化関数の微分通過後のベクトル
		RandomGenerator<double> DropRand;	//!< ドロップアウト用メルセンヌ・ツイスタ
		Matrix<1,P> DropMask;		//!< ドロップアウト用マスクベクトル
		
		static constexpr double DropRate = 0.5;	//!< ドロップアウト率
		double epsilon;			//!< 更新ゲインε (SGD, Momentum, AdaGrad, RMSprop)
		double alpha;
		
		//! @brief 重み行列の正規分布乱数による初期化
		//! @param[in]	sigma	正規分布乱数の標準偏差
		void InitWeightByRandom(const double sigma){
			RandomGenerator RandWl(0, sigma);	// メルセンヌ・ツイスタの生成
			RandomGenerator RandWt(0, sigma);	// メルセンヌ・ツイスタの生成
			RandWl.GetGaussianRandomMatrix(Wl);	// 平均0，標準偏差σのガウシアン乱数行列の取得
			//RandWt.GetGaussianRandomMatrix(Wt);	// 平均0，標準偏差σのガウシアン乱数行列の取得
			//Wl = Wl * 0.01;
			//Wt = Wt * 0.01;
		}
		
		//! @brief ヴァニラ確率的勾配降下法
		//! @param[in]	DiffWl	重み更新差分値
		//! @param[in]	DiffWt	重み更新差分値
		//! @param[in]	diffb	バイアス更新差分値
		void GetUpdatedValue(const Matrix<N,P>& DiffWl, const Matrix<P,P>& DiffWt, const Matrix<1,P>& diffb){
			
			Wl += (-epsilon*DiffWl);	// 重み行列の更新
			Wt += (-epsilon*DiffWt);	// 重み行列の更新
			b  += (-epsilon*diffb);		// バイアスベクトルの更新
			
			//dWl = alpha*dWl - epsilon*DiffWl;	// 更新ゲイン乗算後の重み更新差分値
			//dWt = alpha*dWt - epsilon*DiffWt;	// 更新ゲイン乗算後の重み更新差分値
			//db = alpha*db - epsilon*diffb;	// 更新ゲイン乗算後のバイアス更新差分値
			//Wl += dWl;					// 重み行列の更新
			//Wt += dWt;					// 重み行列の更新
			//b += db;					// バイアスベクトルの更新

		}

};
}

#endif


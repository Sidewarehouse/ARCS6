//! @file RecurrentNeuralNet3.hh
//! @brief 再帰ニューラルネット 3層版
//!
//! 再帰ニューラルネット 3層版
//!
//! @date 2021/08/19
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2021 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef RECURRENTNEURALNET3
#define RECURRENTNEURALNET3

#include <cassert>
#include "RecurrentNeuralLayer.hh"
#include "TimeSeriesDatasets.hh"

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
//! @brief クラステンプレート
//! @tparam 
template <
	size_t Nin,		// ユニット数(入力層)
	size_t NL1,		// ユニット数(内部層1)
	size_t Nout,	// ユニット数(出力層)
	size_t Tlen,	// 時刻データの長さ
	size_t Wind,	// 入力ウィンドウ幅
	size_t Mbat,	// ミニバッチサイズ
	size_t Epch,	// エポック数
	ActvFunc FuncIn = ActvFunc::TANH,		// 入力層の活性化関数
	ActvFunc FuncL1 = ActvFunc::TANH,		// 内部層1の活性化関数
	ActvFunc FuncOut = ActvFunc::TANH,		// 出力層の活性化関数
	NnInitTypes InitType = NnInitTypes::HE,	// 重み初期化のタイプ
	NnDescentTypes GradDesType = NnDescentTypes::SGD,	// 勾配降下法のタイプ
	NnDropout DropOutEnable	= NnDropout::DISABLE		// ドロップアウトイネーブル
>
class RecurrentNeuralNet3 {
	public:
		//! @brief コンストラクタ
		RecurrentNeuralNet3(TimeSeriesDatasets<Nin,Nout,Tlen,Wind,Mbat>& TSD)
			: Dataset(TSD), RNNin(), RNNL1(), RNNout()
		{
			PassedLog();
		}

		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		RecurrentNeuralNet3(RecurrentNeuralNet3&& r)
			: Dataset(r.Dataset)
		{
			
		}

		//! @brief デストラクタ
		~RecurrentNeuralNet3(){
			PassedLog();
		}
		
		//! @brief 誤差逆伝播法による訓練をする関数
		void Train(void){
			// 重み行列の初期化
			RNNin.InitWeight(Nin);		// 入力層の重み行列の乱数による初期化
			//RNNL1.InitWeight(Nin);		// 内部層1の重み行列の乱数による初期化
			RNNout.InitWeight(30);		// 出力層の重み行列の乱数による初期化
			
			//RNNin.DispWeightAndBias();
			//RNNL1.DispWeightAndBias();
			//RNNout.DispWeightAndBias();
			
			// エポック数分のループ
			for(size_t i = 1; i <= Epch; ++i){
				//RNNin.GenerateDropMask();
				//RNNL1.GenerateDropMask();
				//RNNout.GenerateDropMask();
				
				// 順伝播計算
				for(size_t t = 1; t <= Wind; ++t){
					RNNin.PropagateForward(Dataset.InputData, t);		// 入力層の順伝播計算
				}
				RNNout.PropagateForwardForOutput(RNNin.z);				// 出力層の順伝播計算
				
				// 逆伝播計算
				RNNout.PropagateBackwardForOutput(Dataset.TrainData);	// 出力層の逆伝播計算
				for(size_t t = Wind; 1 <= t; --t){
					RNNin.PropagateBackward(RNNout.dLde, t);			// 出力層の逆伝播計算
				}
				
				// 重み・バイアス更新計算
				RNNin.UpdateWeightAndBias(Dataset.InputData);
				RNNout.UpdateWeightAndBiasForOutput(RNNin.z);
				
				if(i % 10 == 1) RNNout.DispError();
				
				RNNin.ClearStateVars();
				RNNout.ClearStateVars();
				
			}
			
			//RNNin.DispWeightAndBias();
			//RNNL1.DispWeightAndBias();
			//RNNout.DispWeightAndBias();
		}
		
	private:
		RecurrentNeuralNet3(const RecurrentNeuralNet3&) = delete;					//!< コピーコンストラクタ使用禁止
		const RecurrentNeuralNet3& operator=(const RecurrentNeuralNet3&) = delete;	//!< 代入演算子使用禁止
		
		static constexpr size_t EpchDisp = 10;					//!< エポックループ表示数
		TimeSeriesDatasets<Nin,Nout,Tlen,Wind,Mbat>& Dataset;	//!< 時系列データセットへの参照
		
		RecurrentNeuralLayer<Nin,30, Tlen,Wind,Mbat,FuncIn, InitType,GradDesType,DropOutEnable> RNNin;	//!< 入力層
		RecurrentNeuralLayer<Nin,NL1, Tlen,Wind,Mbat,FuncL1, InitType,GradDesType,DropOutEnable> RNNL1;	//!< 内部層1
		RecurrentNeuralLayer<30,Nout,Tlen,Wind,Mbat,FuncOut,InitType,GradDesType,DropOutEnable> RNNout;	//!< 出力層
};
}

#endif


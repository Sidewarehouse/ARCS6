//! @file IrisDatasets.hh
//! @brief 機械学習試験用アヤメデータセットクラス(分類番号版)
//!
//! 機械学習の試験用ための，アヤメの「ガク，花びらの計測値」と「アヤメの種類の分類」のデータセットを集めたクラス
//!
//! @date 2020/05/17
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef IRISDATASETS
#define IRISDATASETS

#include <cassert>
#include <array>
#include "Matrix.hh"
#include "Shuffle.hh"
#include "BatchProcessor.hh"

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
//! @brief 機械学習試験用アヤメデータセットクラス(分類番号版)
//! @tparam	M	ミニバッチ数
template <size_t M>
class IrisDatasets {
	public:
		// データセット情報の定義
		static constexpr size_t D = 150;//!< データセット数
		static constexpr size_t N = 4;	//!< 計測データチャネル数
		static constexpr size_t C = 1;	//!< 訓練データチャネル数
		static constexpr size_t FinalMinbatNum = D/M - 1;	//!< 最後のミニバッチ番号

		// アヤメの分類の定義
		static constexpr double IRIS_SETOSA     = 1;	//!< ヒオウギアヤメ
		static constexpr double IRIS_VERSICOLOR = 2;	//!< ブルーフラッグ
		static constexpr double IRIS_VIRGINICA  = 3;	//!< バージニカ
		Matrix<N,D> MeasuredData;	//!< 前処理済みのアヤメの計測データ
		Matrix<C,D> ClassData;		//!< 前処理済みのアヤメの分類データ
		
		//! @brief コンストラクタ
		IrisDatasets()
			: MeasuredData(), ClassData(), Shfl()
		{
			PassedLog();
			MeasuredData.LoadArray(IrisMeasurement);	// 計測データを行列として読み込み
			ClassData.LoadArray(IrisClassification);	// 分類データを行列として読み込み
			ShuffleDatasets();							// シャッフルの実行
		}

		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		IrisDatasets(IrisDatasets&& r)
			: MeasuredData(r.MeasuredData), ClassData(r.ClassData), Shfl()
		{
			
		}
		
		//! @brief デストラクタ
		~IrisDatasets(){
			PassedLog();
		}
		
		//! @brief 計測データを表示する関数
		void DispMeasuredData(void){
			printf("\nIris Mesurement Data:\n");
			PrintMatrix(MeasuredData, "% 3.1f");
		}
		
		//! @brief 分類データを表示する関数
		void DispClassData(void){
			printf("\nIris Classification Data:\n");
			PrintMatrix(ClassData, "% 1.0f");
		}
		
		//! @brief i番目のミニバッチ計測データを取得する関数
		//! @param[in]	i	ミニバッチ番号
		//! @param[out]	Y	計測データ
		void GetMeasuredBatchData(const size_t i, Matrix<M,N>& Y){
			arcs_assert(i <= FinalMinbatNum);	// ミニバッチ番号の範囲チェック
			Y = BatchProcessor::GetMiniBatchData<N,D,M>(MeasuredData, i);	// データセットを切り出してミニバッチにする
		}
		
		//! @brief i番目のミニバッチ分類データを取得する関数
		//! @param[out]	y	訓練データ
		//! @param[in]	i	ミニバッチ番号
		void GetClassBatchData(const size_t i, Matrix<M,C>& y){
			arcs_assert(i <= FinalMinbatNum);	// ミニバッチ番号の範囲チェック
			y = BatchProcessor::GetMiniBatchData<C,D,M>(ClassData, i);		// データセットを切り出してミニバッチにする
		}
		
		//! @brief データセットをシャッフルする関数
		void ShuffleDatasets(void){
			Shfl.ShuffleMatrixRow(MeasuredData, ClassData);	// シャッフルの実行
		}
		
	private:
		IrisDatasets(const IrisDatasets&) = delete;					//!< コピーコンストラクタ使用禁止
		const IrisDatasets& operator=(const IrisDatasets&) = delete;//!< 代入演算子使用禁止
		Shuffle Shfl;	//!< メルセンヌ・ツイスタ・シャッフル
					
		//! @brief アヤメの計測データ
		//! {ガクの長さ [cm], ガクの幅 [cm], 花びらの長さ [cm], 花びらの幅 [cm]}
		static constexpr std::array<std::array<double, N>, D> IrisMeasurement = {{
			{5.1,3.5,1.4,0.2},
			{4.9,3.0,1.4,0.2},
			{4.7,3.2,1.3,0.2},
			{4.6,3.1,1.5,0.2},
			{5.0,3.6,1.4,0.2},
			{5.4,3.9,1.7,0.4},
			{4.6,3.4,1.4,0.3},
			{5.0,3.4,1.5,0.2},
			{4.4,2.9,1.4,0.2},
			{4.9,3.1,1.5,0.1},
			{5.4,3.7,1.5,0.2},
			{4.8,3.4,1.6,0.2},
			{4.8,3.0,1.4,0.1},
			{4.3,3.0,1.1,0.1},
			{5.8,4.0,1.2,0.2},
			{5.7,4.4,1.5,0.4},
			{5.4,3.9,1.3,0.4},
			{5.1,3.5,1.4,0.3},
			{5.7,3.8,1.7,0.3},
			{5.1,3.8,1.5,0.3},
			{5.4,3.4,1.7,0.2},
			{5.1,3.7,1.5,0.4},
			{4.6,3.6,1.0,0.2},
			{5.1,3.3,1.7,0.5},
			{4.8,3.4,1.9,0.2},
			{5.0,3.0,1.6,0.2},
			{5.0,3.4,1.6,0.4},
			{5.2,3.5,1.5,0.2},
			{5.2,3.4,1.4,0.2},
			{4.7,3.2,1.6,0.2},
			{4.8,3.1,1.6,0.2},
			{5.4,3.4,1.5,0.4},
			{5.2,4.1,1.5,0.1},
			{5.5,4.2,1.4,0.2},
			{4.9,3.1,1.5,0.1},
			{5.0,3.2,1.2,0.2},
			{5.5,3.5,1.3,0.2},
			{4.9,3.1,1.5,0.1},
			{4.4,3.0,1.3,0.2},
			{5.1,3.4,1.5,0.2},
			{5.0,3.5,1.3,0.3},
			{4.5,2.3,1.3,0.3},
			{4.4,3.2,1.3,0.2},
			{5.0,3.5,1.6,0.6},
			{5.1,3.8,1.9,0.4},
			{4.8,3.0,1.4,0.3},
			{5.1,3.8,1.6,0.2},
			{4.6,3.2,1.4,0.2},
			{5.3,3.7,1.5,0.2},
			{5.0,3.3,1.4,0.2},
			{7.0,3.2,4.7,1.4},
			{6.4,3.2,4.5,1.5},
			{6.9,3.1,4.9,1.5},
			{5.5,2.3,4.0,1.3},
			{6.5,2.8,4.6,1.5},
			{5.7,2.8,4.5,1.3},
			{6.3,3.3,4.7,1.6},
			{4.9,2.4,3.3,1.0},
			{6.6,2.9,4.6,1.3},
			{5.2,2.7,3.9,1.4},
			{5.0,2.0,3.5,1.0},
			{5.9,3.0,4.2,1.5},
			{6.0,2.2,4.0,1.0},
			{6.1,2.9,4.7,1.4},
			{5.6,2.9,3.6,1.3},
			{6.7,3.1,4.4,1.4},
			{5.6,3.0,4.5,1.5},
			{5.8,2.7,4.1,1.0},
			{6.2,2.2,4.5,1.5},
			{5.6,2.5,3.9,1.1},
			{5.9,3.2,4.8,1.8},
			{6.1,2.8,4.0,1.3},
			{6.3,2.5,4.9,1.5},
			{6.1,2.8,4.7,1.2},
			{6.4,2.9,4.3,1.3},
			{6.6,3.0,4.4,1.4},
			{6.8,2.8,4.8,1.4},
			{6.7,3.0,5.0,1.7},
			{6.0,2.9,4.5,1.5},
			{5.7,2.6,3.5,1.0},
			{5.5,2.4,3.8,1.1},
			{5.5,2.4,3.7,1.0},
			{5.8,2.7,3.9,1.2},
			{6.0,2.7,5.1,1.6},
			{5.4,3.0,4.5,1.5},
			{6.0,3.4,4.5,1.6},
			{6.7,3.1,4.7,1.5},
			{6.3,2.3,4.4,1.3},
			{5.6,3.0,4.1,1.3},
			{5.5,2.5,4.0,1.3},
			{5.5,2.6,4.4,1.2},
			{6.1,3.0,4.6,1.4},
			{5.8,2.6,4.0,1.2},
			{5.0,2.3,3.3,1.0},
			{5.6,2.7,4.2,1.3},
			{5.7,3.0,4.2,1.2},
			{5.7,2.9,4.2,1.3},
			{6.2,2.9,4.3,1.3},
			{5.1,2.5,3.0,1.1},
			{5.7,2.8,4.1,1.3},
			{6.3,3.3,6.0,2.5},
			{5.8,2.7,5.1,1.9},
			{7.1,3.0,5.9,2.1},
			{6.3,2.9,5.6,1.8},
			{6.5,3.0,5.8,2.2},
			{7.6,3.0,6.6,2.1},
			{4.9,2.5,4.5,1.7},
			{7.3,2.9,6.3,1.8},
			{6.7,2.5,5.8,1.8},
			{7.2,3.6,6.1,2.5},
			{6.5,3.2,5.1,2.0},
			{6.4,2.7,5.3,1.9},
			{6.8,3.0,5.5,2.1},
			{5.7,2.5,5.0,2.0},
			{5.8,2.8,5.1,2.4},
			{6.4,3.2,5.3,2.3},
			{6.5,3.0,5.5,1.8},
			{7.7,3.8,6.7,2.2},
			{7.7,2.6,6.9,2.3},
			{6.0,2.2,5.0,1.5},
			{6.9,3.2,5.7,2.3},
			{5.6,2.8,4.9,2.0},
			{7.7,2.8,6.7,2.0},
			{6.3,2.7,4.9,1.8},
			{6.7,3.3,5.7,2.1},
			{7.2,3.2,6.0,1.8},
			{6.2,2.8,4.8,1.8},
			{6.1,3.0,4.9,1.8},
			{6.4,2.8,5.6,2.1},
			{7.2,3.0,5.8,1.6},
			{7.4,2.8,6.1,1.9},
			{7.9,3.8,6.4,2.0},
			{6.4,2.8,5.6,2.2},
			{6.3,2.8,5.1,1.5},
			{6.1,2.6,5.6,1.4},
			{7.7,3.0,6.1,2.3},
			{6.3,3.4,5.6,2.4},
			{6.4,3.1,5.5,1.8},
			{6.0,3.0,4.8,1.8},
			{6.9,3.1,5.4,2.1},
			{6.7,3.1,5.6,2.4},
			{6.9,3.1,5.1,2.3},
			{5.8,2.7,5.1,1.9},
			{6.8,3.2,5.9,2.3},
			{6.7,3.3,5.7,2.5},
			{6.7,3.0,5.2,2.3},
			{6.3,2.5,5.0,1.9},
			{6.5,3.0,5.2,2.0},
			{6.2,3.4,5.4,2.3},
			{5.9,3.0,5.1,1.8}
		}};

		//! @brief アヤメの分類データ
		static constexpr std::array<std::array<double, C>, D> IrisClassification = {{
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_SETOSA},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VERSICOLOR},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA},
			{IRIS_VIRGINICA}
		}};
};
}

#endif


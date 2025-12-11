//! @file MovingDifferentiator.hh
//! @brief 移動微分器
//!
//! 移動窓内の先頭の値と最後尾の値，および時刻計測値から微分値を計算する
//!
//! @date 2025/12/11
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2025 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef MOVINGDIFFERENTIATOR
#define MOVINGDIFFERENTIATOR

#include <cmath>
#include "RingBuffer.hh"
#include "ArcsMatrix.hh"

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
//! @brief 移動微分器
//! @tparam	W	差分のサンプリング数 [-]
//! @tparam T	型名(デフォルトはdouble型)
template <size_t W, typename T = double>
class MovingDifferentiator {
	public:
		//! @brief コンストラクタ
		MovingDifferentiator()
			: VarWindow(), TimeWindow(), FirstTime(true)
		{
			
		}

		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		MovingDifferentiator(MovingDifferentiator&& r)
			: VarWindow(r.VarWindow), TimeWindow(r.TimeWindow), FirstTime(r.FirstTime)
		{
			
		}

		//! @brief デストラクタ
		~MovingDifferentiator(){
			
		}
		
		//! @brief 微分値を計算する関数(double版)
		//! @param[in]	Var		[*] 現在値入力
		//! @param[in]	Time	[s] 現在時刻
		//! @return	微分値 [*/s]
		double GetSignal(const double Var, const double Time){
			if(FirstTime){
				// 初めて呼ばれた場合は，
				VarWindow.FillBuffer(Var);		// 初期位置でバッファを埋めておく
				TimeWindow.FillBuffer(Time);	// 初期時刻でバッファを埋めておく
				FirstTime = false;				// フラグリセット
				return 0;						// 最初の1回目は速度を計算できないのでゼロ
			}else{
				// 2回目以降は，
				const T dx = Var - VarWindow.GetFinalValue();		// 現在位置と位置リングバッファの最後尾との偏差の計算
				const double dt = Time - TimeWindow.GetFinalValue();// 現在時刻と時刻リングバッファの最後尾との偏差の計算
				const T v = dx/dt;				// 微分の計算
				VarWindow.SetFirstValue(Var);	// 位置リングバッファの先頭に現在位置を詰める
				TimeWindow.SetFirstValue(Time);	// 時刻リングバッファの先頭に現在時刻を詰める
				if(std::isfinite(v)){
					return v;	// 微分結果が有限値なら計算結果を返す
				}else{
					return 0;	// ゼロ割のときはゼロにしておく
				}
			}
		}
		
		//! @brief 微分値を計算する関数(行列版)
		//! @param[in]	Var		[*] 現在値入力行列
		//! @param[in]	Time	[s] 現在時刻
		//! @return	微分値 [*/s]
		template<size_t NN, size_t MM>
		Matrix<NN,MM> GetSignal(const Matrix<NN,MM>& Var, const double Time){
			if(FirstTime){
				// 初めて呼ばれた場合は，
				VarWindow.FillBuffer(Var);		// 初期位置でバッファを埋めておく
				TimeWindow.FillBuffer(Time);	// 初期時刻でバッファを埋めておく
				FirstTime = false;				// フラグリセット
				return Matrix<NN,MM>::zeros();	// 最初の1回目は速度を計算できないのでゼロ
			}else{
				// 2回目以降は，
				const T dx = Var - VarWindow.GetFinalValue();		// 現在位置と位置リングバッファの最後尾との偏差の計算
				const double dt = Time - TimeWindow.GetFinalValue();// 現在時刻と時刻リングバッファの最後尾との偏差の計算
				const T v = dx/dt;				// 微分の計算
				VarWindow.SetFirstValue(Var);	// 位置リングバッファの先頭に現在位置を詰める
				TimeWindow.SetFirstValue(Time);	// 時刻リングバッファの先頭に現在時刻を詰める
				return v;						// 計算結果を返す
			}
		}

			//! @brief 微分値を計算する関数(ArcsMat行列版)
		//! @param[in]	Var		[*] 現在値入力行列
		//! @param[in]	Time	[s] 現在時刻
		//! @return	微分値 [*/s]
		template<size_t NN, size_t MM>
		ArcsMat<NN,MM> GetSignal(const ArcsMat<NN,MM>& Var, const double Time){
			if(FirstTime){
				// 初めて呼ばれた場合は，
				VarWindow.FillBuffer(Var);		// 初期位置でバッファを埋めておく
				TimeWindow.FillBuffer(Time);	// 初期時刻でバッファを埋めておく
				FirstTime = false;				// フラグリセット
				return ArcsMat<NN,MM>::zeros();	// 最初の1回目は速度を計算できないのでゼロ
			}else{
				// 2回目以降は，
				const T dx = Var - VarWindow.GetFinalValue();		// 現在位置と位置リングバッファの最後尾との偏差の計算
				const double dt = Time - TimeWindow.GetFinalValue();// 現在時刻と時刻リングバッファの最後尾との偏差の計算
				const T v = dx/dt;				// 微分の計算
				VarWindow.SetFirstValue(Var);	// 位置リングバッファの先頭に現在位置を詰める
				TimeWindow.SetFirstValue(Time);	// 時刻リングバッファの先頭に現在時刻を詰める
				return v;						// 計算結果を返す
			}
		}
		
		//! @brief リセット
		void Reset(void){
			FirstTime = true;
		}
		
	private:
		MovingDifferentiator(const MovingDifferentiator&) = delete;					//!< コピーコンストラクタ使用禁止
		const MovingDifferentiator& operator=(const MovingDifferentiator&) = delete;//!< 代入演算子使用禁止
		RingBuffer<T, W, false> VarWindow;			//!< 変数データ用リングバッファ
		RingBuffer<double, W, false> TimeWindow;	//!< 時間データ用リングバッファ
		bool FirstTime;								//!< 初めて呼ばれたかどうかのフラグ
};
}

#endif


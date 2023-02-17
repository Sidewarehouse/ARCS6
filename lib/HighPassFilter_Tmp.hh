//! @file HighPassFilter.hh
//! @brief 高域通過濾波器クラス
//!
//! 1次高域通過濾波器 G(s)=s/(s+g) (双一次変換)
//! テンプレート形式で挑戦 R020610_川合 動作確認済み
//!
//! @date 2019/02/24
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef HIGHPASSFILTER_TMP
#define HIGHPASSFILTER_TMP

#include <tuple>
#include "Matrix.hh"
#include "Discret.hh"
#include "ARCSeventlog.hh"

namespace ARCS {	// ARCS名前空間
	//! @brief 高域通過濾波器クラス
	template <size_t N = 1>
	class HighPassFilter {
		public:
			//! @brief コンストラクタ(スカラー版)
			//! @param [in] CutFreq 遮断周波数 [rad/s]
			//! @param [in] SmplTime 制御周期 [s]
			HighPassFilter(double CutFreq, double SmplTime)
				: Ts(SmplTime),	// [s]		制御周期の格納
				  g(CutFreq),	// [rad/s]	遮断周波数
				  uZ1(0), yZ1(0), u(0), y(0),
				  gv(), uZ1v(), yZ1v(),uv(),yv()
			{
				PassedLog();
			}

			//! @brief コンストラクタ(ベクトル版)
			//! @param [in] CutFreq 遮断周波数 [rad/s]
			//! @param [in] SmplTime 制御周期 [s]
			HighPassFilter(Matrix<1,N> CutFreq, double SmplTime)
				: Ts(SmplTime),	// [s]		制御周期の格納
				  g(),	// [rad/s]	遮断周波数
				  uZ1(), yZ1(), u(), y(),
				  gv(CutFreq), uZ1v(0), yZ1v(0), uv(0),yv(0)
			{
				PassedLog();
			}

			//! @brief HPF後の出力を取得する関数(スカラー版)
			//! @param [in] u 入力信号
			double GetSignal(double u){
				double y;
				y=1.0/(g*Ts+2.0)*( 2.0*(u-uZ1) - yZ1*(g*Ts-2.0) );

				uZ1=u;
				yZ1=y;

				return y;
			}

			//! @brief HPF後の出力を取得する関数(ベクトル版)
			//! @paramre [in] uv 入力ベクトル
			Matrix<1,N> GetSignal(Matrix<1,N> uv){
				Matrix<1,N> yv;
				// ベクトルの数だけ回す, Matlabの配列仕様と合わせたため、ARCSのMatrix関数はi=1から格納される。
				// i=0は-1 or NA になっているはず。そういう仕様。
				for(size_t i = 1; i <= N; ++i){
					yv[i]=1.0/(gv[i]*Ts+2.0)*( 2.0*(uv[i]-uZ1v[i]) - yZ1v[i]*(gv[i]*Ts-2.0) );

					uZ1v[i]=uv[i];
					yZ1v[i]=yv[i];
				}
				return yv;
			}

			//! @brief ムーブコンストラクタ(スカラー版のみ)
			HighPassFilter(HighPassFilter&& right)
				: Ts(right.Ts),	// [s]		制御周期の格納
				  g(right.g),	// [rad/s]	遮断周波数
				  uZ1(right.uZ1), yZ1(right.yZ1),
				  gv(), uZ1v(), yZ1v(),uv(),yv()
			{
				PassedLog();
			}
			//! @brief デストラクタ
			~HighPassFilter(){
				PassedLog();
			}

		private:
			HighPassFilter(const HighPassFilter&) = delete;					//!< コピーコンストラクタ使用禁止
			const HighPassFilter& operator=(const HighPassFilter&) = delete;//!< 代入演算子使用禁止

			// 共通変数
			double Ts;	//!< [s] 制御周期

			// スカラー版の変数
			double g;	//!< [rad/s]  遮断周波数
			double u=0;	//!<  入力信号
			double y=0;	//!<  出力信号
			double uZ1=0;
			double yZ1=0;

			// ベクトル版の変数
			Matrix<1,N> gv;	//!< [rad/s] 遮断周波数
			Matrix<1,N> uZ1v;	//!< [-] 入力過去値
			Matrix<1,N> yZ1v;	//!< [-] 出力過去値
			Matrix<1,N> uv;	//!< [-] 入力ベクトル
			Matrix<1,N> yv;	//!< [-] 出力ベクトル
	};
}

#endif

//! @file TransferFunction.hh
//! @brief 伝達関数クラス
//!
//! 任意の伝達関数を保持，入力に対する応答を計算して出力するクラス。
//! (MATLABでいうところの「tf」のようなもの)
//! 注意：プロパー(相対次数0以上)のみサポート
//!
//! @date 2022/03/29
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2022 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef TRANSFERFUNCTION
#define TRANSFERFUNCTION

#include <cassert>
#include "Matrix.hh"
#include "StateSpaceSystem.hh"

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
//! @brief 伝達関数クラス
//! @tparam N	分子次数
//! @tparam	D	分母次数
template <size_t N, size_t D>
class TransferFunction {
	public:
		//! @brief 空コンストラクタ
		TransferFunction(void)
			: Sys()
		{
			PassedLog();
		}
		
		//! @brief コンストラクタ
		//! @param[in]	Num	分子の係数ベクトル 例：(b1*s + b0) のとき Matrix<1,2> Num = {b1, b0}
		//! @param[in]	Den	分母の係数ベクトル 例：(a2*s^2 + a1*s + a0) のとき Matrix<1,3> Den = {a2, a1, a0}
		//! @param[in]	SmplTime	サンプリング周期 [s]
		TransferFunction(const Matrix<1,N+1>& Num, const Matrix<1,D+1>& Den, const double SmplTime)
			: Sys()
		{
			static_assert(N <= D);				// プロパーかどうかのチェック
			SetCoefficients(Num, Den, SmplTime);// 伝達関数の係数にセット
			PassedLog();
		}

		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		TransferFunction(TransferFunction&& r)
			: Sys(r.Sys)
		{
			
		}

		//! @brief デストラクタ
		~TransferFunction(){
			PassedLog();
		}

		//! @brief 伝達関数の係数を設定する関数
		//! @param[in]	Num	分子の係数ベクトル 例：(b1*s + b0) のとき Matrix<1,2> Num = {b1, b0}
		//! @param[in]	Den	分母の係数ベクトル 例：(a2*s^2 + a1*s + a0) のとき Matrix<1,3> Den = {a2, a1, a0}
		//! @param[in]	Ts	サンプリング周期 [s]
		void SetCoefficients(const Matrix<1,N+1>& Num, const Matrix<1,D+1>& Den, const double SmplTime){
			// 分母の最上位係数を1に変形
			const Matrix<1,N+1> b_n = Num/Den[1];
			const Matrix<1,D+1> a_d = Den/Den[1];
			
			// 可制御正準系の連続系状態空間モデルの作成
			// A行列の生成
			Matrix<D,D> A;		// 連続系A行列
			for(size_t i = 1; i < D; ++i){
				A.SetElement(i + 1, i, 1);
			}
			for(size_t i = 1; i <= D; ++i){
				A.SetElement(i, D, -a_d[D + 2 - i]);
			}
			
			// B行列の生成
			Matrix<1,D> b;		// 連続系bベクトル
			b[D] = 1;
			
			// C行列とD行列の生成
			if constexpr(N != D){
				// 直達項が無い，相対次数が1以上の場合
				// C行列のみ生成
				Matrix<D,1> c;	// cベクトル
				for(size_t i = 1; i <= N + 1; ++i){
					c.SetElement(i, 1, b_n[N + 2 - i]);
				}
				
				// 状態空間モデルに設定
				Sys.SetContinuous(A, b, c, SmplTime);
			}else{
				// 直達項が有る，相対次数が0の場合
				
				// C行列の生成
				Matrix<D,1> c;	// cベクトル
				for(size_t i = 1; i <= D; ++i){
					c.SetElement(i, 1, b_n[D + 2 - i] - a_d[D + 2 - i]*b_n[1]);
				}
				
				// D行列の生成
				Matrix<1,1> d;	// dベクトル
				d[1] = b_n[1];
				
				// 状態空間モデルに設定
				Sys.SetContinuous(A, b, c, d, SmplTime);
			}
		}
		
		//! @brief 入力信号に対する伝達関数の応答を返す関数(1サンプル遅れ無し)
		//! @param[in]	u	入力信号
		double GetResponse(const double u){
			return Sys.GetNextResponse(u);	// 1サンプル遅れの無い[k+1]の応答を即時に返す
		}
		
		//! @brief 入力信号に対する伝達関数の厳密な応答を返す関数(1サンプル遅れ有り)
		//! @param[in]	u	入力信号
		double GetStrictResponse(const double u){
			return Sys.GetResponse(u);		// 1サンプル遅れの有る[k]の厳密な応答を返す
		}
		
	private:
		TransferFunction(const TransferFunction&) = delete;					//!< コピーコンストラクタ使用禁止
		const TransferFunction& operator=(const TransferFunction&) = delete;//!< 代入演算子使用禁止
		StateSpaceSystem<D> Sys;	//!< SISO状態空間モデル
};
}

#endif


//! @file StateSpaceSystem.hh
//! @brief 状態空間表現によるシステムクラス
//!
//! 線形の状態空間モデルで表現されたシステムを保持，入力信号に対する出力信号を計算する。
//! (MATLABでいうところの「ss」のようなもの)
//!
//! @date 2022/11/11
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2022 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef STATESPACESYSTEM
#define STATESPACESYSTEM

#include <cassert>
#include "Matrix.hh"
#include "Discret.hh"
#include "ArcsControl.hh"

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
//! @brief 状態空間表現によるシステムクラス
//! @tparam	N	次数
//! @tparam I	入力信号の数
//! @tparam O	出力信号の数
template <size_t N, size_t I = 1, size_t O = 1>
class StateSpaceSystem {
	public:
		//! @brief コンストラクタ(空コンストラクタ版)
		StateSpaceSystem(void)
			: Ad(), Bd(), Cd(), Dd(), x(), x_next(), DirectTerm(true)
		{
			PassedLog();
		}
		
		//! @brief コンストラクタ(連続系A,B,C行列設定版)
		//! @param[in]	A	連続系A行列
		//! @param[in]	B	連続系B行列
		//! @param[in]	C	C行列
		//! @param[in]	Ts	サンプリング周期 [s]
		StateSpaceSystem(const Matrix<N,N>& A, const Matrix<I,N>& B, const Matrix<N,O>& C, const double Ts)
			: Ad(), Bd(), Cd(), Dd(), x(), x_next(), DirectTerm(true)
		{
			SetContinuous(A, B, C, Ts);	// 連続系のA行列，B行列，C行列を設定して離散化
			PassedLog();
		}
		
		//! @brief コンストラクタ(連続系A,B,C,D行列設定版)
		//! @param[in]	A	連続系A行列
		//! @param[in]	B	連続系B行列
		//! @param[in]	C	C行列
		//! @param[in]	D	D行列
		//! @param[in]	Ts	サンプリング周期 [s]
		StateSpaceSystem(const Matrix<N,N>& A, const Matrix<I,N>& B, const Matrix<N,O>& C, const Matrix<I,O>& D, const double Ts)
			: Ad(), Bd(), Cd(), Dd(), x(), x_next(), DirectTerm(true)
		{
			SetContinuous(A, B, C, D, Ts);	// 連続系のA行列，B行列，C行列，D行列を設定して離散化
			PassedLog();
		}
		
		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		StateSpaceSystem(StateSpaceSystem&& r)
			: Ad(r.Ad), Bd(r.Bd), Cd(r.Cd), Dd(r.Dd), x(r.x), x_next(r.x_next), DirectTerm(r.DirectTerm)
		{
			
		}

		//! @brief デストラクタ
		~StateSpaceSystem(){
			PassedLog();
		}
		
		//! @brief 連続系のA行列，B行列，C行列を設定して離散化する関数
		//! @param[in]	A	連続系A行列
		//! @param[in]	B	連続系B行列
		//! @param[in]	C	C行列
		//! @param[in]	Ts	サンプリング周期 [s]
		void SetContinuous(const Matrix<N,N>& A, const Matrix<I,N>& B, const Matrix<N,O>& C, const double Ts){
			if constexpr(N == 1){
				// 次数が1次のときは平衡化せずに離散化
				auto Ah = A;
				auto Bh = B;
				auto Ch = C;
				std::tie(Ad, Bd) = ArcsControl::Discretize(Ah, Bh, Ts);	// 離散化
				Cd = Ch;			// C行列はそのまま
			}else{
				// 次数が1次より高いときは平衡化してから離散化
				auto [Ah, Bh, Ch] = ArcsControl::BalanceReal(A, B, C);	// 平衡化
				std::tie(Ad, Bd) = ArcsControl::Discretize(Ah, Bh, Ts);	// 離散化
				Cd = Ch;			// C行列は平衡化後そのまま
			}
			DirectTerm = false;	// 直達項は無し
		}
		
		//! @brief 連続系のA行列，B行列，C行列，D行列を設定して離散化する関数
		//! @param[in]	A	連続系A行列
		//! @param[in]	B	連続系B行列
		//! @param[in]	C	C行列
		//! @param[in]	D	D行列
		//! @param[in]	Ts	サンプリング周期 [s]
		void SetContinuous(const Matrix<N,N>& A, const Matrix<I,N>& B, const Matrix<N,O>& C, const Matrix<I,O>& D, const double Ts){
			SetContinuous(A, B, C, Ts);
			Dd = D;				// D行列は何もせずそのまま
			DirectTerm = true;	// 直達項は有り
		}
		
		//! @brief 離散系のA行列，B行列，C行列を設定する関数
		//! @param[in]	A	離散系A行列
		//! @param[in]	B	離散系B行列
		//! @param[in]	C	C行列
		void SetDiscrete(const Matrix<N,N>& A, const Matrix<I,N>& B, const Matrix<N,O>& C){
			Ad = A;
			Bd = B;
			Cd = C;
			DirectTerm = false;	// 直達項は無し
		}
		
		//! @brief 離散系のA行列，B行列，C行列，D行列を設定する関数
		//! @param[in]	A	離散系A行列
		//! @param[in]	B	離散系B行列
		//! @param[in]	C	C行列
		//! @param[in]	D	D行列
		void SetDiscrete(const Matrix<N,N>& A, const Matrix<I,N>& B, const Matrix<N,O>& C, const Matrix<I,O>& D){
			Ad = A;
			Bd = B;
			Cd = C;
			Dd = D;
			DirectTerm = true;	// 直達項は有り
		}
		
		//! @brief 状態空間モデルの応答を計算して取得する関数(普通版)
		//! @param[in]	u	入力ベクトル
		//! @param[out]	yout	出力ベクトル
		void GetResponses(const Matrix<1,I>& u, Matrix<1,O>& yout){
			x_next = Ad*x + Bd*u;	// 状態方程式
			if(DirectTerm == false){
				yout = Cd*x;		// 直達項なし版の出力方程式
			}else{
				yout = Cd*x + Dd*u;	// 直達項あり版の出力方程式
			}
			x = x_next;				// 状態ベクトルを更新
		}
		
		//! @brief 状態空間モデルの応答を計算して取得する関数(ベクトルで返す版)
		//! @param[in]	u	入力ベクトル
		//! @return	出力ベクトル
		Matrix<1,O> GetResponses(const Matrix<1,I>& u){
			Matrix<1,O> y;			// 出力ベクトル
			GetResponses(u, y);		// 応答計算
			return y;				// 出力ベクトルを返す
		}
		
		//! @brief 状態空間モデルの応答を計算して取得する関数(SISOでスカラーで返す版)
		//! @param[in]	u	入力
		//! @return	出力
		double GetResponse(const double u){
			static_assert(I == 1 && O == 1);// SISOかチェック
			Matrix<1,I> u_vec;				// 入力ベクトル
			Matrix<1,O> y_vec;				// 出力ベクトル
			u_vec[1] = u;
			GetResponses(u_vec, y_vec);		// 応答計算
			return y_vec[1];				// 出力を返す
		}
		
		//! @brief 状態空間モデルの応答を計算して取得する関数(2入力1出力でスカラーで返す版)
		//! @param[in]	u	入力
		//! @return	出力
		double GetResponse(const double u1, const double u2){
			static_assert(I == 2 && O == 1);// 2入力1出力かチェック
			Matrix<1,I> u_vec;				// 入力ベクトル
			Matrix<1,O> y_vec;				// 出力ベクトル
			u_vec[1] = u1;
			u_vec[2] = u2;
			GetResponses(u_vec, y_vec);		// 応答計算
			return y_vec[1];				// 出力を返す
		}
		
		//! @brief 状態空間モデルの応答を計算して取得する関数(次の時刻の出力ベクトルを即時に返す版)
		//! @param[in]	u	入力ベクトル
		//! @param[out]	yout	出力ベクトル
		void GetNextResponses(const Matrix<1,I>& u, Matrix<1,O>& yout){
			x_next = Ad*x + Bd*u;		// 状態方程式
			if(DirectTerm == false){
				yout = Cd*x_next;		// 直達項なし版の出力方程式(次の時刻の出力ベクトルを即時に返す)
			}else{
				yout = Cd*x_next + Dd*u;// 直達項あり版の出力方程式(次の時刻の出力ベクトルを即時に返す)
			}
			x = x_next;					// 状態ベクトルを更新
		}
		
		//! @brief 状態空間モデルの応答を計算して取得する関数(次の時刻の出力ベクトルを即時に返す版)(ベクトルで返す版)
		//! @param[in]	u	入力ベクトル
		//! @return	出力ベクトル
		Matrix<1,O> GetNextResponses(const Matrix<1,I>& u){
			Matrix<1,O> y;			// 出力ベクトル
			GetNextResponses(u, y);	// 応答計算
			return y;				// 出力ベクトルを返す
		}
		
		//! @brief 状態空間モデルの応答を計算して取得する関数(次の時刻の出力ベクトルを即時に返す版)(SISOでスカラーで返す版)
		//! @param[in]	u	入力
		//! @return	出力
		double GetNextResponse(const double u){
			static_assert(I == 1 && O == 1);// SISOかチェック
			Matrix<1,I> u_vec;				// 入力ベクトル
			Matrix<1,O> y_vec;				// 出力ベクトル
			u_vec[1] = u;
			GetNextResponses(u_vec, y_vec);	// 応答計算
			return y_vec[1];				// 出力を返す
		}
		
		//! @brief 状態空間モデルの応答を計算して取得する関数(次の時刻の出力ベクトルを即時に返す版)(2入力1出力でスカラーで返す版)
		//! @param[in]	u	入力
		//! @return	出力
		double GetNextResponse(const double u1, const double u2){
			static_assert(I == 2 && O == 1);// 2入力1出力かチェック
			Matrix<1,I> u_vec;				// 入力ベクトル
			Matrix<1,O> y_vec;				// 出力ベクトル
			u_vec[1] = u1;
			u_vec[2] = u2;
			GetNextResponses(u_vec, y_vec);	// 応答計算
			return y_vec[1];				// 出力を返す
		}
		
		//! @brief 状態ベクトルをクリアする関数
		void ClearStateVector(void){
			x = Matrix<1,N>::zeros();
		}
		
	private:
		StateSpaceSystem(const StateSpaceSystem&) = delete;					//!< コピーコンストラクタ使用禁止
		const StateSpaceSystem& operator=(const StateSpaceSystem&) = delete;//!< 代入演算子使用禁止
		Matrix<N,N> Ad;		//!< 離散系A行列
		Matrix<I,N> Bd;		//!< 離散系B行列
		Matrix<N,O> Cd;		//!< C行列
		Matrix<I,O> Dd;		//!< D行列
		Matrix<1,N> x;		//!< 状態ベクトル
		Matrix<1,N> x_next;	//!< 次の時刻の状態ベクトル
		bool DirectTerm;	//!< 直達項の有無フラグ(true = 直達項あり，false = 直達項なし)
};
}

#endif


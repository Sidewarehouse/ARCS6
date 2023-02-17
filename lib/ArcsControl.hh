//! @file ArcsControl.hh
//! @brief ARCS-Control 制御理論クラス
//!
//! 制御理論に関係する様々なアルゴリズムを詰め合わせた静的クラス
//!
//! @date 2022/07/27
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2022 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef ARCSCONTROL
#define ARCSCONTROL

#include <cassert>
#include <tuple>
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
//! @brief ARCS-Control 制御理論クラス
class ArcsControl {
	public:
		//! @brief 連続リアプノフ方程式 A*X + X*A^T + Q = 0 の解Xを求める関数(実数版, 引数で返す版)
		//! @tparam		NN	行列の幅
		//! @tparam		MM	行列の高さ
		//! @param[in]	A	A行列
		//! @param[in]	Q	Q行列
		//! @param[out]	X	解Xの行列
		template<size_t NN, size_t MM>
		static constexpr void Lyapunov(const Matrix<NN,MM,double>& A, const Matrix<NN,MM,double>& Q, Matrix<NN,MM,double>& X){
			static_assert(NN == MM);					// 正方行列のみ対応
			constexpr auto I = Matrix<NN,NN>::eye();	// 単位行列
			X = Matrix<1,NN*MM>::template vecinv<NN,MM>( -inv( Kronecker(I,A) + Kronecker(A,I) )*vec(Q) );	// A*X + X*A^T + Q = 0 をXについて解く
		}
		
		//! @brief 連続リアプノフ方程式 A*X + X*A^T + Q = 0 の解Xを求める関数(実数版, 戻り値で返す版)
		//! @tparam		NN	行列の幅
		//! @tparam		MM	行列の高さ
		//! @param[in]	A	A行列
		//! @param[in]	Q	Q行列
		//! @return		解Xの行列
		template<size_t NN, size_t MM>
		static constexpr Matrix<NN,MM,double> Lyapunov(const Matrix<NN,MM,double>& A, const Matrix<NN,MM,double>& Q){
			Matrix<NN,MM,double> X;
			Lyapunov(A, Q, X);
			return X;
		}
		
		//! @brief 可制御グラミアンを計算する関数(引数で返す版)
		//! @tparam		NA	A行列の幅
		//! @tparam		MA	A行列の高さ
		//! @tparam		NB	b行列の幅
		//! @param[in]	A	A行列
		//! @param[in]	b	b行列
		//! @param[out]	Wc	可制御グラミアン
		template<size_t NA, size_t MA, size_t NB>
		static constexpr void GramianCtrl(const Matrix<NA,MA>& A, const Matrix<NB,MA>& b, Matrix<NA,MA>& Wc){
			Lyapunov(A, b*tp(b), Wc);
		}
		
		//! @brief 可制御グラミアンを計算する関数(戻り値で返す版)
		//! @tparam		NA	A行列の幅
		//! @tparam		MA	A行列の高さ
		//! @tparam		NB	b行列の幅
		//! @param[in]	A	A行列
		//! @param[in]	b	b行列
		//! @return		可制御グラミアン
		template<size_t NA, size_t MA, size_t NB>
		static constexpr Matrix<NA,MA> GramianCtrl(const Matrix<NA,MA>& A, const Matrix<NB,MA>& b){
			Matrix<NA,MA> Wc;
			GramianCtrl(A, b, Wc);
			return Wc;
		}
		
		//! @brief 可観測グラミアンを計算する関数(引数で返す版)
		//! @tparam		NA	A行列の幅
		//! @tparam		MA	A行列の高さ
		//! @tparam		MC	c行列の高さ
		//! @param[in]	A	A行列
		//! @param[in]	c	c行列
		//! @param[out]	Wo	可観測グラミアン
		template<size_t NA, size_t MA, size_t MC>
		static constexpr void GramianObsrv(const Matrix<NA,MA>& A, const Matrix<NA,MC>& c, Matrix<NA,MA>& Wo){
			Lyapunov(tp(A), tp(c)*c, Wo);
		}
		
		//! @brief 可観測グラミアンを計算する関数(戻り値で返す版)
		//! @tparam		NA	A行列の幅
		//! @tparam		MA	A行列の高さ
		//! @tparam		MC	c行列の高さ
		//! @param[in]	A	A行列
		//! @param[in]	c	c行列
		//! @return		可観測グラミアン
		template<size_t NA, size_t MA, size_t MC>
		static constexpr Matrix<NA,MA> GramianObsrv(const Matrix<NA,MA>& A, const Matrix<NA,MC>& c){
			Matrix<NA,MA> Wo;
			GramianObsrv(A, c, Wo);
			return Wo;
		}
		
		//! @brief 状態空間モデルを平衡化する関数(引数で返す版)
		//! @tparam		NA	A行列の幅
		//! @tparam		MA	A行列の高さ
		//! @tparam		NB	b行列の幅
		//! @tparam		MC	c行列の高さ
		//! @param[in]	A	A行列
		//! @param[in]	b	b行列
		//! @param[in]	c	c行列
		//! @param[out]	Ah	平衡化後のA行列
		//! @param[out]	bh	平衡化後のb行列
		//! @param[out]	ch	平衡化後のc行列
		template<size_t NA, size_t MA, size_t NB, size_t MC>
		static constexpr void BalanceReal(const Matrix<NA,MA>& A, const Matrix<NB,MA>& b, const Matrix<NA,MC>& c, Matrix<NA,MA>& Ah, Matrix<NB,MA>& bh, Matrix<NA,MC>& ch){
			static_assert(NA == MA);	// A行列は正方行列
			
			// 1. グラミアンズの計算
			const Matrix<NA,MA> Wc = GramianCtrl(A, b);		// 可制御グラミアン
			const Matrix<NA,MA> Wo = GramianObsrv(A, c);	// 可観測グラミアン
			
			// 2. 1をコレスキー分解
			Matrix<NA,MA> Lc, Lo;
			Cholesky(Wc, Lc);
			Cholesky(Wo, Lo);
			
			// 3. 2を特異値分解する
			const auto [U, S, V] = SVD( tp(Lo)*Lc );
			
			// 4. 3から変換行列を計算する
			const Matrix<NA,MA> T  = inv(sqrte(S))*tp(U)*tp(Lo);
			const Matrix<NA,MA> Ti = Lc*V*inv(sqrte(S));
			
			// 5. 4を使って状態空間モデルを変換する
			Ah = T*A*Ti;
			bh = T*b;
			ch = c*Ti;
			
			// 下記はデバッグ用
			//PrintMat(T*Wc*tp(T) - tp(Ti)*Wo*Ti);	// 零行列になればOK
			//PrintMat(T*Ti);						// 単位行列になればOK
		}
		
		//! @brief 状態空間モデルを平衡化する関数(タプルで返す版)
		//! @tparam		NA	A行列の幅
		//! @tparam		MA	A行列の高さ
		//! @tparam		NB	b行列の幅
		//! @tparam		MC	c行列の高さ
		//! @param[in]	A	A行列
		//! @param[in]	b	b行列
		//! @param[in]	c	c行列
		//! @return		{Ah, bh, ch} = {平衡化後のA行列, 平衡化後のb行列, 平衡化後のc行列}
		template<size_t NA, size_t MA, size_t NB, size_t MC>
		static constexpr std::tuple<Matrix<NA,MA>, Matrix<NB,MA>, Matrix<NA,MC>> BalanceReal(const Matrix<NA,MA>& A, const Matrix<NB,MA>& b, const Matrix<NA,MC>& c){
			Matrix<NA,MA> Ah;
			Matrix<NB,MA> bh;
			Matrix<NA,MC> ch;
			BalanceReal(A, b, c, Ah, bh, ch);
			return {Ah, bh, ch};
		}
		
		//! @brief 連続系状態空間モデルを離散化する関数(引数で返す版)
		//! @tparam		N	A行列の幅と高さ
		//! @tparam		M	b行列の幅
		//! @param[in]	Ac	連続系のA行列
		//! @param[in]	Bc	連続系のB行列
		//! @param[out]	Ad	離散系のA行列
		//! @param[out]	Bd	離散系のB行列
		//! @param[in]	Ts	[s] サンプリング時間
		//! @param[in]	Npade	[-] パデ近似の次数(省略可：デフォルト値 13次)
		//! @param[in]	Nint	[-] 定積分の分割数(省略可：デフォルト値 10000点)
		template <size_t N, size_t M>
		static constexpr void Discretize(
			const Matrix<N,N>& Ac, const Matrix<M,N>& Bc, Matrix<N,N>& Ad, Matrix<M,N>& Bd,	const double Ts, const size_t Npade = 13, const size_t Nint = 10000
		){
			Ad = expm(Ac*Ts, Npade);					// A行列の離散化
			Bd = integral_expm(Ac, Ts, Nint, Npade)*Bc;	// B行列の離散化
		}
		
		//! @brief 連続系状態空間モデルを離散化する関数(タプルで返す版)
		//! @tparam		N	A行列の幅と高さ
		//! @tparam		M	b行列の幅
		//! @param[in]	Ac	連続系のA行列
		//! @param[in]	Bc	連続系のB行列
		//! @param[in]	Ts	[s] サンプリング時間
		//! @param[in]	Npade	[-] パデ近似の次数(省略可：デフォルト値 13次)
		//! @param[in]	Nint	[-] 定積分の分割数(省略可：デフォルト値 10000点)
		//! @return		{Ad, Bd} = {離散系のA行列, 離散系のB行列}
		template <size_t N, size_t M>
		static constexpr std::tuple<Matrix<N,N>, Matrix<M,N>> Discretize(
			const Matrix<N,N>& Ac, const Matrix<M,N>& Bc, const double Ts, const size_t Npade = 13, const size_t Nint = 10000
		){
			Matrix<N,N> Ad;
			Matrix<M,N> Bd;
			Discretize(Ac, Bc, Ad, Bd, Ts, Npade, Nint);
			return {Ad, Bd};
		}
		
	private:
		ArcsControl() = delete;						//!< コンストラクタ使用禁止
		ArcsControl(ArcsControl&& r) = delete;		//!< ムーブコンストラクタ使用禁止
		~ArcsControl() = delete;					//!< デストラクタ使用禁止
		ArcsControl(const ArcsControl&) = delete;	//!< コピーコンストラクタ使用禁止
		const ArcsControl& operator=(const ArcsControl&) = delete;	//!< 代入演算子使用禁止
};
}

#endif


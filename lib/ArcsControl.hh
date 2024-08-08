//! @file ArcsControl.hh
//! @brief ArcsControl-制御理論ヘッダ
//!
//! 制御理論に関係する様々なアルゴリズムを詰め合わせたヘッダ
//!
//! @date 2024/08/08
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#ifndef ARCSCONTROL
#define ARCSCONTROL

#include <cassert>
#include <tuple>
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

using namespace ARCS::ArcsMatrix;

namespace ARCS {		// ARCS名前空間
namespace ArcsControl {	// ArcsControl名前空間
	//! @brief 連続リアプノフ方程式 AX + XA' + Q = 0 の解Xを求める関数 (引数渡し版)
	//! @tparam		M,MQ,MX	行列の高さ
	//! @tparam		N,NQ,NX	行列の幅
	//! @tparam		T,TQ,TX	行列のデータ型
	//! @param[in]	A	A行列
	//! @param[in]	Q	Q行列
	//! @param[out]	X	解Xの行列
	template<size_t M, size_t N, typename T = double, size_t MQ, size_t NQ, typename TQ = double, size_t MX, size_t NX, typename TX = double>
	static constexpr void Lyapunov(const ArcsMat<M,N,T>& A, const ArcsMat<MQ, NQ, TQ>& Q, ArcsMat<MX,NX,TX>& X){
		static_assert(M == N,   "ArcsCtrl: Size Error");	// 正方行列のみ対応
		static_assert(MQ == NQ, "ArcsCtrl: Size Error");	// 正方行列のみ対応
		static_assert(MX == NX, "ArcsCtrl: Size Error");	// 正方行列のみ対応
		static_assert(MQ == M,  "ArcsCtrl: Size Error");	// サイズチェック
		static_assert(MX == M,  "ArcsCtrl: Size Error");	// サイズチェック
		static_assert(NQ == N,  "ArcsCtrl: Size Error");	// サイズチェック
		static_assert(NX == N,  "ArcsCtrl: Size Error");	// サイズチェック

		// 連続リアプノフ方程式は以下のように線形方程式に変形できるので、
		// 　AX + XA' + Q = 0
		// 　AX + XA' = -Q
		// 　vec(AX + XA') = vec(-Q)
		// 　vec(AX) + vec(XA') = -vec(Q)
		// 　vec(AX) + vec(XA') = -vec(Q)
		// 　(I ox A)vec(X) + (A ox I)vec(X) = -vec(Q)
		// 　{(I ox A) + (A ox I)}vec(X) = -vec(Q)
		// 　Cy = d
		// の線形方程式のyについて解いて、vec作用素の逆を取れば、
		// 　X = vec^(-1)( y )
		// のように解Xが求まる。
		const auto I = ArcsMat<M,N,T>::eye();	// 単位行列
		const auto C = Kron(I, A) + Kron(A, I);	// "ox"はクロネッカー積
		const auto d = -vec(Q);
		X = vecinv<M,N>( linsolve(C, d) );		// 線形方程式 Cy = d を y について解いて逆vec作用素を掛けるだけ
	}

	//! @brief 連続リアプノフ方程式 AX + XA' + Q = 0 の解Xを求める関数 (戻り値返し版)
	//! @tparam		M,MQ	行列の高さ
	//! @tparam		N,NQ	行列の幅
	//! @tparam		T,TQ	行列のデータ型
	//! @param[in]	A	A行列
	//! @param[in]	Q	Q行列
	//! @return	解Xの行列
	template<size_t M, size_t N, typename T = double, size_t MQ, size_t NQ, typename TQ = double>
	static constexpr ArcsMat<M,N,T> Lyapunov(const ArcsMat<M,N,T>& A, const ArcsMat<MQ, NQ, TQ>& Q){
		ArcsMat<M,N,T> X;
		Lyapunov(A, Q, X);
		return X;
	}

	//! @brief 離散リアプノフ方程式 AXA' - X + Q = 0 の解Xを求める関数 (引数渡し版)
	//! @tparam		M,MQ,MX	行列の高さ
	//! @tparam		N,NQ,NX	行列の幅
	//! @tparam		T,TQ,TX	行列のデータ型
	//! @param[in]	A	A行列
	//! @param[in]	Q	Q行列
	//! @param[out]	X	解Xの行列
	template<size_t M, size_t N, typename T = double, size_t MQ, size_t NQ, typename TQ = double, size_t MX, size_t NX, typename TX = double>
	static constexpr void DiscLyapunov(const ArcsMat<M,N,T>& A, const ArcsMat<MQ, NQ, TQ>& Q, ArcsMat<MX,NX,TX>& X){
		static_assert(M == N,   "ArcsCtrl: Size Error");	// 正方行列のみ対応
		static_assert(MQ == NQ, "ArcsCtrl: Size Error");	// 正方行列のみ対応
		static_assert(MX == NX, "ArcsCtrl: Size Error");	// 正方行列のみ対応
		static_assert(MQ == M,  "ArcsCtrl: Size Error");	// サイズチェック
		static_assert(MX == M,  "ArcsCtrl: Size Error");	// サイズチェック
		static_assert(NQ == N,  "ArcsCtrl: Size Error");	// サイズチェック
		static_assert(NX == N,  "ArcsCtrl: Size Error");	// サイズチェック

		// 離散リアプノフ方程式は以下のように線形方程式に変形できるので、
		// 　AXA' - X + Q = 0
		// 　AXA' - X = -Q
		// 　vec(AXA' - X) = vec(-Q)
		// 　vec(AXA') - vec(X) = -vec(Q)
		// 　(A ox A)vec(X) - vec(X) = -vec(Q)
		// 　{(A ox A) - I}vec(X) = -vec(Q)
		// 　Cy = d
		// の線形方程式のyについて解いて、vec作用素の逆を取れば、
		// 　X = vec^(-1)( y )
		// のように解Xが求まる。
		const auto I = ArcsMat<M*M,N*N,T>::eye();	// 単位行列
		const auto C = Kron(A, A) - I;
		const auto d = -vec(Q);
		X = vecinv<M,N>( linsolve(C, d) );			// 線形方程式 Cy = d を y について解いて逆vec作用素を掛けるだけ
	}

	//! @brief 離散リアプノフ方程式 AXA' - X + Q = 0 の解Xを求める関数 (戻り値返し版)
	//! @tparam		M,MQ	行列の高さ
	//! @tparam		N,NQ	行列の幅
	//! @tparam		T,TQ	行列のデータ型
	//! @param[in]	A	A行列
	//! @param[in]	Q	Q行列
	//! @return	解Xの行列
	template<size_t M, size_t N, typename T = double, size_t MQ, size_t NQ, typename TQ = double>
	static constexpr ArcsMat<M,N,T> DiscLyapunov(const ArcsMat<M,N,T>& A, const ArcsMat<MQ, NQ, TQ>& Q){
		ArcsMat<M,N,T> X;
		DiscLyapunov(A, Q, X);
		return X;
	}

	//! @brief 可制御グラミアンを計算する関数 (引数渡し版)
	//! @tparam		M	A行列の高さ
	//! @tparam		N	A行列の幅
	//! @tparam		MB	B行列の高さ
	//! @tparam		NB	B行列の幅
	//! @param[in]	A	A行列
	//! @param[in]	B	B行列
	//! @param[out]	Wc	可制御グラミアン
	template<size_t M, size_t N, typename T = double, size_t MB, size_t NB, typename TB = double, size_t MW, size_t NW, typename TW = double>
	static constexpr void GramianCtrl(const ArcsMat<M,N,T>& A, const ArcsMat<MB,NB,TB>& B, ArcsMat<MW,NW,TW>& Wc){
		static_assert(M == N,   "ArcsCtrl: Size Error");	// 正方行列のみ対応
		static_assert(MB == M , "ArcsCtrl: Size Error");	// サイズチェック
		static_assert(MW == M , "ArcsCtrl: Size Error");	// サイズチェック
		static_assert(NW == N , "ArcsCtrl: Size Error");	// サイズチェック
		Lyapunov(A, B*~B, Wc);
	}
	
	//! @brief 可制御グラミアンを計算する関数 (戻り値返し版)
	//! @tparam		M	A行列の高さ
	//! @tparam		N	A行列の幅
	//! @tparam		MB	B行列の高さ
	//! @tparam		NB	B行列の幅
	//! @param[in]	A	A行列
	//! @param[in]	B	B行列
	//! @return	可制御グラミアン
	template<size_t M, size_t N, typename T = double, size_t MB, size_t NB, typename TB = double>
	static constexpr ArcsMat<M,N,T> GramianCtrl(const ArcsMat<M,N,T>& A, const ArcsMat<MB,NB,TB>& B){
		ArcsMat<M,N,T> Wc;
		GramianCtrl(A, B, Wc);
		return Wc;
	}

	//! @brief 可観測グラミアンを計算する関数 (引数渡し版)
	//! @tparam		M	A行列の高さ
	//! @tparam		N	A行列の幅
	//! @tparam		MC	C行列の高さ
	//! @tparam		NC	C行列の幅
	//! @param[in]	A	A行列
	//! @param[in]	C	C行列
	//! @param[out]	Wo	可観測グラミアン
	template<size_t M, size_t N, typename T = double, size_t MC, size_t NC, typename TC = double, size_t MW, size_t NW, typename TW = double>
	static constexpr void GramianObsrv(const ArcsMat<M,N,T>& A, const ArcsMat<MC,NC,TC>& C, ArcsMat<MW,NW,TW>& Wo){
		static_assert(M == N,   "ArcsCtrl: Size Error");	// 正方行列のみ対応
		static_assert(NC == M , "ArcsCtrl: Size Error");	// サイズチェック
		static_assert(MW == M , "ArcsCtrl: Size Error");	// サイズチェック
		static_assert(NW == N , "ArcsCtrl: Size Error");	// サイズチェック
		Lyapunov(~A, ~C*C, Wo);
	}
	
	//! @brief 可観測グラミアンを計算する関数 (戻り値返し版)
	//! @tparam		M	A行列の高さ
	//! @tparam		N	A行列の幅
	//! @tparam		MC	C行列の高さ
	//! @tparam		NC	C行列の幅
	//! @param[in]	A	A行列
	//! @param[in]	C	C行列
	//! @return	可観測グラミアン
	template<size_t M, size_t N, typename T = double, size_t MC, size_t NC, typename TC = double>
	static constexpr ArcsMat<M,N,T> GramianObsrv(const ArcsMat<M,N,T>& A, const ArcsMat<MC,NC,TC>& C){
		ArcsMat<M,N,T> Wo;
		GramianObsrv(A, C, Wo);
		return Wo;
	}
	
	//! @brief 状態空間モデルを平衡化する関数 (引数渡し版)
	//! @tparam		M,MH	A行列の高さ
	//! @tparam		N,NH	A行列の幅
	//! @tparam		MB,MBH	B行列の高さ
	//! @tparam		NB,NBH	B行列の幅
	//! @tparam		MC,MCH	C行列の高さ
	//! @tparam		NC,NCH	C行列の幅
	//! @param[in]	A	A行列
	//! @param[in]	B	B行列
	//! @param[in]	C	C行列
	//! @param[out]	Ah	平衡化後のA行列
	//! @param[out]	Bh	平衡化後のB行列
	//! @param[out]	Ch	平衡化後のC行列
	template<
		size_t M, size_t N, typename T = double, size_t MB, size_t NB, typename TB = double, size_t MC, size_t NC, typename TC = double,
		size_t MH, size_t NH, typename TH = double, size_t MBH, size_t NBH, typename TBH = double, size_t MCH, size_t NCH, typename TCH = double
	>
	static constexpr void BalanceReal(const ArcsMat<M,N,T>& A, const ArcsMat<MB,NB,TB>& B, const ArcsMat<MC,NC,TC>& C, ArcsMat<MH,NH,TH>& Ah, ArcsMat<MBH,NBH,TBH>& Bh, ArcsMat<MCH,NCH,TCH>& Ch){
		static_assert(M == N,   "ArcsCtrl: Size Error");	// A行列は正方行列
		static_assert(MH == NH, "ArcsCtrl: Size Error");	// Ah行列は正方行列
		static_assert(MB == M , "ArcsCtrl: Size Error");	// サイズチェック
		static_assert(NC == M , "ArcsCtrl: Size Error");	// サイズチェック
		static_assert(MH == M , "ArcsCtrl: Size Error");	// サイズチェック
		static_assert(MBH == MB , "ArcsCtrl: Size Error");	// サイズチェック
		static_assert(NBH == NB , "ArcsCtrl: Size Error");	// サイズチェック
		static_assert(MCH == MC , "ArcsCtrl: Size Error");	// サイズチェック
		static_assert(NCH == NC , "ArcsCtrl: Size Error");	// サイズチェック

		// 参考文献： ALAN J. LAUB, MICHAEL T. HEATH, CHRIS C. PAIGE, AND RO3ERT C. WARD,
		//           "Computation of System Balancing Transformations 115 - and Other Applications of Simultaneous Diagonalization Algorithms,"
		//           IEEE Trans. Autom. Cont. vol. AC-32, no.2, Feb. 1987.
		// 1. グラミアンの計算
		const auto Wc = GramianCtrl(A, B);	// 可制御グラミアン
		const auto Wo = GramianObsrv(A, C);	// 可観測グラミアン

		// 2. 1をコレスキー分解
		const auto Lc = ~Cholesky(Wc);		// 論文では下三角としているので、転置しておく
		const auto Lo = ~Cholesky(Wo);		// 論文では下三角としているので、転置しておく

		// 3. 2を特異値分解する
		const auto [U, S, V] = SVD(~Lo*Lc);

		// 4. 3から変換行列を計算する
		const auto l = ArcsMat<M,1,T>::ones();
		const auto d = l % sqrt(getdiag(S));
		const auto Si = diag(d);
		const auto P  = Si*~U*~Lo;
		const auto Pi = Lc*V*Si;

		// 5. 4を使って状態空間モデルを変換する
		Ah = P*A*Pi;
		Bh = P*B;
		Ch = C*Pi;
		
		// 下記はデバッグ用
		//disp(P*Wc*~P - ~Pi*Wo*Pi);	// 零行列になればOK
		//disp(P*Pi);					// 単位行列になればOK
	}

	//! @brief 状態空間モデルを平衡化する関数 (タプル返し版)
	//! @tparam		M	A行列の高さ
	//! @tparam		N	A行列の幅
	//! @tparam		MB	B行列の高さ
	//! @tparam		NB	B行列の幅
	//! @tparam		MC	C行列の高さ
	//! @tparam		NC	C行列の幅
	//! @param[in]	A	A行列
	//! @param[in]	B	B行列
	//! @param[in]	C	C行列
	//! @return	(Ah, Bh, Ch)	平衡化後のA行列、B行列、C行列のタプル
	template<size_t M, size_t N, typename T = double, size_t MB, size_t NB, typename TB = double, size_t MC, size_t NC, typename TC = double>
	static constexpr std::tuple<ArcsMat<M,N,T>, ArcsMat<MB,NB,TB>, ArcsMat<MC,NC,TC>> BalanceReal(const ArcsMat<M,N,T>& A, const ArcsMat<MB,NB,TB>& B, const ArcsMat<MC,NC,TC>& C){
		ArcsMat<M,N,T> Ah;
		ArcsMat<MB,NB,TB> Bh;
		ArcsMat<MC,NC,TC> Ch;
		BalanceReal(A, B, C, Ah, Bh, Ch);
		return {Ah, Bh, Ch};
	}
	
	//! @brief 連続系状態空間モデルを離散化する関数 (引数渡し版)
	//! @tparam	Npade	[-] パデ近似の次数 (デフォルト値 = 13次)
	//! @tparam	Nint	[-] 定積分の分割数 (デフォルト値 = 10000点)
	//! @tparam	M	Ac行列の高さ
	//! @tparam	N	Ac行列の幅
	//! @tparam	MB	Bc行列の高さ
	//! @tparam	NB	Bc行列の幅
	//! @tparam	MAD	Ad行列の高さ
	//! @tparam	NAD	Ad行列の幅
	//! @tparam	MBD	Bd行列の高さ
	//! @tparam	NBD	Bd行列の幅
	//! @param[in]	Ac	連続系のA行列
	//! @param[in]	Bc	連続系のB行列
	//! @param[out]	Ad	離散系のA行列
	//! @param[out]	Bd	離散系のB行列
	//! @param[in]	Ts	[s] サンプリング時間
	template<
		size_t Npade = 13, size_t Nintg = 10000,
		size_t M, size_t N, typename T = double, size_t MB, size_t NB, typename TB = double,
		size_t MAD, size_t NAD, typename TAD = double, size_t MBD, size_t NBD, typename TBD = double
	>
	static constexpr void Discretize(
		const ArcsMat<M,N,T>& Ac, const ArcsMat<MB,NB,TB>& Bc, ArcsMat<MAD,NAD,TAD>& Ad, ArcsMat<MBD,NBD,TBD>& Bd, const double Ts
	){
		static_assert(M == N,   "ArcsCtrl: Size Error");	// A行列は正方行列
		static_assert(MAD == NAD, "ArcsCtrl: Size Error");	// Ah行列は正方行列
		static_assert(MB == M , "ArcsCtrl: Size Error");	// サイズチェック
		static_assert(MAD == M , "ArcsCtrl: Size Error");	// サイズチェック
		static_assert(MBD == MB , "ArcsCtrl: Size Error");	// サイズチェック
		static_assert(NBD == NB , "ArcsCtrl: Size Error");	// サイズチェック
		static_assert(0 < Npade, "ArcsCtrl: Setting Error");// 設定チェック
		static_assert(0 < Nintg, "ArcsCtrl: Setting Error");// 設定チェック
		arcs_assert(0 < Ts);	// 設定チェック

		expm<Npade>(Ac*Ts, Ad);					// A行列の離散化

		// シンプソン法によるexpmの定積分
		const auto I = ArcsMat<M,N,T>::eye();	// 単位行列
		const T h = Ts/static_cast<T>(2*Nintg);	// [s] 時間ステップ
		T t = 0;								// [s] 時刻
		ArcsMat<M,N,T> S1, S2, intg_expm;
		for(size_t j = 1; j <= Nintg; ++j){
			t = h*static_cast<T>(2*j - 1);		// [s] 奇数時間ステップ
			S1 += expm<Npade>(Ac*t);			// 奇数分の積分
		}
		for(size_t j = 1; j <= Nintg - 1; ++j){
			t = h*static_cast<T>(2*j);			// [s] 偶数時間ステップ
			S2 += expm<Npade>(Ac*t);			// 偶数分の積分
		}
		intg_expm = h/3.0*( I + 4.0*S1 + 2.0*S2 + expm<Npade>(Ac*Ts) );	// 最終的な定積分結果
		Bd = intg_expm*Bc;						// B行列の離散化
	}
	
	//! @brief 連続系状態空間モデルを離散化する関数 (タプル返し版)
	//! @tparam	Npade	[-] パデ近似の次数 (デフォルト値 = 13次)
	//! @tparam	Nint	[-] 定積分の分割数 (デフォルト値 = 10000点)
	//! @tparam	M	Ac行列の高さ
	//! @tparam	N	Ac行列の幅
	//! @tparam	MB	Bc行列の高さ
	//! @tparam	NB	Bc行列の幅
	//! @param[in]	Ac	連続系のA行列
	//! @param[in]	Bc	連続系のB行列
	//! @param[in]	Ts	[s] サンプリング時間
	//! @return	(Ad, Bd)	離散系のA行列とB行列のタプル
	template<
		size_t Npade = 13, size_t Nintg = 10000,
		size_t M, size_t N, typename T = double, size_t MB, size_t NB, typename TB = double
	>
	static constexpr std::tuple<ArcsMat<M,N,T>, ArcsMat<MB,NB,TB>> Discretize(const ArcsMat<M,N,T>& Ac, const ArcsMat<MB,NB,TB>& Bc, const double Ts){
		ArcsMat<M,N,T> Ad;
		ArcsMat<MB,NB,TB> Bd;
		Discretize<Npade, Nintg>(Ac, Bc, Ad, Bd, Ts);
		return {Ad, Bd};
	}

	//! @brief 連続系状態空間モデルのA行列が安定かどうかを返す関数
	//! @param[in]	Ac	A行列
	//! @tparam	M	A行列の高さ
	//! @tparam	N	A行列の幅
	//! @tparam	T	A行列のデータ型
	//! @return	安定性： true = 安定、false = 不安定
	template<size_t M, size_t N, typename T = double>
	static constexpr bool IsStable(const ArcsMat<M,N,T>& Ac){
		static_assert(M == N,   "ArcsCtrl: Size Error");	// A行列は正方行列

		auto v = eig(Ac);	// 固有値を計算
		v.Zeroing();		// 原点にかなり近い極はゼロと扱う
		bool ret = true;
		for(size_t i = 1; i <= M; ++i) ret &= (real(v[i]) < 0);	// すべての固有値が負かどうか調べる
		return ret;
	}

	//! @brief 連続系状態空間モデルのA行列とC行列が可観測かどうかを返す関数
	//! @param[in]	Ac	A行列
	//! @param[in]	C	C行列
	//! @tparam	M	A行列の高さ
	//! @tparam	N	A行列の幅
	//! @tparam	T	A行列のデータ型
	//! @tparam	MC	C行列の高さ
	//! @tparam	NC	C行列の幅
	//! @tparam	TC	C行列のデータ型
	//! @return	可観測性： true = 可観測、false = 不可観測
	template<size_t M, size_t N, typename T = double, size_t MC, size_t NC, typename TC = double>
	static constexpr bool IsObsv(const ArcsMat<M,N,T>& Ac, const ArcsMat<MC,NC,TC>& C){
		static_assert(M == N,   "ArcsCtrl: Size Error");	// A行列は正方行列
		static_assert(NC == M , "ArcsCtrl: Size Error");	// サイズチェック

		ArcsMat<MC,NC,TC> CA;
		ArcsMat<MC*M,NC,TC> Ob;
		for(size_t i = 0; i < M; ++i){
			CA = C*(Ac^i);
			copymatrix(CA,1,MC,1,NC, Ob,i*MC+1,1);	// 可観測性行列を生成
		}
		return M == rank(Ob);	// ランクを計算して状態数と同一なら可観測
	}

	//! @brief 連続系状態空間モデルのA行列とB行列が可制御かどうかを返す関数
	//! @param[in]	Ac	A行列
	//! @param[in]	Bc	B行列
	//! @tparam	M	A行列の高さ
	//! @tparam	N	A行列の幅
	//! @tparam	T	A行列のデータ型
	//! @tparam	MB	B行列の高さ
	//! @tparam	NB	B行列の幅
	//! @tparam	TB	B行列のデータ型
	//! @return	可制御性： true = 可制御、false = 不可制御
	template<size_t M, size_t N, typename T = double, size_t MB, size_t NB, typename TB = double>
	static constexpr bool IsCtrb(const ArcsMat<M,N,T>& Ac, const ArcsMat<MB,NB,TB>& Bc){
		static_assert(M == N,   "ArcsCtrl: Size Error");	// A行列は正方行列
		static_assert(MB == M , "ArcsCtrl: Size Error");	// サイズチェック

		ArcsMat<MB,NB,TB> AB;
		ArcsMat<MB,NB*M,TB> Co;
		for(size_t i = 0; i < M; ++i){
			AB = (Ac^i)*Bc;
			copymatrix(AB,1,MB,1,NB, Co,1,i*NB+1);	// 可制御性行列を生成
		}
		return M == rank(Co);	// ランクを計算して状態数と同一なら可観測
	}

//--------------------- ここから廃止予定

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

//--------------------- ここまで廃止予定

//! @brief 離散系状態空間モデルクラス
//! @tparam	N	次数
//! @tparam I	入力信号の数 (デフォルト値 = 1)
//! @tparam O	出力信号の数 (デフォルト値 = 1)
//! @tparam	T	データ型 (デフォルト値 = double)
template <size_t N, size_t I = 1, size_t O = 1, typename T = double>
class DiscStateSpace {
	public:
		//! @brief コンストラクタ(空コンストラクタ版)
		DiscStateSpace(void) noexcept
			: A(), B(), C(), D(), u(), x(), x_next(), y(), y_next()
		{
			PassedLog();
		}
		
		//! @brief コンストラクタ(離散系A,B,C行列設定版)
		//! @tparam	MA	A行列の高さ
		//! @tparam	NA	A行列の幅
		//! @tparam	TA	A行列のデータ型
		//! @tparam	MB	B行列の高さ
		//! @tparam	NB	B行列の幅
		//! @tparam	TB	B行列のデータ型
		//! @tparam	MC	C行列の高さ
		//! @tparam	NC	C行列の幅
		//! @tparam	TC	C行列のデータ型
		//! @param[in]	Ad	離散系A行列
		//! @param[in]	Bd	離散系B行列
		//! @param[in]	Cd	C行列
		template<
			size_t MA, size_t NA, typename TA = double, size_t MB, size_t NB, typename TB = double,
			size_t MC, size_t NC, typename TC = double
		>
		DiscStateSpace(const ArcsMat<MA,NA,TA>& Ad, const ArcsMat<MB,NB,TB>& Bd, const ArcsMat<MC,NC,TC>& Cd) noexcept
			: A(Ad), B(Bd), C(Cd), D(), u(), x(), x_next(), y(), y_next()
		{
			static_assert(MA == NA, "ArcsCtrl: Size Error");// 正方行列チェック
			static_assert(MA == N, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(MB == N, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(NB == I, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(MC == O, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(NC == N, "ArcsCtrl: Size Error");	// サイズチェック
			PassedLog();
		}
		
		//! @brief コンストラクタ(離散系A,B,C,D行列設定版)
		//! @tparam	MA	A行列の高さ
		//! @tparam	NA	A行列の幅
		//! @tparam	TA	A行列のデータ型
		//! @tparam	MB	B行列の高さ
		//! @tparam	NB	B行列の幅
		//! @tparam	TB	B行列のデータ型
		//! @tparam	MC	C行列の高さ
		//! @tparam	NC	C行列の幅
		//! @tparam	TC	C行列のデータ型
		//! @tparam	MD	D行列の高さ
		//! @tparam	ND	D行列の幅
		//! @tparam	TD	D行列のデータ型
		//! @param[in]	Ad	離散系A行列
		//! @param[in]	Bd	離散系B行列
		//! @param[in]	Cd	C行列
		//! @param[in]	Dd	D行列
		template<
			size_t MA, size_t NA, typename TA = double, size_t MB, size_t NB, typename TB = double,
			size_t MC, size_t NC, typename TC = double, size_t MD, size_t ND, typename TD = double
		>
		DiscStateSpace(const ArcsMat<MA,NA,TA>& Ad, const ArcsMat<MB,NB,TB>& Bd, const ArcsMat<MC,NC,TC>& Cd, const ArcsMat<MD,ND,TD>& Dd) noexcept
			: A(Ad), B(Bd), C(Cd), D(Dd), u(), x(), x_next(), y(), y_next()
		{
			static_assert(MA == NA, "ArcsCtrl: Size Error");// 正方行列チェック
			static_assert(MA == N, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(MB == N, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(NB == I, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(MC == O, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(NC == N, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(MD == O, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(ND == I, "ArcsCtrl: Size Error");	// サイズチェック
			PassedLog();
		}
		
		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		DiscStateSpace(DiscStateSpace&& r) noexcept
			: A(std::move(r.A)), B(std::move(r.B)), C(std::move(r.C)), D(std::move(r.D)),
			  u(std::move(r.u)), x(std::move(r.x)), x_next(std::move(r.x_next)), y(std::move(r.y)), y_next(std::move(r.y_next))
		{
			
		}

		//! @brief ムーブ代入演算子
		//! @param[in]	r	右辺値
		DiscStateSpace& operator=(DiscStateSpace&& r) noexcept {
			A = std::move(r.A);
			B = std::move(r.B);
			C = std::move(r.C);
			D = std::move(r.D);
			u = std::move(r.u);
			x = std::move(r.x);
			x_next = std::move(r.x_next);
			y = std::move(r.y);
			y_next = std::move(r.y_next);
			return *this;
		}

		//! @brief デストラクタ
		~DiscStateSpace() noexcept {
			PassedLog();
		}

		//! @brief 離散系状態空間モデルのA,B,C行列を設定する関数
		//! @tparam	MA	A行列の高さ
		//! @tparam	NA	A行列の幅
		//! @tparam	TA	A行列のデータ型
		//! @tparam	MB	B行列の高さ
		//! @tparam	NB	B行列の幅
		//! @tparam	TB	B行列のデータ型
		//! @tparam	MC	C行列の高さ
		//! @tparam	NC	C行列の幅
		//! @tparam	TC	C行列のデータ型
		//! @param[in]	Ad	離散系A行列
		//! @param[in]	Bd	離散系B行列
		//! @param[in]	Cd	C行列
		template<
			size_t MA, size_t NA, typename TA = double, size_t MB, size_t NB, typename TB = double,
			size_t MC, size_t NC, typename TC = double
		>
		void SetSystem(const ArcsMat<MA,NA,TA>& Ad, const ArcsMat<MB,NB,TB>& Bd, const ArcsMat<MC,NC,TC>& Cd) noexcept {
			static_assert(MA == NA, "ArcsCtrl: Size Error");// 正方行列チェック
			static_assert(MA == N, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(MB == N, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(NB == I, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(MC == O, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(NC == N, "ArcsCtrl: Size Error");	// サイズチェック

			A = Ad;
			B = Bd;
			C = Cd;
		}

		//! @brief 離散系状態空間モデルのA,B,C,D行列を設定する関数
		//! @tparam	MA	A行列の高さ
		//! @tparam	NA	A行列の幅
		//! @tparam	TA	A行列のデータ型
		//! @tparam	MB	B行列の高さ
		//! @tparam	NB	B行列の幅
		//! @tparam	TB	B行列のデータ型
		//! @tparam	MC	C行列の高さ
		//! @tparam	NC	C行列の幅
		//! @tparam	TC	C行列のデータ型
		//! @tparam	MD	D行列の高さ
		//! @tparam	ND	D行列の幅
		//! @tparam	TD	D行列のデータ型
		//! @param[in]	Ad	離散系A行列
		//! @param[in]	Bd	離散系B行列
		//! @param[in]	Cd	C行列
		//! @param[in]	Dd	D行列
		template<
			size_t MA, size_t NA, typename TA = double, size_t MB, size_t NB, typename TB = double,
			size_t MC, size_t NC, typename TC = double, size_t MD, size_t ND, typename TD = double
		>
		void SetSystem(const ArcsMat<MA,NA,TA>& Ad, const ArcsMat<MB,NB,TB>& Bd, const ArcsMat<MC,NC,TC>& Cd, const ArcsMat<MD,ND,TD>& Dd) noexcept {
			static_assert(MA == NA, "ArcsCtrl: Size Error");// 正方行列チェック
			static_assert(MA == N, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(MB == N, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(NB == I, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(MC == O, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(NC == N, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(MD == O, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(ND == I, "ArcsCtrl: Size Error");	// サイズチェック

			A = Ad;
			B = Bd;
			C = Cd;
			D = Dd;
		}

		//! @brief 状態空間モデルへの入力ベクトルを設定する関数
		//! @tparam	MU	入力ベクトルの高さ
		//! @tparam	NU	入力ベクトルの幅
		//! @tparam	TU	入力ベクトルのデータ型
		//! @param[in]	ui	入力ベクトル
		template<size_t MU, size_t NU, typename TU = double>
		void SetInput(const ArcsMat<MU,NU,TU>& ui) noexcept {
			static_assert(MU == I, "ArcsCtrl: Input Size Error");	// サイズチェック
			static_assert(NU == 1, "ArcsCtrl: Vector Error");		// サイズチェック
			u = ui;
		}
		
		//! @brief 状態空間モデルへの入力ベクトルの内の、１つの成分のみを選択して設定する関数
		//! @param[in]	i	入力ベクトルの要素番号(1～N)
		//! @param[in]	ui	入力値
		void SetInput(const size_t i, const T& ui) noexcept {
			arcs_assert(0 < i && i <= I);	// 範囲チェック
			u[i] = ui;
		}

		//! @brief 状態空間モデルへの入力ベクトルを設定する関数(個別1入力版)
		//! @param[in]	ui	入力
		void SetInput1(const T& u1) noexcept {
			static_assert(I == 1, "ArcsCtrl: Input Size Error");	// 1入力系かチェック
			u[1] = u1;
		}
		
		//! @brief 状態空間モデルへの入力ベクトルを設定する関数(個別2入力版)
		//! @param[in]	u1	入力1
		//! @param[in]	u2	入力2
		void SetInput2(const T& u1, const T& u2) noexcept {
			static_assert(I == 2, "ArcsCtrl: Input Size Error");	// 2入力系かチェック
			u[1] = u1;
			u[2] = u2;
		}
		
		//! @brief 状態空間モデルへの入力ベクトルを設定する関数(個別3入力版)
		//! @param[in]	u1	入力1
		//! @param[in]	u2	入力2
		//! @param[in]	u3	入力3
		void SetInput3(const T& u1, const T& u2, const T& u3) noexcept {
			static_assert(I == 3, "ArcsCtrl: Input Size Error");	// 3入力系かチェック
			u[1] = u1;
			u[2] = u2;
			u[3] = u3;
		}
		
		//! @brief 状態空間モデルの応答を計算して状態ベクトルを更新する関数
		void Update(void) noexcept {
			x_next = A*x + B*u;			// 状態方程式
			y = C*x + D*u;				// 出力方程式（教科書通りの正しい出力）
			y_next = C*x_next + D*u;	// 出力方程式（次の時刻の出力ベクトルを即時に得る場合）
			x = x_next;					// 状態ベクトルを更新
		}

		//! @brief 状態空間モデルの出力ベクトルを取得する関数(引数渡し版)
		//! @tparam	MY	出力ベクトルの高さ
		//! @tparam	NY	出力ベクトルの幅
		//! @tparam	TY	出力ベクトルのデータ型
		//! @param[out]	yo	出力ベクトル
		template<size_t MY, size_t NY, typename TY = double>
		void GetOutput(ArcsMat<MY,NY,TY>& yo) noexcept {
			static_assert(MY == O, "ArcsCtrl: Output Size Error");	// サイズチェック
			static_assert(NY == 1, "ArcsCtrl: Vector Error");		// サイズチェック
			yo = y;
		}
		
		//! @brief 状態空間モデルの出力ベクトルを取得する関数(戻り値返し版)
		//! @return	出力ベクトル
		ArcsMat<O,1> GetOutput(void) noexcept {
			return y;
		}
		
		//! @brief 状態空間モデルの出力ベクトルの内の、１つの成分のみを選択して取得する関数
		//! @param[in]	i	出力ベクトルの要素番号(1～N)
		//! @return	出力値
		T GetOutput(const size_t i) noexcept {
			arcs_assert(0 < i && i <= O);	// 範囲チェック
			return y[i];
		}

		//! @brief 状態空間モデルの出力を取得する関数(個別1出力版)
		//! @return	出力値
		T GetOutput1(void) noexcept {
			static_assert(O == 1, "ArcsCtrl: Output Size Error");	// 1出力系かチェック
			return y[1];
		}

		//! @brief 状態空間モデルの出力を取得する関数(個別2出力版)
		//! @return	出力値 (y1, y2) のタプル
		std::tuple<T,T> GetOutput2(void) noexcept {
			static_assert(O == 2, "ArcsCtrl: Output Size Error");	// 2出力系かチェック
			return {y[1], y[2]};
		}

		//! @brief 状態空間モデルの出力を取得する関数(個別3出力版)
		//! @return	出力値 (y1, y2, y3) のタプル
		std::tuple<T,T,T> GetOutput3(void) noexcept {
			static_assert(O == 3, "ArcsCtrl: Output Size Error");	// 3出力系かチェック
			return {y[1], y[2], y[3]};
		}

		//! @brief 状態空間モデルの次サンプルの出力ベクトルを即時に取得する関数(戻り値返し版)
		//! @return	出力ベクトル
		ArcsMat<O,1> GetNextOutput(void) noexcept {
			return y_next;
		}
		
		//! @brief 状態空間モデルの次サンプルの出力ベクトルの内の、１つの成分のみを選択して即時に取得する関数
		//! @param[in]	i	出力ベクトルの要素番号(1～N)
		//! @return	出力値
		T GetNextOutput(const size_t i) noexcept {
			arcs_assert(0 < i && i <= O);	// 範囲チェック
			return y_next[i];
		}

		//! @brief 状態空間モデルの次サンプルの出力を即時に取得する関数(個別1出力版)
		//! @return	出力値
		T GetNextOutput1(void) noexcept {
			static_assert(O == 1, "ArcsCtrl: Output Size Error");	// 1出力系かチェック
			return y_next[1];
		}

		//! @brief 状態空間モデルの次サンプルの出力を即時に取得する関数(個別2出力版)
		//! @return	出力値 (y1, y2) のタプル
		std::tuple<T,T> GetNextOutput2(void) noexcept {
			static_assert(O == 2, "ArcsCtrl: Output Size Error");	// 2出力系かチェック
			return {y_next[1], y_next[2]};
		}

		//! @brief 状態空間モデルの次サンプルの出力を即時に取得する関数(個別3出力版)
		//! @return	出力値
		std::tuple<T,T,T> GetNextOutput3(void) noexcept {
			static_assert(O == 3, "ArcsCtrl: Output Size Error");	// 3出力系かチェック
			return {y_next[1], y_next[2], y_next[3]};
		}
		
		//! @brief 状態ベクトルを任意の値にセットする関数
		//! @tparam	MX	入力される状態ベクトルの高さ
		//! @tparam	NX	入力される状態ベクトルの幅
		//! @tparam	TX	入力される状態ベクトルのデータ型
		//! @param[in]	xi	任意の値を持つ状態ベクトル
		template<size_t MX, size_t NX, typename TX = double>
		void SetStateVector(const ArcsMat<MX,NX,TX>& xi) noexcept {
			static_assert(MX == N, "ArcsCtrl: Size Error");		// サイズチェック
			static_assert(NX == 1, "ArcsCtrl: Vector Error");	// サイズチェック
			x = xi;
			x_next =xi;
		}

		//! @brief 状態ベクトルをクリアする関数
		void ClearStateVector(void) noexcept {
			x = ArcsMat<N,1>::zeros();
			x_next = ArcsMat<N,1>::zeros();
		}

	private:
		DiscStateSpace(const DiscStateSpace&) = delete;					//!< コピーコンストラクタ使用禁止
		const DiscStateSpace& operator=(const DiscStateSpace&) = delete;//!< コピー代入演算子使用禁止
		ArcsMat<N,N,T> A;		//!< 離散系A行列
		ArcsMat<N,I,T> B;		//!< 離散系B行列
		ArcsMat<O,N,T> C;		//!< C行列
		ArcsMat<O,I,T> D;		//!< D行列
		ArcsMat<I,1,T> u;		//!< 入力ベクトル
		ArcsMat<N,1,T> x;		//!< 状態ベクトル
		ArcsMat<N,1,T> x_next;	//!< 次の時刻の状態ベクトル
		ArcsMat<O,1,T> y;		//!< 出力ベクトル
		ArcsMat<O,1,T> y_next;	//!< 次の時刻の出力ベクトル
};

//! @brief 連続系状態空間モデルクラス
//! @tparam	N	次数
//! @tparam I	入力信号の数 (デフォルト値 = 1)
//! @tparam O	出力信号の数 (デフォルト値 = 1)
//! @tparam	T	データ型 (デフォルト値 = double)
template <size_t N, size_t I = 1, size_t O = 1, typename T = double>
class StateSpace {
	public:
		//! @brief コンストラクタ(空コンストラクタ版)
		StateSpace(void) noexcept
			: DiscSys(), A(), B(), C(), D(), Ts(), u(), x(), x_next(), y(), y_next()
		{
			PassedLog();
		}

		//! @brief コンストラクタ(連続系A,B,C行列設定版)
		//! @tparam	MA	A行列の高さ
		//! @tparam	NA	A行列の幅
		//! @tparam	TA	A行列のデータ型
		//! @tparam	MB	B行列の高さ
		//! @tparam	NB	B行列の幅
		//! @tparam	TB	B行列のデータ型
		//! @tparam	MC	C行列の高さ
		//! @tparam	NC	C行列の幅
		//! @tparam	TC	C行列のデータ型
		//! @param[in]	Ac	連続系A行列
		//! @param[in]	Bc	連続系B行列
		//! @param[in]	Cc	C行列
		//! @param[in]	Tsmpl	[s] サンプリング周期
		template<
			size_t MA, size_t NA, typename TA = double, size_t MB, size_t NB, typename TB = double,
			size_t MC, size_t NC, typename TC = double
		>
		StateSpace(const ArcsMat<MA,NA,TA>& Ac, const ArcsMat<MB,NB,TB>& Bc, const ArcsMat<MC,NC,TC>& Cc, const T& Tsmpl) noexcept
			: DiscSys(), A(Ac), B(Bc), C(Cc), D(), Ts(Tsmpl), u(), x(), x_next(), y(), y_next()
		{
			static_assert(MA == NA, "ArcsCtrl: Size Error");// 正方行列チェック
			static_assert(MA == N, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(MB == N, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(NB == I, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(MC == O, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(NC == N, "ArcsCtrl: Size Error");	// サイズチェック
			arcs_assert(0 < Tsmpl);	// サンプリング周期は非負且つ非零
			ConvToDiscSystem();		// 離散系システムに変換
			PassedLog();
		}
		
		//! @brief コンストラクタ(離散系A,B,C,D行列設定版)
		//! @tparam	MA	A行列の高さ
		//! @tparam	NA	A行列の幅
		//! @tparam	TA	A行列のデータ型
		//! @tparam	MB	B行列の高さ
		//! @tparam	NB	B行列の幅
		//! @tparam	TB	B行列のデータ型
		//! @tparam	MC	C行列の高さ
		//! @tparam	NC	C行列の幅
		//! @tparam	TC	C行列のデータ型
		//! @tparam	MD	D行列の高さ
		//! @tparam	ND	D行列の幅
		//! @tparam	TD	D行列のデータ型
		//! @param[in]	Ad	離散系A行列
		//! @param[in]	Bd	離散系B行列
		//! @param[in]	Cd	C行列
		//! @param[in]	Dd	D行列
		//! @param[in]	Tsmpl	[s] サンプリング周期
		template<
			size_t MA, size_t NA, typename TA = double, size_t MB, size_t NB, typename TB = double,
			size_t MC, size_t NC, typename TC = double, size_t MD, size_t ND, typename TD = double
		>
		StateSpace(const ArcsMat<MA,NA,TA>& Ad, const ArcsMat<MB,NB,TB>& Bd, const ArcsMat<MC,NC,TC>& Cd, const ArcsMat<MD,ND,TD>& Dd, const T& Tsmpl) noexcept
			: DiscSys(), A(Ad), B(Bd), C(Cd), D(Dd), Ts(Tsmpl), u(), x(), x_next(), y(), y_next()
		{
			static_assert(MA == NA, "ArcsCtrl: Size Error");// 正方行列チェック
			static_assert(MA == N, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(MB == N, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(NB == I, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(MC == O, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(NC == N, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(MD == O, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(ND == I, "ArcsCtrl: Size Error");	// サイズチェック
			arcs_assert(0 < Tsmpl);	// サンプリング周期は非負且つ非零
			ConvToDiscSystem();		// 離散系システムに変換
			PassedLog();
		}

		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		StateSpace(StateSpace&& r) noexcept
			: DiscSys(std::move(r.DiscSys)), A(std::move(r.A)), B(std::move(r.B)), C(std::move(r.C)), D(std::move(r.D)), Ts(std::move(r.Ts)), 
			  u(std::move(r.u)), x(std::move(r.x)), x_next(std::move(r.x_next)), y(std::move(r.y)), y_next(std::move(r.y_next))
		{
			
		}

		//! @brief ムーブ代入演算子
		//! @param[in]	r	右辺値
		StateSpace& operator=(StateSpace&& r) noexcept {
			DiscSys = std::move(r.DiscSys);
			A = std::move(r.A);
			B = std::move(r.B);
			C = std::move(r.C);
			D = std::move(r.D);
			Ts = std::move(r.Ts);
			u = std::move(r.u);
			x = std::move(r.x);
			x_next = std::move(r.x_next);
			y = std::move(r.y);
			y_next = std::move(r.y_next);
			return *this;
		}
		
		//! @brief デストラクタ
		~StateSpace() noexcept {
			PassedLog();
		}
		
		//! @brief 連続系状態空間モデルのA,B,C行列を設定する関数
		//! @tparam	MA	A行列の高さ
		//! @tparam	NA	A行列の幅
		//! @tparam	TA	A行列のデータ型
		//! @tparam	MB	B行列の高さ
		//! @tparam	NB	B行列の幅
		//! @tparam	TB	B行列のデータ型
		//! @tparam	MC	C行列の高さ
		//! @tparam	NC	C行列の幅
		//! @tparam	TC	C行列のデータ型
		//! @param[in]	Ac	連続系A行列
		//! @param[in]	Bc	連続系B行列
		//! @param[in]	Cc	C行列
		//! @param[in]	Tsmpl	[s] サンプリング周期
		template<
			size_t MA, size_t NA, typename TA = double, size_t MB, size_t NB, typename TB = double,
			size_t MC, size_t NC, typename TC = double
		>
		void SetSystem(const ArcsMat<MA,NA,TA>& Ac, const ArcsMat<MB,NB,TB>& Bc, const ArcsMat<MC,NC,TC>& Cc, const T& Tsmpl) noexcept {
			static_assert(MA == NA, "ArcsCtrl: Size Error");// 正方行列チェック
			static_assert(MA == N, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(MB == N, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(NB == I, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(MC == O, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(NC == N, "ArcsCtrl: Size Error");	// サイズチェック

			A = Ac;
			B = Bc;
			C = Cc;
			Ts = Tsmpl;
			ConvToDiscSystem();		// 離散系システムに変換
		}

		//! @brief 連続系状態空間モデルのA,B,C,D行列を設定する関数
		//! @tparam	MA	A行列の高さ
		//! @tparam	NA	A行列の幅
		//! @tparam	TA	A行列のデータ型
		//! @tparam	MB	B行列の高さ
		//! @tparam	NB	B行列の幅
		//! @tparam	TB	B行列のデータ型
		//! @tparam	MC	C行列の高さ
		//! @tparam	NC	C行列の幅
		//! @tparam	TC	C行列のデータ型
		//! @tparam	MD	D行列の高さ
		//! @tparam	ND	D行列の幅
		//! @tparam	TD	D行列のデータ型
		//! @param[in]	Ac	連続系A行列
		//! @param[in]	Bc	連続系B行列
		//! @param[in]	Cc	C行列
		//! @param[in]	Dc	D行列
		//! @param[in]	Tsmpl	[s] サンプリング周期
		template<
			size_t MA, size_t NA, typename TA = double, size_t MB, size_t NB, typename TB = double,
			size_t MC, size_t NC, typename TC = double, size_t MD, size_t ND, typename TD = double
		>
		void SetSystem(const ArcsMat<MA,NA,TA>& Ac, const ArcsMat<MB,NB,TB>& Bc, const ArcsMat<MC,NC,TC>& Cc, const ArcsMat<MD,ND,TD>& Dc, const T& Tsmpl) noexcept {
			static_assert(MA == NA, "ArcsCtrl: Size Error");// 正方行列チェック
			static_assert(MA == N, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(MB == N, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(NB == I, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(MC == O, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(NC == N, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(MD == O, "ArcsCtrl: Size Error");	// サイズチェック
			static_assert(ND == I, "ArcsCtrl: Size Error");	// サイズチェック

			A = Ac;
			B = Bc;
			C = Cc;
			D = Dc;
			Ts = Tsmpl;
			ConvToDiscSystem();		// 離散系システムに変換
		}

		//! @brief 状態空間モデルへの入力ベクトルを設定する関数
		//! @tparam	MU	入力ベクトルの高さ
		//! @tparam	NU	入力ベクトルの幅
		//! @tparam	TU	入力ベクトルのデータ型
		//! @param[in]	ui	入力ベクトル
		template<size_t MU, size_t NU, typename TU = double>
		void SetInput(const ArcsMat<MU,NU,TU>& ui) noexcept {
			static_assert(MU == I, "ArcsCtrl: Input Size Error");	// サイズチェック
			static_assert(NU == 1, "ArcsCtrl: Vector Error");		// サイズチェック
			DiscSys.SetInput(ui);
		}
		
		//! @brief 状態空間モデルへの入力ベクトルの内の、１つの成分のみを選択して設定する関数
		//! @param[in]	i	入力ベクトルの要素番号(1～N)
		//! @param[in]	ui	入力値
		void SetInput(const size_t i, const T& ui) noexcept {
			arcs_assert(0 < i && i <= I);	// 範囲チェック
			DiscSys.SetInput(i, ui);
		}

		//! @brief 状態空間モデルへの入力ベクトルを設定する関数(個別1入力版)
		//! @param[in]	ui	入力
		void SetInput1(const T& u1) noexcept {
			static_assert(I == 1, "ArcsCtrl: Input Size Error");	// 1入力系かチェック
			DiscSys.SetInput1(u1);
		}
		
		//! @brief 状態空間モデルへの入力ベクトルを設定する関数(個別2入力版)
		//! @param[in]	u1	入力1
		//! @param[in]	u2	入力2
		void SetInput2(const T& u1, const T& u2) noexcept {
			static_assert(I == 2, "ArcsCtrl: Input Size Error");	// 2入力系かチェック
			DiscSys.SetInput2(u1, u2);
		}
		
		//! @brief 状態空間モデルへの入力ベクトルを設定する関数(個別3入力版)
		//! @param[in]	u1	入力1
		//! @param[in]	u2	入力2
		//! @param[in]	u3	入力3
		void SetInput3(const T& u1, const T& u2, const T& u3) noexcept {
			static_assert(I == 3, "ArcsCtrl: Input Size Error");	// 3入力系かチェック
			DiscSys.SetInput3(u1, u2, u3);
		}

		//! @brief 状態空間モデルの応答を計算して状態ベクトルを更新する関数
		void Update(void) noexcept {
			DiscSys.Update();
		}

		//! @brief 状態空間モデルの出力ベクトルを取得する関数(引数渡し版)
		//! @tparam	MY	出力ベクトルの高さ
		//! @tparam	NY	出力ベクトルの幅
		//! @tparam	TY	出力ベクトルのデータ型
		//! @param[out]	yo	出力ベクトル
		template<size_t MY, size_t NY, typename TY = double>
		void GetOutput(ArcsMat<MY,NY,TY>& yo) noexcept {
			static_assert(MY == O, "ArcsCtrl: Output Size Error");	// サイズチェック
			static_assert(NY == 1, "ArcsCtrl: Vector Error");		// サイズチェック
			DiscSys.GetOutput(yo);
		}
		
		//! @brief 状態空間モデルの出力ベクトルを取得する関数(戻り値返し版)
		//! @return	出力ベクトル
		ArcsMat<O,1> GetOutput(void) noexcept {
			return DiscSys.GetOutput();
		}
		
		//! @brief 状態空間モデルの出力ベクトルの内の、１つの成分のみを選択して取得する関数
		//! @param[in]	i	出力ベクトルの要素番号(1～N)
		//! @return	出力値
		T GetOutput(const size_t i) noexcept {
			arcs_assert(0 < i && i <= O);	// 範囲チェック
			return DiscSys.GetOutput(i);
		}

		//! @brief 状態空間モデルの出力を取得する関数(個別1出力版)
		//! @return	出力値
		T GetOutput1(void) noexcept {
			static_assert(O == 1, "ArcsCtrl: Output Size Error");	// 1出力系かチェック
			return DiscSys.GetOutput1();
		}

		//! @brief 状態空間モデルの出力を取得する関数(個別2出力版)
		//! @return	出力値 (y1, y2) のタプル
		std::tuple<T,T> GetOutput2(void) noexcept {
			static_assert(O == 2, "ArcsCtrl: Output Size Error");	// 2出力系かチェック
			return DiscSys.GetOutput2();
		}

		//! @brief 状態空間モデルの出力を取得する関数(個別3出力版)
		//! @return	出力値 (y1, y2, y3) のタプル
		std::tuple<T,T,T> GetOutput3(void) noexcept {
			static_assert(O == 3, "ArcsCtrl: Output Size Error");	// 3出力系かチェック
			return DiscSys.GetOutput3();
		}

	private:
		StateSpace(const StateSpace&) = delete;					//!< コピーコンストラクタ使用禁止
		const StateSpace& operator=(const StateSpace&) = delete;//!< コピー代入演算子使用禁止

		//!< 連続系状態空間モデルを離散系へ変換する関数
		void ConvToDiscSystem(void){
			// 連続系システムを平衡実現
			ArcsMat<N,N,T> Ah;
			ArcsMat<N,I,T> Bh;
			ArcsMat<O,N,T> Ch;
			if(IsStable(A) == true){
				// 安定系であれば、平衡実現
				BalanceReal(A, B, C, Ah, Bh, Ch);
			}else{
				// 不安定系であれば、平衡実現をしないでそのまま
				Ah = A;
				Bh = B;
				Ch = C;
			}

			// 連続系システムを離散化
			ArcsMat<N,N,T> Ad;
			ArcsMat<N,I,T> Bd;
			Discretize(Ah, Bh, Ad, Bd, Ts);	// 離散化

			// 離散系システムに設定
			DiscSys.SetSystem(Ad, Bd, Ch, D);
		}

		DiscStateSpace<N, I, O> DiscSys;	//!< 離散系状態空間モデル
		ArcsMat<N,N,T> A;		//!< 連続系A行列
		ArcsMat<N,I,T> B;		//!< 連続系B行列
		ArcsMat<O,N,T> C;		//!< C行列
		ArcsMat<O,I,T> D;		//!< D行列
		T Ts;					//!< [s] サンプリング周期
		ArcsMat<I,1,T> u;		//!< 入力ベクトル
		ArcsMat<N,1,T> x;		//!< 状態ベクトル
		ArcsMat<N,1,T> x_next;	//!< 次の時刻の状態ベクトル
		ArcsMat<O,1,T> y;		//!< 出力ベクトル
		ArcsMat<O,1,T> y_next;	//!< 次の時刻の出力ベクトル
};


}
}

#endif


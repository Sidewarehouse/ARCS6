//! @file ArcsControl.hh
//! @brief ARCS-Control 制御理論クラス
//!
//! 制御理論に関係する様々なアルゴリズムを詰め合わせた静的クラス
//!
//! @date 2024/07/23
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

namespace ARCS {	// ARCS名前空間
//! @brief ARCS-Control 制御理論クラス
class ArcsControl {
	public:
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

		//---------------------

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


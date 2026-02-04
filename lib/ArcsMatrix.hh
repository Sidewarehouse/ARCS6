//! @file ArcsMatrix.hh
//! @brief ARCS-Matrix 行列演算クラス
//!
//! ArcsMat: 行列に関係する様々な演算を実行するクラス
//! ・スタック領域なのでかなり高速の行列演算が可能
//!
//! ArcsVarMat: 可変サイズの行列に関係する簡単な演算を実行するクラス
//! ・ArcsMatのヒープ領域版
//! ・ヒープ領域を使うためArcsMatに比べて低速
//!
//! @date 2026/01/02
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2026 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.
//
// ・各関数における計算結果はMATLAB/Maximaと比較して合っていることを確認済み。
// ・下記に関して、MATLABと異なる結果を出力する場合がある：
// 　- 行列分解の符号関係
//　 - 複素数のときの行列分解
// 　- 劣決定のときの線形方程式の解
// ・旧来のMatrixクラスのバグ・不具合・不満を一掃。

#ifndef ARCSMATRIX
#define ARCSMATRIX

#include <string>
#include <tuple>
#include <cmath>
#include <cassert>
#include <array>
#include <vector>
#include <complex>
#include <cstdint>
#include <algorithm>

// ARCS組込み用マクロ
#ifdef ARCS_IN
	// ARCSに組み込まれる場合
	#include "ARCSassert.hh"
#else
	// ARCSに組み込まれない場合
	#define arcs_assert(a) (assert(a))
#endif

// 表示用マクロ
#define dispsize(a)  (ArcsMatrix::dispsize_macro((a),#a))	//!< 行列サイズ表示マクロ
#define dispf(a,b) (ArcsMatrix::dispf_macro((a),b,#a))		//!< 行列要素表示マクロ(フォーマット指定あり版)
#define disp(a) (ArcsMatrix::disp_macro((a),#a))			//!< 行列要素表示マクロ(フォーマット指定なし版)

using namespace std::literals::complex_literals;	// 虚数単位リテラル「i」の使用

// ARCS名前空間
namespace ARCS {

// ArcsMatrix設定定義
namespace ArcsMatrix {
	//! @brief 行列の状態の定義
	enum class MatStatus {
		AMT_NA,		//!< 状態定義該当なし
		AMT_LU_ODD,	//!< LU分解したときに並べ替えが奇数回発生
		AMT_LU_EVEN	//!< LU分解したときに並べ替えが偶数回発生
	};

	//! @brief ノルム計算方法の定義
	enum class NormType {
		AMT_L2,		//!< ユークリッドノルム(2-ノルム)
		AMT_L1,		//!< 絶対値ノルム(1-ノルム)
		AMT_LINF	//!< 無限大ノルム(最大値ノルム)
	};

	//! @brief ソート方法の定義
	enum class SortType {
		AMT_ASCENT,	//!< 昇順
		AMT_DESCENT	//!< 降順
	};
}

// ArcsMatrixメタ関数定義
namespace ArcsMatrix {
	// 整数・実数型チェック用メタ関数
	template<typename TT> struct IsIntFloatV {
		static constexpr bool value = std::is_integral<TT>::value | std::is_floating_point<TT>::value;
	};
	template<typename TT> inline constexpr bool IsIntFloat = IsIntFloatV<TT>::value;

	// 複素数型チェック用メタ関数
	template<typename TT> struct IsComplexV : std::false_type {};
	template<> struct IsComplexV<std::complex<double>> : std::true_type {};
	template<> struct IsComplexV<std::complex<float>> : std::true_type {};
	template<> struct IsComplexV<std::complex<long>> : std::true_type {};
	template<> struct IsComplexV<std::complex<int>> : std::true_type {};
	template<typename TT> inline constexpr bool IsComplex = IsComplexV<TT>::value;
	
	// 対応可能型チェック用メタ関数
	template<typename TT> inline constexpr bool IsApplicable = IsIntFloatV<TT>::value | IsComplexV<TT>::value;
}

// ArcsMatrix静的定数と関数（行列処理に必要な定数と関数）
namespace ArcsMatrix {
	static constexpr double EPSILON = 1e-14;		//!< 零とみなす閾値(実数版)
	static constexpr std::complex<double> EPSLCOMP = std::complex(1e-14, 1e-14);	//!< 零とみなす閾値(複素数版)

	//! @brief 符号関数
	//! @tparam	T	データ型
	//! @param[in]	u	入力
	//! @return	符号結果
	template<typename T>
	static constexpr T sgn(T u){
		T ret = 1;
		if constexpr(ArcsMatrix::IsComplex<T>){
			// 複素数型の場合
			if(std::abs(u) < EPSILON){
				ret = std::complex(0.0, 0.0);	// ゼロ割回避
			}else{
				ret = u/std::abs(u);			// 複素数に拡張された符号関数
			}
		}else{
			// 実数型の場合
			ret = std::copysign(1, u);
		}
		return ret;
	}

	//! @brief 二項係数 nCk を計算する関数
	//! @param[in]	n	n個から、
	//! @param[in]	k	k個取り出す組み合わせ
	//! @return	二項係数の結果
	static constexpr size_t nChoosek(const size_t n, const size_t k){
		// Knuth先生の方法
		size_t ret = 0;
		if( k == 0 || k == n){
			ret = 1;
		}else{
			ret = nChoosek(n - 1, k - 1)*n/k;	// 再帰
		}
		return ret;
	}
}

//! @brief ARCS-Matrix 行列演算クラス
//! @tparam M	行列の高さ
//! @tparam	N	行列の幅
//! @tparam T	データ型(デフォルトはdouble型)
template <size_t M, size_t N, typename T = double>
class ArcsMat {
	public:
		//! @brief コンストラクタ
		constexpr ArcsMat(void) noexcept
			: Nindex(0), Mindex(0), Status(ArcsMatrix::MatStatus::AMT_NA), Data({0})
		{
			static_assert(N != 0, "ArcsMat: Size Zero Error");	// サイズゼロの行列は禁止
			static_assert(M != 0, "ArcsMat: Size Zero Error");	// サイズゼロの行列は禁止
			static_assert(ArcsMatrix::IsApplicable<T>, "ArcsMat: Type Error");	// 対応可能型チェック

			FillAll(0);				// すべての要素を零で初期化
		}
		
		//! @brief コンストラクタ(任意初期値版)
		//! @tparam	R	演算子右側の要素の型
		//! @param[in]	InitValue	行列要素の初期値
		template<typename R>
		constexpr explicit ArcsMat(const R InitValue) noexcept
			: Nindex(0), Mindex(0), Status(ArcsMatrix::MatStatus::AMT_NA), Data({0})
		{
			static_assert(N != 0, "ArcsMat: Size Zero Error");	// サイズゼロの行列は禁止
			static_assert(M != 0, "ArcsMat: Size Zero Error");	// サイズゼロの行列は禁止
			static_assert(ArcsMatrix::IsApplicable<T>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック

			FillAll(static_cast<T>(InitValue));			// すべての要素を指定した値で初期化
		}
		
		//! @brief コンストラクタ(初期化リスト版)
		//! @tparam	R	演算子右側の要素の型
		//! @param[in]	InitList	初期化リスト
		//template<typename R = T>
		constexpr ArcsMat(const std::initializer_list<T> InitList) noexcept
			: Nindex(0), Mindex(0), Status(ArcsMatrix::MatStatus::AMT_NA), Data({0})
		{
			static_assert(N != 0, "ArcsMat: Size Zero Error");	// サイズゼロの行列は禁止
			static_assert(M != 0, "ArcsMat: Size Zero Error");	// サイズゼロの行列は禁止
			static_assert(ArcsMatrix::IsApplicable<T>, "ArcsMat: Type Error");	// 対応可能型チェック
			//static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック

			const T* ListVal = InitList.begin();		// 初期化リストの最初のポインタ位置
			size_t Ni = 0;				// 横方向カウンタ
			size_t Mi = 0;				// 縦方向カウンタ
			for(size_t i = 0; i < InitList.size(); ++i){
				// 初期化リストを順番に読み込んでいく
				arcs_assert(Ni < N);	// 横方向カウンタが行の長さ以内かチェック
				arcs_assert(Mi < M);	// 縦方向カウンタが列の高さ以内かチェック
				Data[Ni][Mi] = static_cast<T>(ListVal[i]);	// キャストしてから行列の要素を埋める
				Ni++;					// 横方向カウンタをカウントアップ
				if(Ni == N){			// 横方向カウンタが最後まで行き着いたら，
					Ni = 0;				// 横方向カウンタを零に戻して，
					Mi++;				// その代わりに，縦方向カウンタをカウントアップ
				}
			}
		}
		
		//! @brief コピーコンストラクタ
		//! @param[in]	right	演算子右側
		constexpr ArcsMat(const ArcsMat<M,N,T>& right) noexcept
			: Nindex(0), Mindex(0), Status(right.GetStatus()), Data(right.GetData())
		{
			static_assert(ArcsMatrix::IsApplicable<T>, "ArcsMat: Type Error");	// 対応可能型チェック
			// メンバを取り込む以外の処理は無し
		}
		
		//! @brief コピーコンストラクタ(サイズもしくは型が違う行列の場合の定義)
		//! @tparam	P, Q, R	演算子右側の行列の高さ, 幅, 要素の型
		//! @param[in]	right	演算子右側
		template<size_t P, size_t Q, typename R = double>
		constexpr ArcsMat(const ArcsMat<P,Q,R>& right) noexcept
			: Nindex(0), Mindex(0), Status(right.GetStatus()), Data()
		{
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック

			// 型の種類によってキャスト処理を変更
			if constexpr(ArcsMatrix::IsIntFloat<R> && ArcsMatrix::IsIntFloat<T>){
				// 「整数・浮動小数点型 → 整数・浮動小数点型」のキャスト処理
				Data = static_cast<T>(right.GetData());
			}else if constexpr(ArcsMatrix::IsIntFloat<R> && ArcsMatrix::IsComplex<T>){
				// 「整数・浮動小数点型 → 複素数型」のキャスト処理
				const std::array<std::array<R, M>, N>& RightData = right.ReadOnlyRef();
				for(size_t i = 0; i < N; ++i){
					for(size_t j = 0; j < M; ++j) Data[i][j].real(RightData[i][j]);
				}
			}else if constexpr(ArcsMatrix::IsComplex<R> && ArcsMatrix::IsComplex<T>){
				// 「複素数型 → 複素数型」のキャスト処理
				Data = static_cast<T>(right.GetData());
			}else{
				// キャスト不能の場合の処理
				arcs_assert(false);	// 変換不能、もしここに来たら問答無用でAssertion Failed
			}
		}
		
		//! @brief 行列コピー代入演算子(サイズと型が同じ同士の行列の場合)
		//! @param[in] right 演算子の右側
		//! @return 結果
		constexpr ArcsMat<M,N,T>& operator=(const ArcsMat<M,N,T>& right) noexcept {
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) (*this)(j,i) = right(j,i);
				// 上記は次とほぼ同等→ this->Data[i][j] = right.Data[i][j];
			}
			return (*this);
		}
		
		//! @brief 行列コピー代入演算子(サイズと型が違う行列の場合, エラー検出用の定義)
		//! @tparam	P, Q, R	演算子右側の行列の高さ, 幅, 要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<size_t P, size_t Q, typename R = double>
		constexpr ArcsMat<M,N,T>& operator=(const ArcsMat<P,Q,R>& right) noexcept {
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			arcs_assert(false);	// もしここに来たら問答無用でAssertion Failed
			return (*this);
		}
		
		//! @brief ムーブコンストラクタ
		//! @param[in]	right	演算子右側
		constexpr ArcsMat(ArcsMat<M,N,T>&& right) noexcept
			: Nindex(), Mindex(0), Status(right.GetStatus()), Data(right.GetData())
		{
			// メンバを取り込む以外の処理は無し
		}
		
		//! @brief ムーブコンストラクタ(サイズと型が違う行列の場合, エラー検出用の定義)
		//! @tparam	P, Q, R	演算子右側の行列の高さ, 幅, 要素の型
		//! @param[in]	right	演算子右側
		template<size_t P, size_t Q, typename R = double>
		constexpr ArcsMat(ArcsMat<P,Q,R>&& right) noexcept
			: Nindex(0), Mindex(0), Status(right.GetStatus()), Data(right.GetData())
		{
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			arcs_assert(false);	// もしここに来たら問答無用でAssertion Failed
		}

		//! @brief 行列ムーブ代入演算子(サイズと型が同じ同士の行列の場合) [未実装]
		//! @param[in]	r	演算子右側
		//constexpr ArcsMat<M,N,T>& operator=(ArcsMat<M,N,T>&& right) noexcept {
		//	return *this;
		//}

		//! @brief 行列ムーブ代入演算子(サイズと型が違う行列の場合, エラー検出用の定義) [未実装]
		//! @tparam	P, Q, R	演算子右側の行列の高さ, 幅, 要素の型
		//! @param[in]	r	演算子右側
		//template<size_t P, size_t Q, typename R = double>
		//constexpr ArcsMat<M,N,T>& operator=(ArcsMat<P,Q,R>&& right) noexcept {
		//	return *this;
		//}
		
		//! @brief 縦ベクトル添字演算子(縦ベクトルのm番目の要素の値を返す。x = A(m,1)と同じ意味)
		//!        備考：ArcsMatは縦ベクトル優先なので、横ベクトル添字演算子は無い。
		//! @param[in]	m	縦方向の要素番号( "1" 始まり)
		//! @return	要素の値
		constexpr T operator[](const size_t m) const{
			static_assert(N == 1, "ArcsMat: Vector Error");		// 縦ベクトルチェック
			return Data[0][m - 1];
		}
		
		//! @brief 縦ベクトル添字演算子(縦ベクトルのm番目の要素に値を設定する。A(m,1) = xと同じ意味)
		//!        備考：ArcsMatは縦ベクトル優先なので、横ベクトル添字演算子は無い。
		//! @param[in]	m	縦方向の要素番号( "1" 始まり)
		//! @return	設定後の縦ベクトル
		constexpr T& operator[](const size_t m){
			static_assert(N == 1, "ArcsMat: Vector Error");		// 縦ベクトルチェック
			return Data[0][m - 1];
		}
		
		//! @brief 行列括弧演算子(行列の(m,n)要素の値を返す。サイズチェック無し版)
		//! @param[in]	m	m行目(縦方向の位置)
		//! @param[in]	n	n列目(横方向の位置)
		//! @return	要素の値
		constexpr T operator()(const size_t m, const size_t n) const{
			return Data[n - 1][m - 1];
		}
		
		//! @brief 行列括弧演算子(行列の(m,n)要素に値を設定する。サイズチェック無し版)
		//! @param[in]	m	m行目(縦方向の位置)
		//! @param[in]	n	n列目(横方向の位置)
		//! @return	設定後の行列
		constexpr T& operator()(const size_t m, const size_t n){
			return Data[n - 1][m - 1];
		}
		
		//! @brief 行列括弧演算子(行列の(m,n)要素の値を返す。サイズチェック可能版)
		//! @param[in]	m	m行目(縦方向の位置)
		//! @param[in]	n	n列目(横方向の位置)
		//! @param[in]	chk	サイズチェックフラグ true = サイズチェックする, false = サイズチェックしない
		//! @return	要素の値
		constexpr T operator()(const size_t m, const size_t n, const bool chk) const{
			if(chk == true){
				arcs_assert(0 < m && m <= M);
				arcs_assert(0 < n && n <= N);
			}
			return Data[n - 1][m - 1];
		}
		
		//! @brief 行列括弧演算子(行列の(m,n)要素に値を設定する。サイズチェック可能版)
		//! @param[in]	m	m行目(縦方向の位置)
		//! @param[in]	n	n列目(横方向の位置)
		//! @param[in]	chk	サイズチェックフラグ true = サイズチェックする, false = サイズチェックしない
		//! @return	設定後の行列
		constexpr T& operator()(const size_t m, const size_t n, const bool chk){
			if(chk == true){
				arcs_assert(0 < m && m <= M);
				arcs_assert(0 < n && n <= N);
			}
			return Data[n - 1][m - 1];
		}
		
		//! @brief 単項プラス演算子
		//! @return 結果
		constexpr ArcsMat<M,N,T> operator+(void) const{
			ArcsMat<M,N,T> ret;
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) ret(j,i) = (*this)(j,i);
				// 上記は次とほぼ同等→ ret.Data[i][j] = Data[i][j];
			}
			return ret;
		}
		
		//! @brief 単項マイナス演算子
		//! @return 結果
		constexpr ArcsMat<M,N,T> operator-(void) const{
			ArcsMat<M,N,T> ret;
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) ret(j,i) = -(*this)(j,i);
				// 上記は次とほぼ同等→ ret.Data[i][j] = -Data[i][j];
			}
			return ret;
		}
		
		//! @brief 行列加算演算子(行列＝行列＋行列の場合)
		//! @tparam	P, Q, R	演算子右側の行列の高さ, 幅, 要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<size_t P, size_t Q, typename R = double>
		constexpr ArcsMat<M,N,T> operator+(const ArcsMat<P,Q,R>& right) const{
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			ArcsMat<M,N,T> ret;
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) ret(j,i) = (*this)(j,i) + static_cast<T>( right(j,i) );
				// 上記は次とほぼ同等→ ret.Data[i][j] = Data[i][j] + static_cast<T>(right.Data[i][j]);
			}
			return ret;
		}
		
		//! @brief 行列加算演算子(行列＝行列＋スカラーの場合)
		//! @tparam	R	演算子右側の要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<typename R>
		constexpr ArcsMat<M,N,T> operator+(const R& right) const{
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			ArcsMat<M,N,T> ret;
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) ret(j,i) = (*this)(j,i) + static_cast<T>(right);
				// 上記は次とほぼ同等→ ret.Data[i][j] = Data[i][j] + right;
			}
			return ret;
		}
		
		//! @brief 行列減算演算子(行列＝行列－行列の場合)
		//! @tparam	P, Q, R	演算子右側の行列の高さ, 幅, 要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<size_t P, size_t Q, typename R = double>
		constexpr ArcsMat<M,N,T> operator-(const ArcsMat<P,Q,R>& right) const{
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			ArcsMat<M,N,T> ret;
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) ret(j,i) = (*this)(j,i) - static_cast<T>( right(j,i) );
				// 上記は次とほぼ同等→ ret.Data[i][j] = Data[i][j] - right.Data[i][j];
			}
			return ret;
		}
		
		//! @brief 行列減算演算子(行列＝行列－スカラーの場合)
		//! @tparam	R	演算子右側の要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<typename R>
		constexpr ArcsMat<M,N,T> operator-(const R& right) const{
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			ArcsMat<M,N,T> ret;
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) ret(j,i) = (*this)(j,i) - static_cast<T>(right);
				// 上記は次とほぼ同等→ ret.Data[i][j] = Data[i][j] - right;
			}
			return ret;
		}
		
		//! @brief 行列乗算演算子(行列＝行列＊行列の場合)
		//! @tparam	P, Q, R	演算子右側の行列の高さ, 幅, 要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<size_t P, size_t Q, typename R = double>
		constexpr ArcsMat<M,Q,T> operator*(const ArcsMat<P,Q,R>& right) const{
			static_assert(N == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			ArcsMat<M,Q,T> ret;
			for(size_t k = 1; k <= Q; ++k){
				for(size_t i = 1; i <= N; ++i){
					for(size_t j = 1; j <= M; ++j) ret(j,k) += (*this)(j,i)*static_cast<T>( right(i,k) );
					// 上記は次とほぼ同等→ ret.Data[k][j] += Data[i][j]*right.Data[k][i];
				}
			}
			return ret;
		}
		
		//! @brief 行列乗算演算子(行列＝行列＊スカラーの場合)
		//! @tparam	R	演算子右側の要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<typename R>
		constexpr ArcsMat<M,N,T> operator*(const R& right) const{
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			ArcsMat<M,N,T> ret;
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) ret(j,i) = (*this)(j,i)*static_cast<T>(right);
				// 上記は次とほぼ同等→ ret.Data[i][j] = Data[i][j]*right;
			}
			return ret;
		}
		
		//! @brief 行列除算演算子(行列＝行列／行列の場合)
		//! @tparam	P, Q, R	演算子右側の行列の高さ, 幅, 要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<size_t P, size_t Q, typename R = double>
		constexpr void operator/(const ArcsMat<P,Q,R>& right) const{
			arcs_assert(false);		// この演算子は使用禁止
		}

		//! @brief 行列除算演算子(行列＝行列／スカラーの場合)
		//! @tparam	R	演算子右側の要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<typename R>
		constexpr ArcsMat<M,N,T> operator/(const R& right) const{
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			ArcsMat<M,N,T> ret;
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) ret(j,i) = (*this)(j,i)/static_cast<T>(right);
				// 上記は次とほぼ同等→ ret.Data[i][j] = Data[i][j]/right;
			}
			return ret;
		}

		//! @brief 行列加算代入演算子(行列＝行列＋行列、行列＝行列＋スカラーの場合)
		//! @tparam	R	演算子右側の要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<typename R>
		constexpr ArcsMat<M,N,T>& operator+=(const R& right){
			(*this) = (*this) + right;	// 既に定義済みの加算演算子を利用
			return (*this);
		}
		
		//! @brief 行列減算代入演算子(行列＝行列－行列、行列＝行列－スカラーの場合)
		//! @tparam	R	演算子右側の要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<typename R>
		constexpr ArcsMat<M,N,T>& operator-=(const R& right){
			(*this) = (*this) - right;	// 既に定義済みの減算演算子を利用
			return (*this);
		}
		
		//! @brief 行列乗算代入演算子(行列＝行列＊行列、行列＝行列＊スカラーの場合)
		//! @tparam	R	演算子右側の要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<typename R>
		constexpr ArcsMat<M,N,T>& operator*=(const R& right){
			(*this) = (*this)*right;	// 既に定義済みの乗算演算子を利用
			return (*this);
		}
		
		//! @brief 行列除算代入演算子(行列＝行列／スカラーの場合)
		//! @tparam	R	演算子右側の要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<typename R>
		constexpr ArcsMat<M,N,T>& operator/=(const R& right){
			(*this) = (*this)/right;	// 既に定義済みの除算演算子を利用
			return (*this);
		}
		
		//! @brief 行列べき乗演算子(正方行列のべき乗)
		//! @param[in] right 演算子の右側
		//! @return 結果
		constexpr ArcsMat<M,N,T> operator^(const int& right) const{
			static_assert(M == N, "ArcsMat: Size Error");		// 正方行列チェック
			ArcsMat<M,N,T> ret = ArcsMat<M,N,T>::eye();
			if(0 <= right){
				// 非負のときはべき乗を計算
				for(size_t k = 1; k <= static_cast<size_t>(right); ++k) ret *= (*this);
			}else if(right == -1){
				// -1乗のときは逆行列を返す
				ret = ArcsMat<M,N,T>::inv((*this));
			}else{
				// -1より下の負数は未対応
				arcs_assert(false);
			}
			return ret;
		}
		
		//! @brief 行列アダマール積演算子(行列の要素ごとの乗算)
		//! @tparam	P, Q, R	演算子右側の行列の高さ, 幅, 要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<size_t P, size_t Q, typename R = double>
		constexpr ArcsMat<M,N,T> operator&(const ArcsMat<P,Q,R>& right) const{
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			ArcsMat<M,N,T> ret;
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) ret(j,i) = (*this)(j,i)*static_cast<T>( right(j,i) );
				// 上記は次とほぼ同等→ ret.Data[i][j] = Data[i][j]*right.Data[i][j];
			}
			return ret;
		}
		
		//! @brief 行列アダマール除算演算子 (行列の要素ごとの除算)
		//! @tparam	P, Q, R	演算子右側の行列の高さ, 幅, 要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<size_t P, size_t Q, typename R = double>
		constexpr ArcsMat<M,N,T> operator%(const ArcsMat<P,Q,R>& right) const{
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			ArcsMat<M,N,T> ret;
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) ret(j,i) = (*this)(j,i)/static_cast<T>( right(j,i) );
				// 上記は次とほぼ同等→ ret.Data[i][j] = Data[i][j]/right.Data[i][j];
			}
			return ret;
		}
		
		//! @brief 行列加算演算子 (スカラー＋行列の場合)
		//! @param[in] left		左側のスカラー値
		//! @param[in] right	右側の行列
		constexpr friend ArcsMat<M,N,T> operator+(const T& left, const ArcsMat<M,N,T>& right){
			return right + left;	// 既に定義済みの加算演算子を利用
		}
		
		//! @brief 行列減算演算子 (スカラー－行列の場合)
		//! @param[in] left		左側のスカラー値
		//! @param[in] right	右側の行列
		constexpr friend ArcsMat<M,N,T> operator-(const T& left, const ArcsMat<M,N,T>& right){
			return -right + left;	// 既に定義済みの単項マイナスと加算演算子を利用
		}
		
		//! @brief 行列乗算演算子 (スカラー＊行列の場合)
		//! @param[in] left		左側のスカラー値
		//! @param[in] right	右側の行列
		constexpr friend ArcsMat<M,N,T> operator*(const T& left, const ArcsMat<M,N,T>& right){
			return right*left;		// 既に定義済みの乗算演算子を利用
		}
		
		//! @brief 行列除算演算子 (スカラー／行列の場合)
		//! @param[in] left		左側のスカラー値
		//! @param[in] right	右側の行列
		constexpr friend void operator/(const T& left, const ArcsMat<M,N,T>& right){
			arcs_assert(false);		// この演算子は使用禁止
		}
		
		//! @brief 行列要素の各メモリアドレスを表示する関数
		constexpr void DispAddress(void) const{
			if(__builtin_constant_p(Data) == true) return;	// コンパイル時には処理を行わない

			printf("Mem addr of matrix:\n");
			for(size_t j = 0; j < M; ++j){
				printf("[ ");
				for(size_t i = 0; i < N; ++i){
					printf("%p ", &(Data[i][j]));
				}
				printf("]\n");
			}
			printf("\n");
		}
		
		//! @brief 行列の要素を表示
		//! @tparam	M, N, T	行列の高さ, 幅, 要素の型
		//! @param[in]	format	表示形式 (デフォルト値 "% g", "% 6.3e" とか "% 6.3f" などprintfの引数と同一)
		//! @param[in]	name	変数名 (デフォルト値 "")
		constexpr void Disp(const std::string& format = "% g", const std::string& name = "") const{
			if(__builtin_constant_p(Data) == true) return;	// コンパイル時には処理を行わない

			// 変数名の表示
			if(name != ""){
				printf("%s = \n", name.c_str());
			}

			// 各要素の表示
			for(size_t j = 0; j < M; ++j){
				printf("[ ");
				for(size_t i = 0; i < N; ++i){
					// データ型によって表示方法を変える
					if constexpr(ArcsMatrix::IsComplex<T>){
						// 複素数型の場合
						// 実数部の表示
						printf(format.c_str(), Data[i][j].real());
						// 虚数部の表示
						if(0.0 <= Data[i][j].imag()){
							printf(" +");
						}else{
							printf(" -");
						}
						printf( format.c_str(), std::abs(Data[i][j].imag()) );
						printf("%c", CMPLX_UNIT);
					}else{
						// それ以外の場合
						printf(format.c_str(), Data[i][j]);
					}
					printf(" ");
				}
				printf("]\n");
			}
			printf("\n");
		}

		//! @brief 行列のサイズを表示
		constexpr void DispSize(void) const{
			if(__builtin_constant_p(Data) == true) return;	// コンパイル時には処理を行わない

			printf("[ Height: %zu  x  Width: %zu ]\n\n", M, N);
		}
		
		//! @brief 行列の高さ(行数)を返す関数
		//! @return 行列の幅
		constexpr size_t GetHeight(void) const{
			return M;
		}
		
		//! @brief 行列の幅(列数)を返す関数
		//! @return 行列の幅
		constexpr size_t GetWidth(void) const{
			return N;
		}
		
		//! @brief 非ゼロ要素の数を返す関数
		//! @param[in]	eps	許容誤差 (デフォルト値 = EPSILON)
		//! @return 結果
		constexpr size_t GetNumOfNonZero(const T eps = ArcsMatrix::EPSILON) const{
			size_t ret = 0;
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j){
					if( std::abs(eps) < std::abs( (*this)(j,i) )) ++ret;
				}
			}
			return ret;
		}
		
		//! @brief 非ゼロの虚数部の要素の数を返す関数
		//! @param[in]	eps	許容誤差 (デフォルト値 = EPSILON)
		//! @return 結果
		constexpr size_t GetNumOfNonZeroImag(const T eps = ArcsMatrix::EPSILON) const{
			static_assert(ArcsMatrix::IsComplex<T>, "ArcsMat: Type Error (Need Complex)");
			size_t ret = 0;
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j){
					if( std::abs(eps) < std::abs( ((*this)(j,i)).imag() )) ++ret;
				}
			}
			return ret;
		}
		
		//! @brief 行列要素に値を設定する関数
		//! @tparam	T1, T2	要素の型
		//! @param[in]	u1	要素1の値
		//! @param[in]	u2	要素2以降の値
		template<typename T1, typename... T2>				// 可変長引数テンプレート
		constexpr void Set(const T1& u1, const T2&... u2){	// 再帰で順番に可変長引数を読み込んでいく
			static_assert(ArcsMatrix::IsApplicable<T1>, "ArcsMat: Type Error");	// 対応可能型チェック
			arcs_assert(Mindex < M);		// 縦方向カウンタが高さ以内かチェック
			arcs_assert(Nindex < N);		// 横方向カウンタが幅以内かチェック
			Data[Nindex][Mindex] = static_cast<T>(u1);	// キャストしてから行列の要素を埋める
			Nindex++;			// 横方向カウンタをインクリメント
			if(Nindex == N){	// 横方向カウンタが最後まで行き着いたら，
				Nindex = 0;		// 横方向カウンタを零に戻して，
				Mindex++;		// その代わりに，縦方向カウンタをインクリメント
			}
			Set(u2...);			// 自分自身を呼び出す(再帰)
		}
		constexpr void Set(){
			// 再帰の最後に呼ばれる関数
			Nindex = 0;	// すべての作業が終わったので，
			Mindex = 0;	// 横方向と縦方向カウンタを零に戻しておく
		}
		
		//! @brief 行列要素から値を読み込む関数
		//! @tparam	T1, T2	要素の型
		//! @param[in]	u1	要素1の値
		//! @param[in]	u2	要素2以降の値
		template<typename T1, typename... T2>	// 可変長引数テンプレート
		constexpr void Get(T1& u1, T2&... u2){	// 再帰で順番に可変長引数を読み込んでいく
			static_assert(ArcsMatrix::IsApplicable<T1>, "ArcsMat: Type Error");	// 対応可能型チェック
			arcs_assert(Mindex < M);		// 縦方向カウンタが高さ以内かチェック
			arcs_assert(Nindex < N);		// 横方向カウンタが幅以内かチェック
			u1 = static_cast<T1>(Data[Nindex][Mindex]);	// 行列の要素からキャストして読み込み
			Nindex++;			// 横方向カウンタをインクリメント
			if(Nindex == N){	// 横方向カウンタが最後まで行き着いたら，
				Nindex = 0;		// 横方向カウンタを零に戻して，
				Mindex++;		// その代わりに，縦方向カウンタをインクリメント
			}
			Get(u2...);			// 自分自身を呼び出す(再帰)
		}
		constexpr void Get(){
			// 再帰の最後に呼ばれる関数
			Nindex = 0;	// すべての作業が終わったので，
			Mindex = 0;	// 横方向と縦方向カウンタを零に戻しておく
		}
		
		//! @brief すべての要素を指定した値で埋める関数
		//! @tparam	R	要素の型
		//! @param[in] u 埋める値
		template<typename R>
		constexpr void FillAll(const R& u){
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j){
					// 型の種類によって値埋め処理を変更
					if constexpr(ArcsMatrix::IsIntFloat<R> && ArcsMatrix::IsIntFloat<T>){
						// 「整数・浮動小数点型 → 整数・浮動小数点型」の値埋め処理
						Data[i][j] = static_cast<T>(u);
					}else if constexpr(ArcsMatrix::IsIntFloat<R> && ArcsMatrix::IsComplex<T>){
						// 「整数・浮動小数点型 → 複素数型」の値埋め処理
						Data[i][j].real(u);
						Data[i][j].imag(0);
					}else if constexpr(ArcsMatrix::IsComplex<R> && ArcsMatrix::IsComplex<T>){
						// 「複素数型 → 複素数型」の値埋め処理
						Data[i][j] = static_cast<T>(u);
					}else{
						// 値埋め不能の場合の処理
						arcs_assert(false);	// 変換不能、もしここに来たら問答無用でAssertion Failed
					}
				}
			}
		}
		
		//! @brief すべての要素を指定したゼロで埋める関数
		constexpr void FillAllZero(void){
			FillAll(0);
		}
		
		//! @brief 1次元std::array配列を縦ベクトルとして読み込む関数
		//! @tparam	P, R	配列の長さ, 要素の型
		//! @param[in]	Array	std::array配列(縦MM×横1)
		template<size_t P, typename R = double>
		constexpr void LoadArray(const std::array<R, P>& Array){
			static_assert(N == 1, "ArcsMat: Vector Error");	// 縦ベクトルチェック
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			for(size_t j = 0; j < M; ++j) Data[0][j] = Array[j];
		}
		
		//! @brief 縦ベクトルを1次元std::array配列に書き込む関数
		//! @tparam	P, R	配列の長さ, 要素の型
		//! @param[out]	Array	std::array配列(縦MM×横1)
		template<size_t P, typename R = double>
		constexpr void StoreArray(std::array<R, P>& Array) const{
			static_assert(N == 1, "ArcsMat: Vector Error");	// 縦ベクトルチェック
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			for(size_t j = 0; j < M; ++j) Array[j] = Data[0][j];
		}
		
		//! @brief 2次元std::array配列を行列として読み込む関数
		//! @tparam	P, Q, R	配列の高さ, 幅, 要素の型
		//! @param[in]	Array	std::array配列(縦MM×横NN)
		template<size_t P, size_t Q, typename R = double>
		constexpr void LoadArray(const std::array<std::array<R, P>, Q>& Array){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			Data = Array;
		}
		
		//! @brief 行列を2次元std::array配列に書き込む関数
		//! @tparam	P, Q, R	配列の高さ, 幅, 要素の型
		//! @param[in]	Array	std::array配列(縦MM×横NN)
		template<size_t P, size_t Q, typename R = double>
		constexpr void StoreArray(std::array<std::array<R, P>, Q>& Array) const{
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			Array = Data;
		}
		
		//! @brief 行列の状態をそのまま返す関数
		//! @return 行列の状態
		constexpr ArcsMatrix::MatStatus GetStatus(void) const{
			return Status;
		}

		//! @brief std:arrayの2次元配列データをそのまま返す関数
		//! @return	2次元配列データ
		constexpr std::array<std::array<T, M>, N> GetData(void) const{
			return Data;
		}

		//! @brief std:arrayの2次元配列データの読み込み専用の参照を返す関数
		//! @return	2次元配列データの参照
		constexpr const std::array<std::array<T, M>, N>& ReadOnlyRef(void) const{
			return Data;
		}
		
		//! @brief 指定した先頭位置から縦ベクトルを抜き出して返す関数 (引数渡し版)
		//! @tparam	P, Q, R	縦ベクトルの高さ, 幅, 要素の型
		//! @param[out]	v	縦ベクトル
		//! @param[in]	m	先頭位置 m行目
		//! @param[in]	n	先頭位置 n列目
		template<size_t P, size_t Q, typename R = double>
		constexpr void GetVerticalVec(ArcsMat<P,Q,R>& v, const size_t m, const size_t n) const{
			static_assert(Q == 1, "ArcsMat: Vector Error");			// 縦ベクトルチェック
			static_assert(P <= M, "ArcsMat: Vector Size Error");	// ベクトルサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			arcs_assert(0 < m && P + m - 1 <= M);	// はみ出しチェック
			arcs_assert(0 < n && n <= N);			// サイズチェック
			for(size_t j = 1; j <= P; ++j) v(j,1) = (*this)(m + j - 1, n);
		}
		
		//! @brief 指定した先頭位置から縦ベクトルを抜き出して返す関数 (戻り値渡し版)
		//! @tparam	P	縦ベクトルの高さ
		//! @param[in]	m	先頭位置 m行目
		//! @param[in]	n	先頭位置 n列目
		//! @return	縦ベクトル
		template<size_t P>
		constexpr ArcsMat<P,1,T> GetVerticalVec(const size_t m, const size_t n) const{
			ArcsMat<P,1,T> ret;
			GetVerticalVec(ret, m, n);
			return ret;
		}
		
		//! @brief 指定した先頭位置から横ベクトルを抜き出して返す関数 (引数渡し版)
		//! @tparam	P, Q, R	横ベクトルの高さ, 幅, 要素の型
		//! @param[out]	w	横ベクトル
		//! @param[in]	m	先頭位置 m行目
		//! @param[in]	n	先頭位置 n列目
		template<size_t P, size_t Q, typename R = double>
		constexpr void GetHorizontalVec(ArcsMat<P,Q,R>& w, const size_t m, const size_t n) const{
			static_assert(P == 1, "ArcsMat: Vector Error");			// 横ベクトルチェック
			static_assert(Q <= N, "ArcsMat: Vector Size Error");	// ベクトルサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			arcs_assert(0 < m && m <= M);			// サイズチェック
			arcs_assert(0 < n && Q + n - 1 <= N);	// はみ出しチェック
			for(size_t i = 1; i <= Q; ++i) w(1,i) = (*this)(m, n + i - 1);
		}
		
		//! @brief 指定した先頭位置から横ベクトルを抜き出して返す関数 (戻り値渡し版)
		//! @tparam	Q	横ベクトルの幅
		//! @param[in]	m	先頭位置 m行目
		//! @param[in]	n	先頭位置 n列目
		//! @return	横ベクトル
		template<size_t Q>
		constexpr ArcsMat<1,Q,T> GetHorizontalVec(const size_t m, const size_t n) const{
			ArcsMat<1,Q,T> ret;
			GetHorizontalVec(ret, m, n);
			return ret;
		}
		
		//! @brief 指定した先頭位置に縦ベクトルを埋め込む関数
		//! @tparam	P, Q, R	縦ベクトルの高さ, 幅, 要素の型
		//! @param[in]	v	縦ベクトル
		//! @param[in]	m	先頭位置 m行目
		//! @param[in]	n	先頭位置 n列目
		template<size_t P, size_t Q, typename R = double>
		constexpr void SetVerticalVec(const ArcsMat<P,Q,R>& v, const size_t m, const size_t n){
			static_assert(Q == 1, "ArcsMat: Vector Error");			// 縦ベクトルチェック
			static_assert(P <= M, "ArcsMat: Vector Size Error");	// ベクトルサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			arcs_assert(0 < m && P + m - 1 <= M);	// はみ出しチェック
			arcs_assert(0 < n && n <= N);			// サイズチェック
			for(size_t j = 1; j <= P; ++j) (*this)(m + j - 1, n) = v(j,1);
		}
		
		//! @brief 指定した先頭位置に横ベクトルを埋め込む関数
		//! @tparam	P, Q, R	縦ベクトルの高さ, 幅, 要素の型
		//! @param[in]	w	横ベクトル
		//! @param[in]	m	先頭位置 m行目
		//! @param[in]	n	先頭位置 n列目
		template<size_t P, size_t Q, typename R = double>
		constexpr void SetHorizontalVec(const ArcsMat<P,Q,R>& w, size_t m, size_t n){
			static_assert(P == 1, "ArcsMat: Vector Error");		// 横ベクトルチェック
			static_assert(Q <= N, "ArcsMat: Vector Size Error");	// ベクトルサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			arcs_assert(0 < m && m <= M);			// サイズチェック
			arcs_assert(0 < n && Q + n - 1 <= N);	// はみ出しチェック
			for(size_t i = 1; i <= Q; ++i) (*this)(m, n + i - 1) = w(1,i);
		}

		//! @brief 縦方向の要素を逆順に入れ替える関数
		constexpr void FlipVertical(void){
			ArcsMat<M,N,T> buff;
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j){
					buff(j, i) = (*this)(M - j + 1, i);
				}
			}
			(*this) = buff;
		}

		//! @brief 横方向の要素を逆順に入れ替える関数
		constexpr void FlipHorizontal(void){
			ArcsMat<M,N,T> buff;
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j){
					buff(j, i) = (*this)(j, N - i + 1);
				}
			}
			(*this) = buff;
		}
		
		//! @brief ゼロに近い要素を完全にゼロにする関数
		//! @param[in]	eps	許容誤差 (デフォルト値 = EPSILON)
		constexpr void Zeroing(const T eps = ArcsMatrix::EPSILON){
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j){
					if(std::abs( (*this)(j,i) ) < std::real(eps)) (*this)(j,i) = 0;
				}
			}
		}

		//! @brief ゼロに近い虚数部の要素を完全にゼロにする関数
		//! @param[in]	eps	許容誤差 (デフォルト値 = EPSILON)
		constexpr void ZeroingImag(const T eps = ArcsMatrix::EPSILON){
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j){
					if(std::abs( ((*this)(j,i)).imag() ) < std::real(eps)){
						((*this)(j,i)).imag(0);
					}
				}
			}
		}
		
		//! @brief 下三角(主対角除く)に限定して、ゼロに近い要素を完全にゼロにする関数
		//! @param[in]	eps	許容誤差 (デフォルト値 = EPSILON)
		constexpr void ZeroingTriLo(const T eps = ArcsMatrix::EPSILON){
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = i + 1; j <= M; ++j){
					if constexpr(ArcsMatrix::IsComplex<T>){
						//if(std::abs( std::real((*this)(j,i)) ) < std::real(eps)) ((*this)(j,i)).real(0);
						//if(std::abs( std::imag((*this)(j,i)) ) < std::real(eps)) ((*this)(j,i)).imag(0);
						if(std::abs( (*this)(j,i) ) < std::real(eps)) (*this)(j,i) = 0;
					}else{
						if(std::abs( (*this)(j,i) ) < eps) (*this)(j,i) = 0;
					}
				}
			}
		}

	public:
		// ここから下は静的メンバ関数

		//! @brief m行n列の零行列を返す関数
		//! @return 零行列
		static constexpr ArcsMat<M,N,T> zeros(void){
			ArcsMat<M,N,T> ret;
			return ret;
		}

		//! @brief m行n列の要素がすべて1の行列を返す関数
		//! @return 1行列
		static constexpr ArcsMat<M,N,T> ones(void){
			ArcsMat<M,N,T> ret;
			ret.FillAll(1);
			return ret;
		}
		
		//! @brief n行n列の単位行列を返す関数
		//! @return 単位行列
		static constexpr ArcsMat<M,N,T> eye(void){
			static_assert(M == N, "ArcsMat: Size Error");	// 正方行列チェック
			ArcsMat<M,N,T> ret;
			if constexpr(ArcsMatrix::IsComplex<T>){
				// 複素数型の場合
				for(size_t i = 1; i <= N; ++i) ret(i,i) = std::complex(1.0, 0.0);	// 対角成分を 1 + j0 で埋める
			}else{
				// それ以外の場合
				for(size_t i = 1; i <= N; ++i) ret(i,i) = static_cast<T>(1);		// 対角成分を 1 で埋める
			}
			return ret;
		}

		//! @brief 単調増加の縦ベクトルを返す関数
		//! @return 1～MMまでの単調増加の縦ベクトル
		static constexpr ArcsMat<M,N,T> ramp(void){
			static_assert(N == 1, "ArcsMat: Vector Error");// 行列のサイズチェック
			ArcsMat<M,N,T> ret;
			for(size_t j = 1; j <= M; ++j) ret[j] = j;		// 単調増加を書き込む
			return ret;
		}
		
		//! @brief 指定した列から縦ベクトルとして抽出する関数 (引数渡し版)
		//! @tparam	P, Q, R	出力ベクトルの高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	y	出力ベクトル
		//! @param[in]	n	抽出したい列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void getcolumn(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& y, const size_t n){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			y.SetVerticalVec(U.GetVerticalVec<M>(1,n) ,1, 1);
		}
		
		//! @brief 指定した列から縦ベクトルとして抽出する関数 (戻り値渡し版)
		//! @param[in]	U	入力行列
		//! @param[in]	n	抽出したい列
		//! @return	出力ベクトル
		static constexpr ArcsMat<M,1,T> getcolumn(const ArcsMat<M,N,T>& U, const size_t n){
			return U.GetVerticalVec<M>(1, n);
		}
		
		//! @brief 指定した列を縦ベクトルで上書きする関数 (引数渡し版)
		//! @tparam	P, Q, R	入力ベクトルの高さ, 幅, 要素の型
		//! @param[in,out]	UY	入出力行列
		//! @param[in]	u	入力ベクトル
		//! @param[in]	n	上書きしたい列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void setcolumn(ArcsMat<M,N,T>& UY, const ArcsMat<P,Q,R>& u, const size_t n){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			UY.SetVerticalVec(u, 1, n);
		}
		
		//! @brief 指定した列を縦ベクトルで上書きする関数 (戻り値渡し版)
		//! @tparam	P, Q, R	入力ベクトルの高さ, 幅, 要素の型
		//! @param[in]	u	入力ベクトル
		//! @param[in]	n	上書きしたい列
		//! @param[in]	U	入力行列
		//! @return	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr ArcsMat<M,N,T> setcolumn(const ArcsMat<P,Q,R>& u, const size_t n, const ArcsMat<M,N,T>& U){
			ArcsMat<M,N,T> Y = U;
			setcolumn(Y, u, n);
			return Y;
		}
		
		//! @brief 指定した列と列を入れ替える関数 (引数渡し版)
		//! @param[in,out]	UY	入出力行列
		//! @param[in]	n1	指定列1
		//! @param[in]	n2	指定列2
		static constexpr void swapcolumn(ArcsMat<M,N,T>& UY, const size_t n1, const size_t n2){
			ArcsMat<M,1,T> p, q;		// バッファ用ベクトル
			p = getcolumn(UY, n1);		// n1列目を抽出
			q = getcolumn(UY, n2);		// n2列目を抽出
			setcolumn(UY, p, n2);		// m2列目に書き込み
			setcolumn(UY, q, n1);		// m1列目に書き込み
		}
		
		//! @brief 指定した列と列を入れ替える関数 (戻り値渡し版)
		//! @param[in]	n1	指定列1
		//! @param[in]	n2	指定列2
		//! @param[in]	U	入力行列
		//! @return	Y	出力行列
		static constexpr ArcsMat<M,N,T> swapcolumn(const size_t n1, const size_t n2, const ArcsMat<M,N,T>& U){
			ArcsMat<M,N,T> Y = U;
			swapcolumn(Y, n1, n2);
			return Y;
		}
		
		//! @brief n列目のm1行目からm2行目までを数値aで埋める関数 (m1 <= m2 であること) (引数渡し版)
		//! @tparam	R	埋める値の型
		//! @param[in,out]	UY	入出力行列
		//! @param[in]	a	埋める値
		//! @param[in]	n	指定列
		//! @param[in]	m1	開始行
		//! @param[in]	m2	終了行
		template<typename R = double>
		static constexpr void fillcolumn(ArcsMat<M,N,T>& UY, const R a, const size_t n, const size_t m1, const size_t m2){
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			arcs_assert(0 < n && n <= N);	// 列が幅以内かチェック
			arcs_assert(0 < m1 && m1 <= M);	// 行が高さ以内かチェック
			arcs_assert(0 < m2 && m2 <= M);	// 行が高さ以内かチェック
			arcs_assert(m1 <= m2);			// 開始行と終了行が入れ替わらないかチェック
			for(size_t j = m1; j <= m2; ++j) UY(j, n) = (T)a;
		}
		
		//! @brief n列目のm1行目からm2行目までを数値aで埋める関数 (n1 <= n2 であること) (戻り値渡し版)
		//! @tparam	R	埋める値の型
		//! @param[in]	a	埋める値
		//! @param[in]	n	指定列
		//! @param[in]	m1	開始行
		//! @param[in]	m2	終了行
		//! @param[in]	U	入力行列
		//! @return	出力行列
		template<typename R = double>
		static constexpr ArcsMat<M,N,T> fillcolumn(const R a, const size_t n, const size_t m1, const size_t m2, const ArcsMat<M,N,T>& U){
			ArcsMat<M,N,T> Y = U;
			fillcolumn(Y, a, n, m1, m2);
			return Y;
		}
		
		//! @brief 並び替え指定横ベクトルuが昇順になるように，行列Uの列を並び替える関数 (戻り値渡し版のみ)
		//! @tparam	P, Q, R	並び替え指定横ベクトルの高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[in]	u	並び替え指定横ベクトル
		//! @return	出力行列
		template<size_t P, size_t Q, typename R = size_t>
		static constexpr ArcsMat<M,N,T> ordercolumn(const ArcsMat<M,N,T>& U, const ArcsMat<P,Q,R>& u){
			static_assert(P == 1, "ArcsMat: Vector Error");	// 横ベクトルチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_integral_v<R>, "ArcsMat: Input u should be integer type.");	// 整数型チェック
			ArcsMat<M,N,T> Y;
			for(size_t i = 1; i <= N; ++i){
				setcolumn(Y,  getcolumn(U, static_cast<size_t>( u(1,i) )), i);
			}
			return Y;
		}
		
		//! @brief 並び替え指定横ベクトルuが昇順になるように，行列Uの列と指定横ベクトルの両方を並び替える関数 (引数渡し版のみ)
		//! @tparam	P, Q, R	並び替え指定横ベクトルの高さ, 幅, 要素の型
		//! @param[in,out]	UY	入出力行列
		//! @param[in,out]	uy	並び替え指定横ベクトル
		template<size_t P, size_t Q, typename R = size_t>
		static constexpr void ordercolumn_and_vec(ArcsMat<M,N,T>& UY, ArcsMat<P,Q,R>& uy){
			static_assert(P == 1, "ArcsMat: Vector Error");	// 横ベクトルチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_integral_v<R>, "ArcsMat: Input u should be integer type.");	// 整数型チェック
			for(size_t i = 1; i <= N; ++i){
				swapcolumn(UY, i, uy(1,i));
				ArcsMat<P,Q,R>::swapcolumn(uy, i, static_cast<size_t>( uy(1,i) ));
			}
		}
		
		//! @brief 各列の総和を計算して横ベクトルを出力する関数 (引数渡し版)
		//! @tparam	P, Q, R	出力ベクトルの高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	y	出力ベクトル
		template<size_t P, size_t Q, typename R = double>
		static constexpr void sumcolumn(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& y){
			static_assert(P == 1, "ArcsMat: Vector Error");	// 横ベクトルチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			y.FillAllZero();
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) y(1,i) += static_cast<R>( U(j,i) );
			}
		}
		
		//! @brief 各列の総和を計算して横ベクトルを出力する関数 (戻り値渡し版)
		//! @param[in]	U	入力行列
		//! @return	出力ベクトル
		static constexpr ArcsMat<1,N,T> sumcolumn(const ArcsMat<M,N,T>& U){
			ArcsMat<1,N,T> y;
			sumcolumn(U, y);
			return y;
		}
		
		//! @brief 指定した行から横ベクトルとして抽出する関数 (引数渡し版)
		//! @tparam	P, Q, R	出力ベクトルの高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	y	出力ベクトル
		//! @param[in]	m	抽出したい行
		template<size_t P, size_t Q, typename R = double>
		static constexpr void getrow(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& y, const size_t m){
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			y.SetHorizontalVec(U.GetHorizontalVec<N>(m,1) ,1, 1);
		}
		
		//! @brief 指定した行から横ベクトルとして抽出する関数 (戻り値渡し版)
		//! @param[in]	U	入力行列
		//! @param[in]	m	抽出したい行
		//! @return	出力ベクトル
		static constexpr ArcsMat<1,N,T> getrow(const ArcsMat<M,N,T>& U, const size_t m){
			return U.GetHorizontalVec<N>(m, 1);
		}
		
		//! @brief 指定した行を横ベクトルで上書きする関数 (引数渡し版)
		//! @tparam	P, Q, R	入力ベクトルの高さ, 幅, 要素の型
		//! @param[in,out]	UY	入出力行列
		//! @param[in]	u	入力ベクトル
		//! @param[in]	m	上書きしたい行
		template<size_t P, size_t Q, typename R = double>
		static constexpr void setrow(ArcsMat<M,N,T>& UY, const ArcsMat<P,Q,R>& u, const size_t m){
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			UY.SetHorizontalVec(u, m, 1);
		}
		
		//! @brief 指定した行を横ベクトルで上書きする関数 (戻り値渡し版)
		//! @tparam	P, Q, R	入力ベクトルの高さ, 幅, 要素の型
		//! @param[in]	u	入力ベクトル
		//! @param[in]	m	上書きしたい行
		//! @param[in]	U	入力行列
		//! @return	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr ArcsMat<M,N,T> setrow(const ArcsMat<P,Q,R>& u, const size_t m, const ArcsMat<M,N,T>& U){
			ArcsMat<M,N,T> Y = U;
			setrow(Y, u, m);
			return Y;
		}
		
		//! @brief 指定した行と行を入れ替える関数 (引数渡し版)
		//! @param[in,out]	UY	入出力行列
		//! @param[in]	m1	指定行1
		//! @param[in]	m2	指定行2
		static constexpr void swaprow(ArcsMat<M,N,T>& UY, const size_t m1, const size_t m2){
			ArcsMat<1,N,T> p, q;	// バッファ用ベクトル
			p = getrow(UY, m1);		// m1行目を抽出
			q = getrow(UY, m2);		// m2行目を抽出
			setrow(UY, p, m2);		// m2行目に書き込み
			setrow(UY, q, m1);		// m1行目に書き込み
		}
		
		//! @brief 指定した行と行を入れ替える関数 (戻り値渡し版)
		//! @param[in]	m1	指定行1
		//! @param[in]	m2	指定行2
		//! @param[in]	U	入力行列
		//! @return	Y	出力行列
		static constexpr ArcsMat<M,N,T> swaprow(const size_t m1, const size_t m2, const ArcsMat<M,N,T>& U){
			ArcsMat<M,N,T> Y = U;
			swaprow(Y, m1, m2);
			return Y;
		}
		
		//! @brief m行目のn1列目からn2列目までを数値aで埋める関数 (n1 <= n2 であること) (引数渡し版)
		//! @tparam	R	埋める値の型
		//! @param[in,out]	UY	入出力行列
		//! @param[in]	a	埋める値
		//! @param[in]	m	指定行
		//! @param[in]	n1	開始列
		//! @param[in]	n2	終了列
		template<typename R = double>
		static constexpr void fillrow(ArcsMat<M,N,T>& UY, const R a, const size_t m, const size_t n1, const size_t n2){
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			arcs_assert(0 < m && m <= M);		// 行が高さ以内かチェック
			arcs_assert(0 < n1 && n1 <= N);	// 列が幅以内かチェック
			arcs_assert(0 < n2 && n2 <= N);	// 列が幅以内かチェック
			arcs_assert(n1 <= n2);				// 開始列と終了列が入れ替わらないかチェック
			for(size_t i = n1; i <= n2; ++i) UY(m, i) = static_cast<T>(a);
		}
		
		//! @brief m行目のn1列目からn2列目までを数値aで埋める関数 (n1 <= n2 であること) (戻り値渡し版)
		//! @tparam	R	埋める値の型
		//! @param[in]	a	埋める値
		//! @param[in]	m	指定行
		//! @param[in]	n1	開始列
		//! @param[in]	n2	終了列
		//! @param[in]	U	入力行列
		//! @return	出力行列
		template<typename R = double>
		static constexpr ArcsMat<M,N,T> fillrow(const R a, const size_t m, const size_t n1, const size_t n2, const ArcsMat<M,N,T>& U){
			ArcsMat<M,N,T> Y = U;
			fillrow(Y, a, m, n1, n2);
			return Y;
		}
		
		//! @brief 並び替え指定縦ベクトルuが昇順になるように，行列Uの行を並び替える関数 (戻り値渡し版のみ)
		//! @tparam	P, Q, R	並び替え指定縦ベクトルの高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[in]	u	並び替え指定縦ベクトル
		//! @return	出力行列
		template<size_t P, size_t Q, typename R = size_t>
		static constexpr ArcsMat<M,N,T> orderrow(const ArcsMat<M,N,T>& U, const ArcsMat<P,Q,R>& u){
			static_assert(Q == 1, "ArcsMat: Vector Error");	// 縦ベクトルチェック
			static_assert(M == P, "ArcsMat: Size Error");		// 行列のサイズチェック
			static_assert(std::is_integral_v<R>, "ArcsMat: Input u should be integer type.");	// 整数型チェック
			ArcsMat<M,N,T> Y;
			for(size_t j = 1; j <= M; ++j){
				setrow(Y,  getrow(U, static_cast<size_t>( u(j,1) )), j);
			}
			return Y;
		}
		
		//! @brief 並び替え指定縦ベクトルuが昇順になるように，行列Uの行と指定縦ベクトルの両方を並び替える関数 (引数渡し版のみ)
		//! @tparam	P, Q, R	並び替え指定縦ベクトルの高さ, 幅, 要素の型
		//! @param[in,out]	UY	入出力行列
		//! @param[in,out]	uy	並び替え指定縦ベクトル
		template<size_t P, size_t Q, typename R = size_t>
		static constexpr void orderrow_and_vec(ArcsMat<M,N,T>& UY, ArcsMat<P,Q,R>& uy){
			static_assert(Q == 1, "ArcsMat: Vector Error");	// 縦ベクトルチェック
			static_assert(M == P, "ArcsMat: Size Error");		// 行列のサイズチェック
			static_assert(std::is_integral_v<R>, "ArcsMat: Input u should be integer type.");	// 整数型チェック
			for(size_t j = 1; j <= M; ++j){
				swaprow(UY, j, uy(j,1));
				ArcsMat<P,Q,R>::swaprow(uy, j, static_cast<size_t>( uy(j,1) ));
			}
		}
		
		//! @brief 各行の総和を計算して縦ベクトルを出力する関数 (引数渡し版)
		//! @tparam	P, Q, R	出力ベクトルの高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	y	出力ベクトル
		template<size_t P, size_t Q, typename R = double>
		static constexpr void sumrow(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& y){
			static_assert(Q == 1, "ArcsMat: Vector Error");	// 縦ベクトルチェック
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			y.FillAllZero();
			for(size_t j = 1; j <= M; ++j){
				for(size_t i = 1; i <= N; ++i) y(j,1) += static_cast<R>( U(j,i) );
			}
		}
		
		//! @brief 各行の総和を計算して縦ベクトルを出力する関数 (戻り値渡し版)
		//! @param[in]	U	入力行列
		//! @return	出力ベクトル
		static constexpr ArcsMat<M,1,T> sumrow(const ArcsMat<M,N,T>& U){
			ArcsMat<M,1,T> y;
			sumrow(U, y);
			return y;
		}
		
		//! @brief 指定位置から縦ベクトルを抽出する関数 (引数渡し版)
		//! @tparam	P, Q, R	出力ベクトルの高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[in]	y	出力ベクトル
		//! @param[in]	m	先頭位置 m行目
		//! @param[in]	n	先頭位置 n列目
		template<size_t P, size_t Q, typename R = double>
		static constexpr void getvvector(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& y, const size_t m, const size_t n){
			U.GetVerticalVec(y, m, n);
		}
		
		//! @brief 指定位置から縦ベクトルを抽出する関数 (戻り値渡し版)
		//! @tparam	P	出力ベクトルの高さ
		//! @param[in]	U	入力行列
		//! @param[in]	m	先頭位置 m行目
		//! @param[in]	n	先頭位置 n列目
		//! @return	縦ベクトル
		template<size_t P = M>
		static constexpr ArcsMat<P,1,T> getvvector(const ArcsMat<M,N,T>& U, const size_t m, const size_t n){
			return U.GetVerticalVec<P>(m, n);
		}
		
		//! @brief 指定位置に縦ベクトルで上書きする関数 (引数渡し版)
		//! @tparam	P, Q, R	入力ベクトルの高さ, 幅, 要素の型
		//! @param[in,out]	UY	入出力行列
		//! @param[in]	u	入力ベクトル
		//! @param[in]	m	先頭位置 m行目
		//! @param[in]	n	先頭位置 n列目
		template<size_t P, size_t Q, typename R = double>
		static constexpr void setvvector(ArcsMat<M,N,T>& UY, const ArcsMat<P,Q,R>& u, const size_t m, const size_t n){
			UY.SetVerticalVec(u, m, n);
		}
		
		//! @brief 指定位置に縦ベクトルで上書きする関数 (戻り値渡し版)
		//! @tparam	P, Q, R	入力ベクトルの高さ, 幅, 要素の型
		//! @param[in]	u	入力ベクトル
		//! @param[in]	m	先頭位置 m行目
		//! @param[in]	n	先頭位置 n列目
		//! @param[in]	U	入力行列
		//! @return	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr ArcsMat<M,N,T> setvvector(const ArcsMat<P,Q,R>& u, const size_t m, const size_t n, const ArcsMat<M,N,T>& U){
			ArcsMat<M,N,T> Y = U;
			Y.SetVerticalVec(u, m, n);
			return Y;
		}
		
		//! @brief 指定位置から横ベクトルを抽出する関数 (引数渡し版)
		//! @tparam	P, Q, R	出力ベクトルの高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[in]	y	出力ベクトル
		//! @param[in]	m	先頭位置 m行目
		//! @param[in]	n	先頭位置 n列目
		template<size_t P, size_t Q, typename R = double>
		static constexpr void gethvector(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& y, const size_t m, const size_t n){
			U.GetHorizontalVec(y, m, n);
		}
		
		//! @brief 指定位置から横ベクトルを抽出する関数 (戻り値渡し版)
		//! @tparam	Q	出力ベクトルの幅
		//! @param[in]	U	入力行列
		//! @param[in]	m	先頭位置 m行目
		//! @param[in]	n	先頭位置 n列目
		//! @return	横ベクトル
		template<size_t Q>
		static constexpr ArcsMat<1,Q,T> gethvector(const ArcsMat<M,N,T>& U, const size_t m, const size_t n){
			return U.GetHorizontalVec<Q>(m, n);
		}
		
		//! @brief 指定位置に横ベクトルで上書きする関数 (引数渡し版)
		//! @tparam	P, Q, R	入力ベクトルの高さ, 幅, 要素の型
		//! @param[in,out]	UY	入出力行列
		//! @param[in]	u	入力ベクトル
		//! @param[in]	m	先頭位置 m行目
		//! @param[in]	n	先頭位置 n列目
		template<size_t P, size_t Q, typename R = double>
		static constexpr void sethvector(ArcsMat<M,N,T>& UY, const ArcsMat<P,Q,R>& u, const size_t m, const size_t n){
			UY.SetHorizontalVec(u, m, n);
		}
		
		//! @brief 指定位置に横ベクトルで上書きする関数 (戻り値渡し版)
		//! @tparam	P, Q, R	入力ベクトルの高さ, 幅, 要素の型
		//! @param[in]	u	入力ベクトル
		//! @param[in]	m	先頭位置 m行目
		//! @param[in]	n	先頭位置 n列目
		//! @param[in]	U	入力行列
		//! @return	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr ArcsMat<M,N,T> sethvector(const ArcsMat<P,Q,R>& u, const size_t m, const size_t n, const ArcsMat<M,N,T>& U){
			ArcsMat<M,N,T> Y = U;
			Y.SetHorizontalVec(u, m, n);
			return Y;
		}
		
		//! @brief 行列から指定位置の小行列を抽出する関数 (引数渡し版)
		//! @tparam	P, Q, R	小行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	抽出した小行列
		//! @param[in]	m	抽出する縦位置（小行列の上）
		//! @param[in]	n	抽出する横位置（小行列の左）
		template<size_t P, size_t Q, typename R = double>
		static constexpr void getsubmatrix(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y, const size_t m, const size_t n){
			static_assert(P <= M, "ArcsMat: Size Error");	// 小行列の方が高さが小さいかチェック
			static_assert(Q <= N, "ArcsMat: Size Error");	// 小行列の方が幅が小さいかチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			arcs_assert(P + m - 1 <= M);	// 下側がハミ出ないかチェック
			arcs_assert(Q + n - 1 <= N);	// 右側がハミ出ないかチェック
			for(size_t i = 1; i <= Q; ++i) ArcsMat<P,Q,R>::setcolumn(Y, getvvector<P>(U, m, n + i - 1), i);
		}
		
		//! @brief 行列から指定位置の小行列を抽出する関数 (戻り値渡し版)
		//! @tparam	P, Q, 出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[in]	m	抽出する縦位置（小行列の上）
		//! @param[in]	n	抽出する横位置（小行列の左）
		//! @return	抽出した小行列
		template<size_t P, size_t Q>
		static constexpr ArcsMat<P,Q,T> getsubmatrix(const ArcsMat<M,N,T>& U, const size_t m, const size_t n){
			ArcsMat<P,Q,T> Y;
			getsubmatrix(U, Y, m, n);
			return Y;
		}
		
		//! @brief 小行列を行列の指定位置に上書きする関数 (引数渡し版)
		//! @tparam	P, Q, R	小行列の高さ, 幅, 要素の型
		//! @param[in,out]	UY	入出力行列
		//! @param[in]	U	書き込む小行列
		//! @param[in]	m	上書きしたい縦位置（小行列の上）
		//! @param[in]	n	上書きしたい横位置（小行列の左）
		template<size_t P, size_t Q, typename R = double>
		static constexpr void setsubmatrix(ArcsMat<M,N,T>& UY, const ArcsMat<P,Q,R>& U, const size_t m, const size_t n){
			static_assert(P <= M, "ArcsMat: Size Error");	// 小行列の方が高さが小さいかチェック
			static_assert(Q <= N, "ArcsMat: Size Error");	// 小行列の方が幅が小さいかチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			arcs_assert(P + m - 1 <= M);	// 下側がハミ出ないかチェック
			arcs_assert(Q + n - 1 <= N);	// 右側がハミ出ないかチェック
			for(size_t i = 1; i <= Q; ++i) setvvector(UY, ArcsMat<P,Q,R>::getcolumn(U, i), m, i + n - 1);
		}
		
		//! @brief 小行列を行列の指定位置に上書きする関数 (戻り値渡し版)
		//! @tparam	P, Q, R	小行列の高さ, 幅, 要素の型
		//! @param[in]	Us	書き込む小行列
		//! @param[in]	m	上書きしたい縦位置（小行列の上）
		//! @param[in]	n	上書きしたい横位置（小行列の左）
		//! @param[in]	U	入力行列
		//! @return	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr ArcsMat<M,N,T> setsubmatrix(const ArcsMat<P,Q,R>& Us, const size_t m, const size_t n, const ArcsMat<M,N,T>& U){
			ArcsMat<M,N,T> Y = U;
			setsubmatrix(Y, Us, m, n);
			return Y;
		}
		
		//! @brief 行列Uから別の行列Yへ位置とサイズを指定してコピーする関数(引数渡し版のみ)
		//! 等価なMATLABコード: Y(my:my+(m2-m1), ny:ny+(n2-n1)) = U(m1:m2, n1:n2)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型, 入力行列の高さ, 幅, 要素の型
		//! @param[in]	U	コピー元の行列
		//! @param[in]	m1	コピー元の縦方向の抽出開始位置
		//! @param[in]	m2	コピー元の縦方向の抽出終了位置
		//! @param[in]	n1	コピー元の横方向の抽出開始位置
		//! @param[in]	n2	コピー元の横方向の抽出終了位置
		//! @param[in,out]	Y	コピー先の行列
		//! @param[in]	my	コピー先の縦方向の書き込み開始位置
		//! @param[in]	ny	コピー先の横方向の書き込み開始位置
		template<size_t P, size_t Q, typename R = double>
		static constexpr void copymatrix(
			const ArcsMat<M,N,T>& U, const size_t m1, const size_t m2, const size_t n1, const size_t n2,
			ArcsMat<P,Q,R>& Y, const size_t my, const size_t ny
		){
			static_assert(ArcsMatrix::IsApplicable<T>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			arcs_assert(m1 <= m2);	// サイズチェック
			arcs_assert(n1 <= n2);	// サイズチェック
			arcs_assert(m2 <= M);	// サイズチェック
			arcs_assert(n2 <= N);	// サイズチェック
			arcs_assert(my+(m2-m1) <= P);	// サイズチェック
			arcs_assert(ny+(n2-n1) <= Q);	// サイズチェック

			for(size_t i = 0; i <= (n2 - n1); ++i){
				for(size_t j = 0; j <= (m2 - m1); ++j){
					Y(my + j, ny + i) = U(m1 + j, n1 + i);
				}
			}
		}

		//! @brief 行列の各要素を上にm行分シフトする関数(下段の行はゼロになる)(引数渡し版)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		//! @param[in]	m	シフトする行数(デフォルト値 = 1)
		template<size_t P, size_t Q, typename R = double>
		static constexpr void shiftup(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y, const size_t m = 1){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			
			// 行を上にシフトする
			for(size_t j = 1; j <= P - m; ++j){
				setrow(Y, getrow(U, j + m), j);	// 横ベクトルを抽出して行に書き込み
			}
			
			// 残りの行をゼロにする
			for(size_t i = 1; i <= Q; ++i){
				for(size_t j = P - m + 1; j <= P; ++j){
					Y(j,i) = (R)0;
				}
			}
		}
		
		//! @brief 行列の各要素を上にm行分シフトする関数(下段の行はゼロになる)(戻り値渡し版)
		//! @param[in]	U	入力行列
		//! @param[in]	m	シフトする行数(デフォルト値 = 1)
		//! @return	Y	出力行列
		static constexpr ArcsMat<M,N,T> shiftup(const ArcsMat<M,N,T>& U, const size_t m = 1){
			ArcsMat<M,N,T> Y;
			shiftup(U, Y, m);
			return Y;
		}
		
		//! @brief 行列の各要素を下にm行分シフトする関数(上段の行はゼロになる)(引数渡し版)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		//! @param[in]	m	シフトする行数(デフォルト値 = 1)
		template<size_t P, size_t Q, typename R = double>
		static constexpr void shiftdown(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y, const size_t m = 1){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			
			// 行を下にシフトする
			for(size_t j = m + 1; j <= P; ++j){
				setrow(Y, getrow(U, j - m), j);	// 横ベクトルを抽出して行に書き込み
			}
			
			// 残りの行をゼロにする
			for(size_t i = 1; i <= Q; ++i){
				for(size_t j = 1; j <= m; ++j){
					Y(j,i) = (R)0;
				}
			}
		}
		
		//! @brief 行列の各要素を下にm行分シフトする関数(上段の行はゼロになる)(戻り値渡し版)
		//! @param[in]	U	入力行列
		//! @param[in]	m	シフトする行数(デフォルト値 = 1)
		//! @return	出力行列
		static constexpr ArcsMat<M,N,T> shiftdown(const ArcsMat<M,N,T>& U, const size_t m = 1){
			ArcsMat<M,N,T> Y;
			shiftdown(U, Y, m);
			return Y;
		}

		//! @brief 行列の各要素を左にn列分シフトする関数(右段の列はゼロになる)(引数渡し版)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		//! @param[in]	n	シフトする列数(デフォルト値 = 1)
		template<size_t P, size_t Q, typename R = double>
		static constexpr void shiftleft(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y, const size_t n = 1){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			
			// 列を左にシフトする
			for(size_t i = 1; i <= Q - n; ++i){
				setcolumn(Y, getcolumn(U, i + n), i);	// 縦ベクトルを抽出して列に書き込み
			}
			
			// 残りの列をゼロにする
			for(size_t i = Q - n + 1; i <= Q; ++i){
				for(size_t j = 1; j <= P; ++j){
					Y(j,i) = (R)0;
				}
			}
		}
		
		//! @brief 行列の各要素を左にn列分シフトする関数(右段の列はゼロになる)(戻り値渡し版)
		//! @param[in]	U	入力行列
		//! @param[in]	n	シフトする列数(デフォルト値 = 1)
		//! @return	Y	出力行列
		static constexpr ArcsMat<M,N,T> shiftleft(const ArcsMat<M,N,T>& U, const size_t n = 1){
			ArcsMat<M,N,T> Y;
			shiftleft(U, Y, n);
			return Y;
		}
		
		//! @brief 行列の各要素を右にn列分シフトする関数(左段の列はゼロになる)(引数渡し版)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		//! @param[in]	n	シフトする列数(デフォルト値 = 1)
		template<size_t P, size_t Q, typename R = double>
		static constexpr void shiftright(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y, const size_t n = 1){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			
			// 列を右にシフトする
			for(size_t i = n + 1; i <= Q; ++i){
				setcolumn(Y, getcolumn(U, i - n), i);	// 縦ベクトルを抽出して列に書き込み
			}
			
			// 残りの列をゼロにする
			for(size_t i = 1; i <= n; ++i){
				for(size_t j = 1; j <= P; ++j){
					Y(j,i) = (R)0;
				}
			}
		}
		
		//! @brief 行列の各要素を右にn列分シフトする関数(左段の列はゼロになる)(戻り値渡し版)
		//! @param[in]	U	入力行列
		//! @param[in]	n	シフトする列数(デフォルト値 = 1)
		//! @return	出力行列
		static constexpr ArcsMat<M,N,T> shiftright(const ArcsMat<M,N,T>& U, const size_t n = 1){
			ArcsMat<M,N,T> Y;
			shiftright(U, Y, n);
			return Y;
		}
		
		//! @brief 行列を縦方向にソートする関数 (引数渡し版)
		//! @tparam	ST	ソート方法： AMT_ASCENT = 昇順(デフォルト), AMT_DESCENT = 降順
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<ArcsMatrix::SortType ST = ArcsMatrix::SortType::AMT_ASCENT, size_t P, size_t Q, typename R = double>
		static constexpr void sortcolumn(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			std::array<T,M> v = {0};
			ArcsMat<M,1,T> w;
			for(size_t n = 1; n <= N; ++n){
				(ArcsMat<M,N,T>::getcolumn(U, n)).StoreArray(v);// std::arrayに読み込んでから、
				if constexpr(ST == ArcsMatrix::SortType::AMT_ASCENT){
					std::sort(v.begin(), v.end());				// 昇順ソートを実行
				}else if constexpr(ST == ArcsMatrix::SortType::AMT_DESCENT){
					std::sort(v.begin(), v.end(), [](const T& v1, const T& v2){ return v2 < v1;  } );	// 降順ソートを実行
				}else{
					arcs_assert(false);							// ここには来ない
				}
				w.LoadArray(v);									// ArcsMatに結果を読み込んで、
				ArcsMat<M,N,T>::setcolumn(Y, w, n);				// 結果を書き込み
			}
		}

		//! @brief 行列を縦方向にソートする関数 (戻り値返し版)
		//! @tparam	ST	ソート方法： AMT_ASCENT = 昇順(デフォルト), AMT_DESCENT = 降順
		//! @param[in]	U	入力行列
		//! @return	出力行列
		template<ArcsMatrix::SortType ST = ArcsMatrix::SortType::AMT_ASCENT>
		static constexpr ArcsMat<M,N,T> sortcolumn(const ArcsMat<M,N,T>& U){
			ArcsMat<M,N,T> Y;
			ArcsMat<M,N,T>::sortcolumn(U, Y);
			return Y;
		}

		//! @brief 行列を横方向にソートする関数 (引数渡し版)
		//! @tparam	ST	ソート方法： AMT_ASCENT = 昇順(デフォルト), AMT_DESCENT = 降順
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<ArcsMatrix::SortType ST = ArcsMatrix::SortType::AMT_ASCENT, size_t P, size_t Q, typename R = double>
		static constexpr void sortrow(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			ArcsMat<N,M,T> Ut, Yt;
			ArcsMat<M,N,T>::tp(U, Ut);						// 転置して、
			ArcsMat<N,M,T>::template sortcolumn<ST>(Ut, Yt);// 縦方向にソートして、
			ArcsMat<N,M,T>::tp(Yt, Y);						// また転置してもとに戻す
		}

		//! @brief 行列を横方向にソートする関数 (戻り値返し版)
		//! @tparam	ST	ソート方法： AMT_ASCENT = 昇順(デフォルト), AMT_DESCENT = 降順
		//! @param[in]	U	入力行列
		//! @return	出力行列
		template<ArcsMatrix::SortType ST = ArcsMatrix::SortType::AMT_ASCENT>
		static constexpr ArcsMat<M,N,T> sortrow(const ArcsMat<M,N,T>& U){
			ArcsMat<M,N,T> Y;
			ArcsMat<M,N,T>::template sortrow<ST>(U, Y);
			return Y;
		}

		//! @brief 行列を縦に連結する関数(引数渡し版)
		//! @tparam	P, Q, R, D, E, F	入力行列2と出力行列の高さ, 幅, 要素の型
		//! @param[in]	U1	入力行列1
		//! @param[in]	U2	入力行列2
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double, size_t D, size_t E, typename F = double>
		static constexpr void concatv(const ArcsMat<M,N,T>& U1, const ArcsMat<P,Q,R>& U2, ArcsMat<D,E,F>& Y){
			static_assert(N == Q, "ArcsMat: Size Error");				// 行列のサイズチェック
			static_assert(M + P == D, "ArcsMat: Output Size Error");	// 出力行列のサイズチェック
			static_assert(N == E, "ArcsMat: Output Size Error");		// 出力行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<F>, "ArcsMat: Type Error");	// 対応可能型チェック
			ArcsMat<D,E,F>::setsubmatrix(Y, U1, 1, 1);
			ArcsMat<D,E,F>::setsubmatrix(Y, U2, M + 1, 1);
		}
		
		//! @brief 行列を縦に連結する関数(戻り値渡し版)
		//! @tparam	P, Q, R	入力行列2の高さ, 幅, 要素の型
		//! @param[in]	U1	入力行列1
		//! @param[in]	U2	入力行列2
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr ArcsMat<M+P,N,T> concatv(const ArcsMat<M,N,T>& U1, const ArcsMat<P,Q,R>& U2){
			ArcsMat<M+P,N,T> Y;
			concatv(U1, U2, Y);
			return Y;
		}
		
		//! @brief 行列を横に連結する関数(引数渡し版)
		//! @tparam	P, Q, R, D, E, F	入力行列2と出力行列の高さ, 幅, 要素の型
		//! @param[in]	U1	入力行列1
		//! @param[in]	U2	入力行列2
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double, size_t D, size_t E, typename F = double>
		static constexpr void concath(const ArcsMat<M,N,T>& U1, const ArcsMat<P,Q,R>& U2, ArcsMat<D,E,F>& Y){
			static_assert(M == P, "ArcsMat: Size Error");				// 行列のサイズチェック
			static_assert(N + Q == E, "ArcsMat: Output Size Error");	// 出力行列のサイズチェック
			static_assert(M == D, "ArcsMat: Output Size Error");		// 出力行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<F>, "ArcsMat: Type Error");	// 対応可能型チェック
			ArcsMat<D,E,F>::setsubmatrix(Y, U1, 1, 1);
			ArcsMat<D,E,F>::setsubmatrix(Y, U2, 1, N + 1);
		}
		
		//! @brief 行列を横に連結する関数(戻り値渡し版)
		//! @tparam	P, Q, R	入力行列2の高さ, 幅, 要素の型
		//! @param[in]	U1	入力行列1
		//! @param[in]	U2	入力行列2
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr ArcsMat<M,N+Q,T> concath(const ArcsMat<M,N,T>& U1, const ArcsMat<P,Q,R>& U2){
			ArcsMat<M,N+Q,T> Y;
			concath(U1, U2, Y);
			return Y;
		}
		
		//! @brief 4つの行列を1つに連結する関数(引数渡し版)
		//! @tparam	P, Q, R, D, E, F, G, H, L, V, W, X	入力行列12-22と出力行列の高さ, 幅, 要素の型
		//! @param[in]	U11	入力行列11(左上)
		//! @param[in]	U12	入力行列12(右上)
		//! @param[in]	U21	入力行列11(左下)
		//! @param[in]	U22	入力行列12(右下)
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double, size_t D, size_t E, typename F = double, size_t G, size_t H, typename L = double, size_t V, size_t W, typename X = double>
		static constexpr void concat4(const ArcsMat<M,N,T>& U11, const ArcsMat<P,Q,R>& U12, const ArcsMat<D,E,F>& U21, const ArcsMat<G,H,L>& U22, ArcsMat<V,W,X>& Y){
			static_assert(M + D == P + G, "ArcsMat: Size Error");		// 行列のサイズチェック
			static_assert(N + Q == E + H, "ArcsMat: Size Error");		// 行列のサイズチェック
			static_assert(M + D == V, "ArcsMat: Output Size Error");	// 出力行列のサイズチェック
			static_assert(E + H == W, "ArcsMat: Output Size Error");	// 出力行列のサイズチェック
			static_assert(P + G == V, "ArcsMat: Output Size Error");	// 出力行列のサイズチェック
			static_assert(N + Q == W, "ArcsMat: Output Size Error");	// 出力行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<F>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<L>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<X>, "ArcsMat: Type Error");	// 対応可能型チェック
			if constexpr(N == E){
				// 縦長行列ベースで連結する場合
				static_assert(Q == H, "ArcsMat: Size Error");			// 行列のサイズチェック
				const ArcsMat<M+D,N> Ul = ArcsMat<M,N,T>::concatv(U11, U21);
				const ArcsMat<P+G,Q> Ur = ArcsMat<P,Q,R>::concatv(U12, U22);
				ArcsMat<M+D,N>::concath(Ul, Ur, Y);
			}else{
				// 横長行列ベースで連結する場合
				static_assert(M == P, "ArcsMat: Size Error");			// 行列のサイズチェック
				static_assert(D == G, "ArcsMat: Size Error");			// 行列のサイズチェック
				const ArcsMat<M,N+Q> Ut = ArcsMat<M,N,T>::concath(U11, U12);
				const ArcsMat<D,E+H> Ub = ArcsMat<D,E,F>::concath(U21, U22);
				ArcsMat<M,N+Q>::concatv(Ut, Ub, Y);
			}
		}
		
		//! @brief 4つの行列を1つに連結する関数(戻り値渡し版)
		//! @tparam	P, Q, R, D, E, F, G, H, L	入力行列12-22と出力行列の高さ, 幅, 要素の型
		//! @param[in]	U11	入力行列11(左上)
		//! @param[in]	U12	入力行列12(右上)
		//! @param[in]	U21	入力行列11(左下)
		//! @param[in]	U22	入力行列12(右下)
		//! @return	出力行列
		template<size_t P, size_t Q, typename R = double, size_t D, size_t E, typename F = double, size_t G, size_t H, typename L = double>
		static constexpr ArcsMat<M+D,N+Q,T> concat4(const ArcsMat<M,N,T>& U11, const ArcsMat<P,Q,R>& U12, const ArcsMat<D,E,F>& U21, const ArcsMat<G,H,L>& U22){
			ArcsMat<M+D,N+Q,T> Y;
			concat4(U11, U12, U21, U22, Y);
			return Y;
		}
		
		//! @brief 縦ベクトルの各要素を対角要素に持つ正方行列を生成する関数(引数渡し版)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	u	入力ベクトル
		//! @param[in]	Y	出力行列
		//! @return	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void diag(const ArcsMat<M,N,T>& u, ArcsMat<P,Q,R>& Y){
			static_assert(N == 1, "ArcsMat: Vector Error");	// 縦ベクトルチェック
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(P == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			for(size_t j = 1; j <= M; ++j) Y(j,j) = static_cast<R>( u(j,1) );
		}
		
		//! @brief 縦ベクトルの各要素を対角要素に持つ正方行列を生成する関数(戻り値渡し版)
		//! @param[in]	u	入力ベクトル
		//! @return	出力行列
		static constexpr ArcsMat<M,M,T> diag(const ArcsMat<M,N,T>& u){
			ArcsMat<M,M,T> Y;
			diag(u, Y);
			return Y;
		}
		
		//! @brief 行列の対角要素を縦ベクトルとして取得する関数(引数渡し版)
		//! @tparam	P, Q, R, L	出力行列の高さ, 幅, 要素の型, 対角要素の数
		//! @param[in]	U	入力行列
		//! @param[in]	k	k番目の対角, k=0で主対角、K<0で主対角より下、0<kで主対角より上 (デフォルト値 = 0)
		//! @param[out]	y	出力ベクトル
		template<size_t P, size_t Q, typename R = double, size_t L = std::min(M,N)>
		static constexpr void getdiag(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& y, const ssize_t k = 0){
			static_assert(Q == 1, "ArcsMat: Vector Error");	// 縦ベクトルチェック
			static_assert(P == L, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック

			// オフセット条件設定
			size_t m = 1, j = 0, i = 0;
			if(0 <= k){
				// 主対角 or 主対角より上を取る場合
				j = 1;
				i = static_cast<size_t>( 1 + k );
			}else if(k < 0){
				// 主対角より下の対角を取る場合
				j = static_cast<size_t>( 1 - k );
				i = 1;
			}else{
				arcs_assert(false);	// ここには来ない
			}

			// 対角要素を抽出
			while((j <= M) && (i <= N)){
				y(m,1) = static_cast<R>( U(j,i) );
				++m;
				++j;
				++i;
			}
		}
				
		//! @brief 行列の対角要素を縦ベクトルとして取得する関数(戻り値渡し版)
		//! @tparam	L	対角要素の数
		//! @param[in]	U	入力行列
		//! @param[in]	k	k番目の対角, k=0で主対角、K<0で主対角より下、0<kで主対角より上 (デフォルト値 = 0)
		//! @return	出力ベクトル
		template<size_t L = std::min(M,N)>
		static constexpr ArcsMat<L,1,T> getdiag(const ArcsMat<M,N,T>& U, const ssize_t k = 0){
			ArcsMat<L,1,T> y;
			getdiag(U, y, k);
			return y;
		}
		
		//! @brief 行列のトレースを返す関数(戻り値渡し版のみ)
		//! @param[in]	U	入力行列
		//! @return	結果
		static constexpr T trace(const ArcsMat<M,N,T>& U){
			const ArcsMat<1,1,T> y = ArcsMat<M,1,T>::sumcolumn(ArcsMat<M,N,T>::getdiag(U));
			return y[1];
		}
		
		//! @brief 行列の対角要素の総積を返す関数(戻り値渡し版のみ)
		//! @tparam	L	対角要素の数
		//! @param[in]	U	入力行列
		//! @return	結果
		template<size_t L = std::min(M,N)>
		static constexpr T multdiag(const ArcsMat<M,N,T>& U){
			T y = 1;
			for(size_t j = 1; j <= L; ++j) y *= U(j,j);
			return y;
		}
		
		//! @brief 行列要素の最大値の要素番号を返す関数(タプル返し版のみ)
		//! @param[in]	U	入力行列
		//! @return	結果(タプルで返す)
		static constexpr std::tuple<size_t,size_t> maxidx(const ArcsMat<M,N,T>& U){
			size_t k = 1, l = 1;
			for(size_t i = 1; i <= N; ++i){		// 横方向探索
				for(size_t j = 1; j <= M; ++j){	// 縦方向探索
					if constexpr(ArcsMatrix::IsComplex<T>){
						// 複素数版のとき
						if( std::abs(U(k,l)) < std::abs(U(j,i)) ){
							k = j;
							l = i;
						}
					}else{
						// 実数版のとき
						if( U(k,l) < U(j,i) ){
							k = j;
							l = i;
						}
					}
				}
			}
			return {k, l};
		}
		
		//! @brief 行列要素の最大値を返す関数(戻り値渡し版のみ)
		//! @param[in]	U	入力行列
		//! @return	結果
		static constexpr T max(const ArcsMat<M,N,T>& U){
			const auto ji = maxidx(U);
			return U(std::get<0>(ji), std::get<1>(ji));
		}
		
		//! @brief 行列要素の最小値の要素番号を返す関数(タプル返し版のみ)
		//! @param[in]	U	入力行列
		//! @return	結果(タプルで返す)
		static constexpr std::tuple<size_t,size_t> minidx(const ArcsMat<M,N,T>& U){
			size_t k = 1, l = 1;
			for(size_t i = 1; i <= N; ++i){		// 横方向探索
				for(size_t j = 1; j <= M; ++j){	// 縦方向探索
					if constexpr(ArcsMatrix::IsComplex<T>){
						// 複素数版のとき
						if( std::abs(U(k,l)) > std::abs(U(j,i)) ){
							k = j;
							l = i;
						}
					}else{
						// 実数版のとき
						if( U(k,l) > U(j,i) ){
							k = j;
							l = i;
						}
					}
				}
			}
			return {k, l};
		}
		
		//! @brief 行列要素の最小値を返す関数(戻り値渡し版のみ)
		//! @param[in]	U	入力行列
		//! @return	結果
		static constexpr T min(const ArcsMat<M,N,T>& U){
			const auto ji = minidx(U);
			return U(std::get<0>(ji), std::get<1>(ji));
		}

		//! @brief 行列要素の総和を返す関数(戻り値渡し版のみ)
		//! @param[in]	U	入力行列
		//! @return	結果
		static constexpr T sum(const ArcsMat<M,N,T>& U){
			ArcsMat<1,1,T> y = ArcsMat<1,N,T>::sumrow( ArcsMat<M,N,T>::sumcolumn(U) );
			return y[1];
		}

		//! @brief 行列要素の指数関数を計算する関数(引数渡し版)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void exp(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_floating_point_v<T>, "ArcsMat: Type Error (Floating Point)");	// 型チェック
			static_assert(std::is_floating_point_v<R>, "ArcsMat: Type Error (Floating Point)");	// 型チェック
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) Y(j,i) = static_cast<R>( std::exp( U(j,i) ) );
			}
		}

		//! @brief 行列要素の指数関数を計算する関数(戻り値渡し版)
		//! @param[in]	U	入力行列
		//! @return	Y	出力行列
		static constexpr ArcsMat<M,N,T> exp(const ArcsMat<M,N,T>& U){
			ArcsMat<M,N,T> Y;
			ArcsMat<M,N,T>::exp(U, Y);
			return Y;
		}

		//! @brief 行列要素の対数関数(底e版)を計算する関数(引数渡し版)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void log(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_floating_point_v<T>, "ArcsMat: Type Error (Floating Point)");	// 型チェック
			static_assert(std::is_floating_point_v<R>, "ArcsMat: Type Error (Floating Point)");	// 型チェック
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) Y(j,i) = static_cast<R>( std::log( U(j,i) ) );
			}
		}

		//! @brief 行列要素の対数関数(底e版)を計算する関数(戻り値渡し版)
		//! @param[in]	U	入力行列
		//! @return	Y	出力行列
		static constexpr ArcsMat<M,N,T> log(const ArcsMat<M,N,T>& U){
			ArcsMat<M,N,T> Y;
			ArcsMat<M,N,T>::log(U, Y);
			return Y;
		}

		//! @brief 行列要素の対数関数(底10版)を計算する関数(引数渡し版)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void log10(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_floating_point_v<T>, "ArcsMat: Type Error (Floating Point)");	// 型チェック
			static_assert(std::is_floating_point_v<R>, "ArcsMat: Type Error (Floating Point)");	// 型チェック
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) Y(j,i) = static_cast<R>( std::log10( U(j,i) ) );
			}
		}

		//! @brief 行列要素の対数関数(底10版)を計算する関数(戻り値渡し版)
		//! @param[in]	U	入力行列
		//! @return	Y	出力行列
		static constexpr ArcsMat<M,N,T> log10(const ArcsMat<M,N,T>& U){
			ArcsMat<M,N,T> Y;
			ArcsMat<M,N,T>::log10(U, Y);
			return Y;
		}

		//! @brief 行列要素の正弦関数を計算する関数(引数渡し版)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void sin(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_floating_point_v<T>, "ArcsMat: Type Error (Floating Point)");	// 型チェック
			static_assert(std::is_floating_point_v<R>, "ArcsMat: Type Error (Floating Point)");	// 型チェック
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) Y(j,i) = static_cast<R>( std::sin( U(j,i) ) );
			}
		}

		//! @brief 行列要素の正弦関数を計算する関数(戻り値渡し版)
		//! @param[in]	U	入力行列
		//! @return	Y	出力行列
		static constexpr ArcsMat<M,N,T> sin(const ArcsMat<M,N,T>& U){
			ArcsMat<M,N,T> Y;
			ArcsMat<M,N,T>::sin(U, Y);
			return Y;
		}

		//! @brief 行列要素の余弦関数を計算する関数(引数渡し版)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void cos(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_floating_point_v<T>, "ArcsMat: Type Error (Floating Point)");	// 型チェック
			static_assert(std::is_floating_point_v<R>, "ArcsMat: Type Error (Floating Point)");	// 型チェック
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) Y(j,i) = static_cast<R>( std::cos( U(j,i) ) );
			}
		}

		//! @brief 行列要素の余弦を計算する関数(戻り値渡し版)
		//! @param[in]	U	入力行列
		//! @return	Y	出力行列
		static constexpr ArcsMat<M,N,T> cos(const ArcsMat<M,N,T>& U){
			ArcsMat<M,N,T> Y;
			ArcsMat<M,N,T>::cos(U, Y);
			return Y;
		}

		//! @brief 行列要素の正接関数を計算する関数(引数渡し版)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void tan(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_floating_point_v<T>, "ArcsMat: Type Error (Floating Point)");	// 型チェック
			static_assert(std::is_floating_point_v<R>, "ArcsMat: Type Error (Floating Point)");	// 型チェック
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) Y(j,i) = static_cast<R>( std::tan( U(j,i) ) );
			}
		}

		//! @brief 行列要素の正接関数を計算する関数(戻り値渡し版)
		//! @param[in]	U	入力行列
		//! @return	Y	出力行列
		static constexpr ArcsMat<M,N,T> tan(const ArcsMat<M,N,T>& U){
			ArcsMat<M,N,T> Y;
			ArcsMat<M,N,T>::tan(U, Y);
			return Y;
		}

		//! @brief 行列要素の双曲線正接関数を計算する関数(引数渡し版)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void tanh(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_floating_point_v<T>, "ArcsMat: Type Error (Floating Point)");	// 型チェック
			static_assert(std::is_floating_point_v<R>, "ArcsMat: Type Error (Floating Point)");	// 型チェック
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) Y(j,i) = static_cast<R>( std::tanh( U(j,i) ) );
			}
		}

		//! @brief 行列要素の双曲線正接関数を計算する関数(戻り値渡し版)
		//! @param[in]	U	入力行列
		//! @return	Y	出力行列
		static constexpr ArcsMat<M,N,T> tanh(const ArcsMat<M,N,T>& U){
			ArcsMat<M,N,T> Y;
			ArcsMat<M,N,T>::tanh(U, Y);
			return Y;
		}
		
		//! @brief 行列要素の平方根を計算する関数(引数渡し版)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void sqrt(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) Y(j,i) = static_cast<R>( std::sqrt( U(j,i) ) );
			}
		}
		
		//! @brief 行列要素の平方根を計算する関数(戻り値渡し版)
		//! @param[in]	U	入力行列
		//! @return	Y	出力行列
		static constexpr ArcsMat<M,N,T> sqrt(const ArcsMat<M,N,T>& U){
			ArcsMat<M,N,T> Y;
			ArcsMat<M,N,T>::sqrt(U, Y);
			return Y;
		}
		
		//! @brief 行列要素の符号関数を計算する関数(引数渡し版)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void sign(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) Y(j,i) = static_cast<R>( ArcsMatrix::sgn( U(j,i) ) );
			}
		}
		
		//! @brief 行列要素の符号関数を計算する関数(戻り値渡し版)
		//! @param[in]	U	入力行列
		//! @return	Y	出力行列
		static constexpr ArcsMat<M,N,T> sign(const ArcsMat<M,N,T>& U){
			ArcsMat<M,N,T> Y;
			ArcsMat<M,N,T>::sign(U, Y);
			return Y;
		}
		
		//! @brief 行列要素の絶対値を計算する関数(引数渡し版)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void abs(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) Y(j,i) = static_cast<R>( std::abs( U(j,i) ) );
			}
		}

		//! @brief 行列要素の絶対値を計算する関数(戻り値渡し版)
		//! @tparam	R	出力行列の要素の型
		//! @param[in]	U	入力行列
		//! @return	Y	出力行列
		template<typename R = double>
		static constexpr ArcsMat<M,N,R> abs(const ArcsMat<M,N,T>& U){
			ArcsMat<M,N,R> Y;
			ArcsMat<M,N,T>::abs(U, Y);
			return Y;
		}
		
		//! @brief 複素数行列要素の偏角を計算する関数(引数渡し版)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void arg(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsComplex<T>, "ArcsMat: Type Error (Need Complex)");
			static_assert(ArcsMatrix::IsIntFloat<R>, "ArcsMat: Type Error (Need Integer or Float)");
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) Y(j,i) = static_cast<R>( std::arg( U(j,i) ) );
			}
		}
		
		//! @brief 複素数行列要素の実数部を取得する関数(引数渡し版)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void real(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsComplex<T>, "ArcsMat: Type Error (Need Complex)");
			static_assert(ArcsMatrix::IsIntFloat<R>, "ArcsMat: Type Error (Need Integer or Float)");
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) Y(j,i) = static_cast<R>( std::real( U(j,i) ) );
			}
		}

		//! @brief 複素数行列要素の実数部を計算する関数(戻り値渡し版, 複素数版特殊化)
		//! @tparam	R	入力行列の複素数要素の型
		//! @param[in]	U	入力行列
		//! @return	Y	出力行列
		template<typename R = double>
		static constexpr ArcsMat<M,N,R> real(const ArcsMat<M,N,std::complex<R>>& U){
			ArcsMat<M,N,R> Y;							// 実数返し用
			ArcsMat<M,N,std::complex<R>>::real(U, Y);	// 入力は複素数、出力は実数
			return Y;	// 実数で返す
		}
		
		//! @brief 複素数行列要素の虚数部を取得する関数(引数渡し版)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void imag(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsComplex<T>, "ArcsMat: Type Error (Need Complex)");
			static_assert(ArcsMatrix::IsIntFloat<R>, "ArcsMat: Type Error (Need Integer or Float)");
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) Y(j,i) = static_cast<R>( std::imag( U(j,i) ) );
			}
		}
		
		//! @brief 複素数行列要素の虚数部を計算する関数(戻り値渡し版, 複素数版特殊化)
		//! @tparam	R	入力行列の複素数要素の型
		//! @param[in]	U	入力行列
		//! @return	Y	出力行列
		template<typename R = double>
		static constexpr ArcsMat<M,N,R> imag(const ArcsMat<M,N,std::complex<R>>& U){
			ArcsMat<M,N,R> Y;							// 実数返し用
			ArcsMat<M,N,std::complex<R>>::imag(U, Y);	// 入力は複素数、出力は実数
			return Y;	// 実数で返す
		}
		
		//! @brief 複素数行列要素の複素共役を取得する関数(引数渡し版)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = std::complex<double>>
		static constexpr void conj(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsComplex<T>, "ArcsMat: Type Error (Need Complex)");
			static_assert(ArcsMatrix::IsComplex<R>, "ArcsMat: Type Error (Need Complex)");
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) Y(j,i) = static_cast<R>( std::conj( U(j,i) ) );
			}
		}
		
		//! @brief 複素数行列要素の複素共役を取得する関数(戻り値渡し版)
		//! @param[in]	U	入力行列
		//! @return	出力行列
		static constexpr ArcsMat<M,N,T> conj(const ArcsMat<M,N,T>& U){
			ArcsMat<M,N,T> ret;
			ArcsMat<M,N,T>::conj(U, ret);
			return ret;
		}
		
		//! @brief 転置行列を返す関数 (引数渡し版)
		//! @tparam	P, Q, R	出力の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void tp(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(M == Q, "ArcsMat: Transpose Size Error");	// 転置サイズチェック
			static_assert(N == P, "ArcsMat: Transpose Size Error");	// 転置サイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) Y(i,j) = static_cast<R>( U(j,i) );
			}
		}
		
		//! @brief 転置行列を返す関数 (戻り値渡し版)
		//! @param[in]	U	入力行列
		//! @return	Y	出力行列
		static constexpr ArcsMat<N,M,T> tp(const ArcsMat<M,N,T>& U){
			ArcsMat<N,M,T> ret;
			ArcsMat<M,N,T>::tp(U, ret);
			return ret;
		}
		
		//! @brief エルミート転置行列を返す関数 (引数渡し版)
		//! @tparam	P, Q, R	出力の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = std::complex<double>>
		static constexpr void Htp(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(ArcsMatrix::IsComplex<T>, "ArcsMat: Type Error (Need Complex)");
			static_assert(ArcsMatrix::IsComplex<R>, "ArcsMat: Type Error (Need Complex)");
			ArcsMat<M,N,T>::tp(ArcsMat<M,N,T>::conj(U), Y);	// 複素共役して転置
		}
		
		//! @brief エルミート転置行列を返す関数 (戻り値渡し版)
		//! @param[in]	U	入力行列
		//! @return	出力行列
		static constexpr ArcsMat<N,M,T> Htp(const ArcsMat<M,N,T>& U){
			ArcsMat<N,M,T> ret;
			ArcsMat<M,N,T>::Htp(U, ret);
			return ret;
		}

		//! @brief 転置演算子
		//! @return 結果
		constexpr ArcsMat<N,M,T> operator~(void) const{
			ArcsMat<N,M,T> ret;
			if constexpr(ArcsMatrix::IsComplex<T>){
				// 複素数型の場合
				ret = Htp(*this);	// エルミート転置して返す
			}else{
				// それ以外の型の場合
				ret = tp(*this);	// 普通に転置して返す
			}
			return ret;
		}

		//! @brief 行列のノルムを返す関数(戻り値渡し版のみ)
		//! @tparam	NRM	ノルムのタイプ (デフォルト値 = AMT_L2)
		//! @param[in]	U	入力行列
		//! @return	結果
		template<ArcsMatrix::NormType NRM = ArcsMatrix::NormType::AMT_L2, typename R = double>
		static constexpr R norm(const ArcsMat<M,N,T>& U){
			R ret = 0;
			if constexpr(NRM == ArcsMatrix::NormType::AMT_L2){
				// ユークリッドノルム(2-ノルム)が指定されたとき
				if constexpr(M == 1 || N == 1){
					// ベクトル版
					ArcsMat<M,N,R> v = ArcsMat<M,N,T>::abs(U);
					ret = std::sqrt( ArcsMat<M,N,R>::sum( v & v ) );
				}else{
					// 行列版
					const auto [W, S, V] = ArcsMat<M,N,T>::SVD(U);
					ret = std::abs( ArcsMat<M,N,T>::max(S) );
				}
			}else if constexpr(NRM == ArcsMatrix::NormType::AMT_L1){
				// 絶対値ノルム(1-ノルム)が指定されたとき
				if constexpr(M == 1 || N == 1){
					// ベクトル版
					ret = ArcsMat<M,N,R>::sum( ArcsMat<M,N,T>::abs(U) );
				}else{
					// 行列版
					ret = ArcsMat<1,N,R>::max( ArcsMat<M,N,R>::sumcolumn( ArcsMat<M,N,T>::abs(U) ) );
				}
			}else if constexpr(NRM == ArcsMatrix::NormType::AMT_LINF){
				// 無限大ノルム(最大値ノルム)が指定されたとき
				if constexpr(M == 1 || N == 1){
					// ベクトル版
					ret = ArcsMat<M,N,R>::max( ArcsMat<M,N,T>::abs(U) );
				}else{
					// 行列版
					ret = ArcsMat<M,1,R>::max( ArcsMat<M,N,R>::sumrow( ArcsMat<M,N,T>::abs(U) ) );
				}
			}else{
				arcs_assert(false);	// ここには来ない
			}
			return ret;
		}

		//! @brief n列目を左端として右上の上三角部分のみを返す関数(下三角部分はゼロ)(引数渡し版)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		//! @param[in]	n	切り出す左端位置 n列目(デフォルト値 = 1)
		template<size_t P, size_t Q, typename R = double>
		static constexpr void gettriup(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y, const size_t n = 1){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック

			for(size_t j = 1; j <= M; ++j){
				for(size_t i = j + n - 1; i <= N; ++i){
					Y(j,i) = static_cast<R>( U(j,i) );
				}
			}
		}

		//! @brief n列目を左端として右上の上三角部分のみを返す関数(下三角部分はゼロ)(戻り値渡し版)
		//! @param[in]	U	入力行列
		//! @param[in]	n	切り出す左端位置 n列目(デフォルト値 = 1)
		//! @return	出力行列
		static constexpr ArcsMat<M,N,T> gettriup(const ArcsMat<M,N,T>& U, const size_t n = 1){
			ArcsMat<M,N,T> Y;
			gettriup(U, Y, n);
			return Y;
		}

		//! @brief m行目を上端として左下の下三角部分のみを返す関数(上三角部分はゼロ)(引数渡し版)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		//! @param[in]	m	切り出す上端位置 m行目(デフォルト値 = 1)
		template<size_t P, size_t Q, typename R = double>
		static constexpr void gettrilo(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y, const size_t m = 1){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック

			for(size_t i = 1; i <= N; ++i){
				for(size_t j = i + m - 1; j <= M; ++j){
					Y(j,i) = static_cast<R>( U(j,i) );
				}
			}
		}

		//! @brief m行目を上端として左下の下三角部分のみを返す関数(上三角部分はゼロ)(戻り値渡し版)
		//! @param[in]	U	入力行列
		//! @param[in]	m	切り出す上端位置 m行目(デフォルト値 = 1)
		//! @return	出力行列
		static constexpr ArcsMat<M,N,T> gettrilo(const ArcsMat<M,N,T>& U, const size_t m = 1){
			ArcsMat<M,N,T> Y;
			gettrilo(U, Y, m);
			return Y;
		}

		//! @brief LU分解の結果と置換行列を返す関数(引数渡し版)
		//! @tparam	ML, NL, TL, MU, NU, TU, MP, NP, TP	L,U,P行列の高さ, 幅, 要素の型
		//! @param[in]	A	入力行列
		//! @param[out]	L	下三角行列
		//! @param[out]	U	上三角行列
		//! @param[out]	P	置換行列
		template<size_t ML, size_t NL, typename TL = double, size_t MU, size_t NU, typename TU = double, size_t MP, size_t NP, typename TP = double>
		static constexpr void LUP(const ArcsMat<M,N,T>& A, ArcsMat<ML,NL,TL>& L, ArcsMat<MU,NU,TU>& U, ArcsMat<MP,NP,TP>& P){
			static_assert(M == N, "ArcsMat: Size Error");	// 正方行列チェック
			static_assert(ML == NL, "ArcsMat: Size Error");	// 正方行列チェック
			static_assert(MU == NU, "ArcsMat: Size Error");	// 正方行列チェック
			static_assert(MP == NP, "ArcsMat: Size Error");	// 正方行列チェック
			static_assert(M == ML, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == NL, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(M == MU, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == NU, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(M == MP, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == NP, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<TL>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TU>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TP>, "ArcsMat: Type Error");	// 対応可能型チェック
			
			// 中間変数
			ArcsMat<M,N,T> X = A;		// 行列操作用にコピー
			size_t perm_count = 0;		// 入れ替えカウンタ
			double max_buff = 0;		// 最大値バッファ
			P = ArcsMat<M,N,T>::eye();	// 置換行列の準備

			// 列ごとに処理を実施
			size_t k = 0;
			for(size_t i = 1; i <= N; ++i){
				// i列目の中での最大値を探す
				k = i;									// 対角要素の縦位置で初期化
				max_buff = static_cast<double>( std::abs(X(i,i)) );			// 対角要素の値で初期化
				for(size_t j = i + 1; j <= M; ++j){
					if(max_buff < static_cast<double>( std::abs(X(j,i)) )){	// スキャン中における最大値か判定
						k = j;							// 最大値の候補の縦位置
						max_buff = static_cast<double>( std::abs(X(j,i)) );	// 最大値の候補であればその値を記憶
					}
				}

				// 対角要素が最大値でなければ、対角要素が最大となるように行丸ごと入れ替え
				if(k != i){
					ArcsMat<M,N,T>::swaprow(X, i, k);	// 対角要素の行と最大値の行を入れ替え
					ArcsMat<M,N,T>::swaprow(P, i, k);	// 置換行列も同じ様に入れ替え
					perm_count++;						// 入れ替えカウンタ
				}

				// LU分解のコア部分
				if( std::abs( X(i,i) ) < ArcsMatrix::EPSILON ){
					// 対角要素が零なら，i列目においてはLU分解は完了
				}else{
					for(size_t j = i + 1; j <= M; ++j){
						X(j,i) /= X(i,i);					// 対角要素で除算
						for(size_t l = i + 1; l <= N; ++l){
							X(j,l) -= X(j,i)*X(i,l);
						}
					}
				}
			}
			
			// 下三角行列と上三角行列に分離する
			ArcsMat<M,N,T>::gettrilo(X, L);	// 下三角のみを抽出
			ArcsMat<M,N,T>::gettriup(X, U);	// 上三角のみを抽出
			for(size_t j = 1; j <= M; ++j) L(j,j) = 1;	// 下三角行列の対角要素はすべて1
			
			// 入れ替え回数の判定と判定結果の保持
			if(perm_count % 2 == 0){
				P.Status = ArcsMatrix::MatStatus::AMT_LU_EVEN;	// 偶数のとき
			}else{
				P.Status = ArcsMatrix::MatStatus::AMT_LU_ODD;	// 奇数のとき
			}
		}

		//! @brief LU分解の結果と置換行列を返す関数(タプル返し版)
		//! @param[in]	A	入力行列
		//! @return	(L, U, P)	(下三角行列, 上三角行列, 置換行列)のタプル
		static constexpr std::tuple<ArcsMat<M,N,T>, ArcsMat<M,N,T>, ArcsMat<M,N,T>> LUP(const ArcsMat<M,N,T>& A){
			ArcsMat<M,N,T> L, U, P;
			ArcsMat<M,N,T>::LUP(A, L, U, P);
			return {L, U, P};
		}
		
		//! @brief LU分解の結果のみ返す関数(引数渡し版)
		//! @tparam	ML, NL, TL, MU, NU, TU	L,U行列の高さ, 幅, 要素の型
		//! @param[in]	A	入力行列
		//! @param[out]	L	下三角行列
		//! @param[out]	U	上三角行列
		template<size_t ML, size_t NL, typename TL = double, size_t MU, size_t NU, typename TU = double>
		static constexpr void LU(const ArcsMat<M,N,T>& A, ArcsMat<ML,NL,TL>& L, ArcsMat<MU,NU,TU>& U){
			ArcsMat<M,N,T> P;
			ArcsMat<M,N,T>::LUP(A, L, U, P);
			L = tp(P)*L;
		}

		//! @brief LU分解の結果のみ返す関数(タプル返し版)
		//! @param[in]	A	入力行列
		//! @return	(L, U)	(下三角行列, 上三角行列)のタプル
		static constexpr std::tuple<ArcsMat<M,N,T>, ArcsMat<M,N,T>> LU(const ArcsMat<M,N,T>& A){
			ArcsMat<M,N,T> L, U;
			ArcsMat<M,N,T>::LU(A, L, U);
			return {L, U};
		}
		
		//! @brief 行列式の値を返す関数(戻り値返し版のみ)
		//! @param[in]	A	入力行列
		//! @return	結果
		static constexpr T det(const ArcsMat<M,N,T>& A){
			static_assert(M == N, "ArcsMat: Size Error");		// 正方行列チェック

			// |A| = |L||U| でしかも |L|と|U|は対角要素の総積に等しく，さらにLの対角要素は1なので|L|は省略可。
			// LU分解の際に内部で並べ替える毎に符号が変わるので、並べ替え回数によって符号反転をする。
			const auto [L, U, P] = LUP(A);	// LU分解
			if(P.Status == ArcsMatrix::MatStatus::AMT_LU_ODD){	// LU分解の符号判定
				return -ArcsMat<M,N,T>::multdiag(U);	// 奇数のとき
			}else if(P.Status == ArcsMatrix::MatStatus::AMT_LU_EVEN){
				return  ArcsMat<M,N,T>::multdiag(U);	// 偶数のとき
			}else{
				arcs_assert(false);	// ここには来ない
			}
		}
		
		//! @brief Householder行列を生成する関数(引数渡し版)
		//! @tparam	MH, NH, TH	H行列の高さ, 幅, 要素の型
		//! @param[in]	v	入力縦ベクトル
		//! @param[out]	H	ハウスホルダー行列
		//! @param[in]	k	次元 (デフォルト値 = 1)
		template<size_t MH, size_t NH, typename TH = double>
		static constexpr void Householder(const ArcsMat<M,N,T>& v, ArcsMat<MH,NH,TH>& H, const size_t k = 1){
			static_assert(N == 1, "ArcsMat: Vector Error");	// 縦ベクトルチェック
			static_assert(M == MH,  "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(MH == NH, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TH>, "ArcsMat: Type Error");	// 対応可能型チェック
			arcs_assert(k <= M);	// 次元チェック
			
			const auto I = ArcsMat<MH,NH,TH>::eye();	// [C++20移行時にconstexprに改修]
			auto vHv = ~v*v;
			if( std::abs(vHv[1]) != 0 ){
				H = I - 2.0*v*~v/vHv[1];
			}else{
				H = I - 2.0*v*~v;
			}
			H = ArcsMat<M,M,T>::shiftdown(H, k);
			H = ArcsMat<M,M,T>::shiftright(H, k);
			for(size_t i = 1; i <= k; ++i) H(i,i) = 1;
		}
		
		//! @brief Householder行列を生成する関数(戻り値返し版)
		//! @param[in]	v	入力縦ベクトル
		//! @param[in]	k	次元 (デフォルト値 = 1)
		//! @return	ハウスホルダー行列
		static constexpr ArcsMat<M,M,T> Householder(const ArcsMat<M,N,T>& v, const size_t k = 1){
			ArcsMat<M,M,T> H;
			ArcsMat<M,N,T>::Householder(v, H, k);
			return H;
		}
		
		//! @brief QR分解(引数渡し版)
		//! 注意：複素数で縦長行列の場合ではMATLABとは異なる解を出力する
		//! @tparam	MQ, NQ, TQ, MQ, NQ, TQ	Q,R行列の高さ, 幅, 要素の型
		//! @param[in]	A	入力行列
		//! @param[out]	Q	ユニタリ行列 Q行列
		//! @param[out] R	上三角行列 R行列
		template<size_t MQ, size_t NQ, typename TQ = double, size_t MR, size_t NR, typename TR = double>
		static constexpr void QR(const ArcsMat<M,N,T>& A, ArcsMat<MQ,NQ,TQ>& Q, ArcsMat<MR,NR,TR>& R){
			static_assert(MQ == NQ, "ArcsMat: Size Error");	// 正方行列チェック
			static_assert(M == MQ, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(M == MR, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == NR, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<TQ>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TR>, "ArcsMat: Type Error");	// 対応可能型チェック

			// 事前準備
			constexpr size_t K = std::min(M,N);	// 縦長か横長か、短い方をKとする
			constexpr ArcsMat<M,1,T> e = {1};	// 単位ベクトルを生成
			ArcsMat<M,1,T> a, v;				// "Householder Reflections"のベクトル
			ArcsMat<M,M,T> H;					// Householder行列
			Q = ArcsMat<MQ,NQ,TQ>::eye();
			R = A;

			// "Householder Reflections"を使ったQR分解アルゴリズム(複素数・実数両対応版)
			for(size_t k = 1; k <= K; ++k){
				a = ArcsMat<M,N,T>::getcolumn(R, k);
				a = ArcsMat<M,1,T>::shiftup(a, k - 1);
				
				// 単位行列部分の生成を実数と複素数で変える
				if constexpr(ArcsMatrix::IsComplex<T>){
					// 複素数の場合
					v = a + std::exp( std::complex( 0.0, std::arg(a[1])) )*std::sqrt( (~a*a)[1] )*e;
				}else{
					// 実数の場合
					v = a + ArcsMatrix::sgn(a[1])*std::sqrt( (~a*a)[1] )*e;
				}
				
				ArcsMat<M,1,T>::Householder(v, H, k - 1);	// ハウスホルダー行列を生成
				R = H*R;
				Q = Q*~H;
			}
			
			// Aが縦長を除き、且つ複素数の場合にはRの対角項を実数に変換
			/* (複素Schur分解で不具合が出たのでPending)
			if constexpr( (M <= N) && ArcsMatrix::IsComplex<T>){
				const auto l  = ArcsMat<K,1,T>::ones();	// [C++20移行時にconstexprに改修]
				const auto d  = ArcsMat<K,1,T>::sign( ArcsMat<MR,NR,T>::getdiag(R) );
				const auto di = l % d;
				const auto D  = ArcsMat<K,1,T>::diag(d);
				const auto Di = ArcsMat<K,1,T>::diag(di);
				Q = Q*D;
				R = Di*R;
			}
			*/
			R.ZeroingTriLo();	// 主対角より下の下三角のゼロイング
		}

		//! @brief QR分解(タプル返し版)
		//! @param[in]	A	入力行列
		//! @return	(Q, R)	(直交行列, 上三角行列)のタプル
		static constexpr std::tuple<ArcsMat<M,M,T>, ArcsMat<M,N,T>> QR(const ArcsMat<M,N,T>& A){
			ArcsMat<M,M,T> Q;
			ArcsMat<M,N,T> R;
			ArcsMat<M,N,T>::QR(A, Q, R);
			return {Q, R};
		}

		//! @brief SVD特異値分解(引数渡し版)
		//! 注意：複素数で非正方行列の場合ではMATLABとは異なる解を出力する
		//! @tparam LoopMax	ループ打ち切り最大回数 デフォルト値 = 100*max(M,N)
		//! @tparam	MU, NU, TU, MS, NS, TS, MV, NV, TV	U,Σ,V行列の高さ, 幅, 要素の型
		//! @param[in]	A	入力行列
		//! @param[out]	U	U行列
		//! @param[out]	S	Σ行列
		//! @param[out]	V	V行列
		template<
			size_t LoopMax = 100*std::max(M,N),
			size_t MU, size_t NU, typename TU = double, size_t MS, size_t NS, typename TS = double,
			size_t MV, size_t NV, typename TV = double
		>
		static constexpr void SVD(const ArcsMat<M,N,T>& A, ArcsMat<MU,NU,TU>& U, ArcsMat<MS,NS,TS>& S, ArcsMat<MV,NV,TV>& V){
			static_assert(MU == NU, "ArcsMat: Size Error");	// 正方行列チェック
			static_assert(MV == NV, "ArcsMat: Size Error");	// 正方行列チェック
			static_assert(M == MU, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(M == MS, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == NS, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == NV, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<TU>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TS>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TV>, "ArcsMat: Type Error");	// 対応可能型チェック
			
			// 初期化
			ArcsMat<N,M,T> Snm;
			ArcsMat<M,M,T> Qm;
			ArcsMat<N,N,T> Qn;
			U = ArcsMat<M,M,T>::eye();
			V = ArcsMat<N,N,T>::eye();
			S = A;
			T E = 0, F = 0;
			
			// ループ打ち切り最大回数に達するまでループ
			for(size_t i = 0; i < LoopMax; ++i){
				// MATLABに準拠するために縦長か横長かで漸化式の順序を変える
				if constexpr( N <= M ){
					// 正方および縦長の場合
					ArcsMat<M,N,T>::QR(S, Qm, S);
					ArcsMat<N,M,T>::QR(~S, Qn, Snm);
					S = ~Snm;
				}else{
					// 横長の場合
					ArcsMat<N,M,T>::QR(~S, Qn, Snm);
					S = ~Snm;
					ArcsMat<M,N,T>::QR(S, Qm, S);
				}
				U = U*Qm;
				V = V*Qn;
				E = ArcsMat<N,M,T>::template norm<ArcsMatrix::NormType::AMT_LINF>( ArcsMat<N,M,T>::gettriup(Snm) );
				F = ArcsMat<std::min(M,N),1,T>::template norm<ArcsMatrix::NormType::AMT_LINF>( ArcsMat<N,M,T>::getdiag(Snm)  );
				if(std::abs(E - F) < ArcsMatrix::EPSILON) break;		// 誤差がイプシロンを下回ったらループ打ち切り
				//printf("%zu: %f\n", i, std::abs(E - F));	// 誤差確認用コード
			}
			
			// 符号修正
			for(size_t k = 1; k <= std::min(M,N); ++k){
				if( std::real( ArcsMatrix::sgn( S(k,k) )) < 0){
					S(k,k) = -S(k,k);
					ArcsMat<M,M,T>::setcolumn(U, -ArcsMat<M,M,T>::getcolumn(U, k), k);
				}
			}

			// Aが正方で且つ複素数の場合には、MATLABに準拠させるための等価変換を実行
			if constexpr((M == N) && ArcsMatrix::IsComplex<T>){
				const auto l  = ArcsMat<1,NV,T>::ones();	// [C++20移行時にconstexprに改修]
				const ArcsMat<1,NV,TV> d = l % ArcsMat<1,NV,TV>::sign( ArcsMat<MV,NV,TV>::getrow(V, 1) );	// V行列の一番上の横ベクトルの符号の逆数を計算
				for(size_t i = 1; i <= NV; ++i){
					ArcsMat<MU,NU,TU>::setcolumn(U, d(1,i)*ArcsMat<MU,NU,TU>::getcolumn(U, i), i);	// U行列の等価変換
					ArcsMat<MV,NV,TV>::setcolumn(V, d(1,i)*ArcsMat<MV,NV,TV>::getcolumn(V, i), i);	// V行列の等価変換
				}
			}
		}

		//! @brief SVD特異値分解(タプル返し版)
		//! @param[in]	A	入力行列
		//! @return	(U, S, V)	(U行列, Σ行列, V行列)のタプル
		static constexpr std::tuple<ArcsMat<M,M,T>, ArcsMat<M,N,T>, ArcsMat<N,N,T>> SVD(const ArcsMat<M,N,T>& A){
			ArcsMat<M,M,T> U;
			ArcsMat<M,N,T> S;
			ArcsMat<N,N,T> V;
			SVD(A, U, S, V);
			return {U, S, V};
		}
		
		//! @brief 行列の階数を返す関数(戻り値返し版のみ)
		//! @param[in]	A	入力行列
		//! @param[in]	eps	ランク許容誤差(デフォルト値 = EPSILON)
		//! @return	結果
		static constexpr size_t rank(const ArcsMat<M,N,T>& A, const T eps = ArcsMatrix::EPSILON){
			static_assert(ArcsMatrix::IsApplicable<T>, "ArcsMat: Type Error");	// 対応可能型チェック
			const auto [U, S, V] = ArcsMat<M,N,T>::SVD(A);	// まず、特異値分解して，
			const auto s = ArcsMat<M,N,T>::getdiag(S);		// 次に、S行列の対角要素を抜き出して、
			return s.GetNumOfNonZero(eps);					// 最後に、非ゼロ要素の数をカウントするとそれが階数と等価
		}
		
		//! @brief 修正コレスキー分解(LDL分解) (引数渡し版)
		//! @tparam	ML, NL, TL, MD, ND, TD	A, L, D行列の高さ, 幅, 要素の型
		//! @param[in]	A	入力行列
		//! @param[out]	L	下三角行列
		//! @param[out]	D	対角行列
		template<size_t ML, size_t NL, typename TL = double, size_t MD, size_t ND, typename TD = double>
		static constexpr void LDL(const ArcsMat<M,N,T>& A, ArcsMat<ML,NL,TL>& L, ArcsMat<MD,ND,TD>& D){
			static_assert(M == N,  "ArcsMat: Size Error");	// 正方行列チェック
			static_assert(ML == M, "ArcsMat: Size Error");	// 正方行列チェック
			static_assert(NL == N, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(MD == M, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ND == N, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TL>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TD>, "ArcsMat: Type Error");	// 対応可能型チェック
			L(1,1) = A(1,1);
			D(1,1) = static_cast<T>(1)/L(1,1);
			
			for(size_t i = 2; i <= N; ++i){
				for(size_t j = 1; j <= i; ++j){
					T lld = A(i,j);
					for(size_t k = 1; k < j; ++k){
						lld -= L(i,k)*L(j,k)*D(k,k);
					}
					L(i,j) = lld;
				}
				D(i,i) = static_cast<T>(1)/L(i,i);
			}
			
			// MATLABに準拠させるための等価変換を実行
			const auto l = ArcsMat<M,1,T>::ones();						// [C++20移行時にconstexprに改修]
			const ArcsMat<M,1,T> d = l % ArcsMat<M,N,T>::getdiag(D);	// 対角要素の逆数を抽出
			L = L*D;
			D = ArcsMat<M,1,T>::diag(d);
		}
		
		//! @brief 修正コレスキー分解(LDL分解) (タプル返し版)
		//! @param[in]	A	入力行列
		//! @return	(L, D)	(L行列, D行列)のタプル
		static constexpr std::tuple<ArcsMat<M,N,T>, ArcsMat<M,N,T>> LDL(const ArcsMat<M,N,T>& A){
			ArcsMat<M,N,T> L, D;
			ArcsMat<M,N,T>::LDL(A, L, D);
			return {L, D};
		}

		//! @brief 修正コレスキー分解(引数渡し版)
		//! @tparam	ML, NL, TL	L行列の高さ, 幅, 要素の型
		//! @param[in]	A	入力行列
		//! @param[out]	R	上三角行列
		template<size_t ML, size_t NL, typename TL = double>
		static constexpr void Cholesky(const ArcsMat<M,N,T>& A, ArcsMat<ML,NL,TL>& R){
			static_assert(M == N,  "ArcsMat: Size Error");	// 正方行列チェック
			static_assert(ML == M, "ArcsMat: Size Error");	// 正方行列チェック
			static_assert(NL == N, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TL>, "ArcsMat: Type Error");	// 対応可能型チェック
			ArcsMat<M,N,T> L, D;
			ArcsMat<M,N,T>::LDL(A, L, D);	// まず、LDL分解してから、
			L = L*ArcsMat<M,N,T>::sqrt(D);	// 次に、対角行列の平方根を取って、下三角行列に掛けたものを出力
			R = ArcsMat<M,N,T>::tp(L);		// MATLABに準拠させるための等価変換を実行
		}

		//! @brief 修正コレスキー分解(戻り値返し版)
		//! @param[in]	A	入力行列
		//! @return	上三角行列
		static constexpr ArcsMat<M,N,T> Cholesky(const ArcsMat<M,N,T>& A){
			ArcsMat<M,N,T> R;
			ArcsMat<M,N,T>::Cholesky(A, R);
			return R;
		}
		
		//! @brief Ax = bの形の線形方程式をxについて解く関数(正方行列A・ベクトルx,b版) (内部用引数渡し版のみ)
		//! @tparam	MB, NB, TB, MX, NX, TX	bとxの高さ, 幅, 要素の型
		//! @param[in]	A	係数行列(正方行列)
		//! @param[in]	b	係数ベクトル
		//! @param[out]	x	解ベクトル
		template<size_t MB, size_t NB, typename TB = double, size_t MX, size_t NX, typename TX = double>
		static constexpr void linsolve_vec(const ArcsMat<M,N,T>& A, const ArcsMat<MB,NB,TB>& b, ArcsMat<MX,NX,TX>& x){
			static_assert(M == N,  "ArcsMat: Size Error");		// 正方行列チェック
			static_assert(MB == M,  "ArcsMat: Size Error");		// 行列のサイズチェック
			static_assert(MX == N,  "ArcsMat: Size Error");		// 行列のサイズチェック
			static_assert(NB == 1, "ArcsMat: Vector Error");	// 縦ベクトルチェック
			static_assert(NX == 1, "ArcsMat: Vector Error");	// 縦ベクトルチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TB>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TX>, "ArcsMat: Type Error");	// 対応可能型チェック
			
			// まず、Ax = b の A をLUP分解すると、
			// 　Ax = b  →  P'LUx = b  →  PP'LUx = Pb  →  LUx = Pb
			// となり、Ux = d、Pb = c と置き換えると、
			// 　Ld = c
			// と書けるので、そこまでの計算をすると下記のコードとなる。
			const auto [L, U, P] = ArcsMat<M,N,T>::LUP(A);
			const auto c = P*b;
			ArcsMat<M,1,T> d;

			// 次に、Ld = c の d について解く。Lは下三角行列で対角が1なので、下記の前進代入法で求まる。
			T buff = 0;			// 行方向加算用のバッファ
			d[1] = c[1];		// 1行目は答えがそのまま現れる
			for(size_t i = 2; i <= M; ++i){
				for(size_t j = 1; j < i; ++j) buff += L(i,j)*d[j];
				d[i] = c[i] - buff;	// i行目の答え
				buff = 0;
			}

			// さらに、Ux = d の x について解く。Uは上三角行列なので、下記の後退代入法で求まる。
			x[M] = d[M]/U(M,M);	// 最後のM行目はUの対角項で割るだけ
			for(ssize_t i = M - 1; 0 < i; --i){
				for(size_t j = static_cast<size_t>(i) + 1; j <= N; ++j){
					buff += U(i,j)*x[j];
				}
				x[i] = (d[i] - buff)/U(i,i);
				buff = 0;
			}
		}
		
		//! @brief AX = Bの形の線形方程式をxについて解く関数(正方行列A・行列X,B版) (内部用引数渡し版のみ)
		//! @tparam	MB, NB, TB, MX, NX, TX	BとXの高さ, 幅, 要素の型
		//! @param[in]	A	係数行列(正方行列)
		//! @param[in]	B	係数行列(正方行列)
		//! @param[out]	X	解行列
		template<size_t MB, size_t NB, typename TB = double, size_t MX, size_t NX, typename TX = double>
		static constexpr void linsolve_mat(const ArcsMat<M,N,T>& A, const ArcsMat<MB,NB,TB>& B, ArcsMat<MX,NX,TX>& X){
			static_assert(M == N,   "ArcsMat: Size Error");		// 正方行列チェック
			static_assert(MB == M,  "ArcsMat: Size Error");		// 行列のサイズチェック
			static_assert(MX == N,  "ArcsMat: Size Error");		// 行列のサイズチェック
			static_assert(NX == NB, "ArcsMat: Size Error");		// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TB>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TX>, "ArcsMat: Type Error");	// 対応可能型チェック

			// 列ごとに線形方程式を解く
			ArcsMat<MX,1,TX> x;
			for(size_t i = 1; i <= NB; ++i){
				ArcsMat<M,N,T>::linsolve_vec(A, ArcsMat<MB,NB,TB>::getcolumn(B, i), x);
				ArcsMat<MX,NX,TX>::setcolumn(X, x, i);
			}
		}
		
		//! @brief Ax = bの形の線形方程式をxについて解く関数(非正方縦長行列A・ベクトルx,b版) (内部用引数渡し版のみ)
		//! @tparam	MB, NB, TB, MX, NX, TX	bとxの高さ, 幅, 要素の型
		//! @param[in]	A	係数行列(非正方縦長行列)
		//! @param[in]	b	係数ベクトル
		//! @param[out]	x	解ベクトル
		template<size_t MB, size_t NB, typename TB = double, size_t MX, size_t NX, typename TX = double>
		static constexpr void linsolve_vec_nsqv(const ArcsMat<M,N,T>& A, const ArcsMat<MB,NB,TB>& b, ArcsMat<MX,NX,TX>& x){
			static_assert(N < M,   "ArcsMat: Size Error");		// 非正方縦長行列チェック
			static_assert(MB == M,  "ArcsMat: Size Error");		// 行列のサイズチェック
			static_assert(MX == N,  "ArcsMat: Size Error");		// 行列のサイズチェック
			static_assert(NB == 1, "ArcsMat: Vector Error");	// 縦ベクトルチェック
			static_assert(NX == 1, "ArcsMat: Vector Error");	// 縦ベクトルチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TB>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TX>, "ArcsMat: Type Error");	// 対応可能型チェック

			// まず、Ax = b の A をQR分解すると、
			// 　Ax = b  →  QRx = b
			// となり、Rx = d と置き換えると、
			// 　Qd = b
			// 　Q'Qd = Q'b  (← Qは直交行列なので Q'Q = I)
			// 　d = Q'b
			// と書けるので、そこまでの計算をすると下記のコードとなる。
			const auto [Q, R] = ArcsMat<M,N>::QR(A);
			const auto d = ~Q*b;

			// 次に、Rx = d の x について解く。Rは縦長上三角行列なので、下記の後退代入法で求まる。
			T buff = 0;			// 行方向加算用のバッファ
			x[N] = d[N]/R(N,N);	// xの最後のN行目は対角項で割るだけ。
			for(ssize_t i = N - 1; 0 < i; --i){
				for(size_t j = static_cast<size_t>(i) + 1; j <= N; ++j){
					buff += R(i,j)*x[j];
				}
				x[i] = (d[i] - buff)/R(i,i);
				buff = 0;
			}
		}
		
		//! @brief AX = Bの形の線形方程式をxについて解く関数(非正方縦長行列A・行列X,B版) (内部用引数渡し版のみ)
		//! @tparam	MB, NB, TB, MX, NX, TX	BとXの高さ, 幅, 要素の型
		//! @param[in]	A	係数行列(非正方縦長行列)
		//! @param[in]	B	係数行列(非正方縦長行列)
		//! @param[out]	X	解行列
		template<size_t MB, size_t NB, typename TB = double, size_t MX, size_t NX, typename TX = double>
		static constexpr void linsolve_mat_nsqv(const ArcsMat<M,N,T>& A, const ArcsMat<MB,NB,TB>& B, ArcsMat<MX,NX,TX>& X){
			static_assert(N < M,   "ArcsMat: Size Error");		// 非正方縦長行列チェック
			static_assert(MB == M,  "ArcsMat: Size Error");		// 行列のサイズチェック
			static_assert(MX == N,  "ArcsMat: Size Error");		// 行列のサイズチェック
			static_assert(NX == NB, "ArcsMat: Size Error");		// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TB>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TX>, "ArcsMat: Type Error");	// 対応可能型チェック

			// 列ごとに線形方程式を解く
			ArcsMat<MX,1,TX> x;
			for(size_t i = 1; i <= NB; ++i){
				ArcsMat<M,N,T>::linsolve_vec_nsqv(A, ArcsMat<MB,NB,TB>::getcolumn(B, i), x);
				ArcsMat<MX,NX,TX>::setcolumn(X, x, i);
			}
		}
		
		//! @brief Ax = bの形の線形方程式をxについて解く関数(非正方横長行列A・ベクトルx,b版) (内部用引数渡し版のみ)
		//!        この関数はMATLABとは異なる解を出力する、ただしもちろん、Ax = b は成立
		//! @tparam	MB, NB, TB, MX, NX, TX	bとxの高さ, 幅, 要素の型
		//! @param[in]	A	係数行列(非正方横長行列)
		//! @param[in]	b	係数ベクトル
		//! @param[out]	x	解ベクトル
		template<size_t MB, size_t NB, typename TB = double, size_t MX, size_t NX, typename TX = double>
		static constexpr void linsolve_vec_nsqh(const ArcsMat<M,N,T>& A, const ArcsMat<MB,NB,TB>& b, ArcsMat<MX,NX,TX>& x){
			static_assert(M < N,   "ArcsMat: Size Error");		// 非正方横長行列チェック
			static_assert(MB == M,  "ArcsMat: Size Error");		// 行列のサイズチェック
			static_assert(MX == N,  "ArcsMat: Size Error");		// 行列のサイズチェック
			static_assert(NB == 1, "ArcsMat: Vector Error");	// 縦ベクトルチェック
			static_assert(NX == 1, "ArcsMat: Vector Error");	// 縦ベクトルチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TB>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TX>, "ArcsMat: Type Error");	// 対応可能型チェック

			// まず、Ax = b の両辺を転置すると、
			// 　(Ax)' = b'  →  x'A' = b'
			// となり、A'をQR分解すると、
			// 　x'QR = b'
			// のように表せて、x'Q = d と置き換えると、
			// 　dR = b'
			// と書けるので、dについて解く。Rは縦長上三角行列なので、下記の前進代入法で求まる。
			const auto [Q, R] = ArcsMat<N,M,T>::QR(~A);
			ArcsMat<1,MX,T> d;
			T buff = 0;				// 行方向加算用のバッファ
			d(1,1) = b[1]/R(1,1);	// 1列目は1列目の対角項で割るだけ
			for(size_t i = 2; i <= M; ++i){
				for(size_t j = 1; j < i; ++j) buff += R(j,i)*d(1,j);
				d(1,i) = (b[i] - buff)/R(i,i);	// i列目の答え
				buff = 0;
			}

			// 次に、
			// 　x'Q = d  →  (x'Q)' = d'  →  Q'x = d'  →  QQ'x = Qd'  → x = Qd'
			// となるので、下記でxが求まる。
			x = Q*~d;
		}
		
		//! @brief AX = Bの形の線形方程式をXについて解く関数(非正方横行列A・行列X,B版) (内部用引数渡し版のみ)
		//!        この関数はMATLABとは異なる解を出力する、ただしもちろん、AX = B は成立
		//! @tparam	MB, NB, TB, MX, NX, TX	BとXの高さ, 幅, 要素の型
		//! @param[in]	A	係数行列(非正方横長行列)
		//! @param[in]	B	係数行列(非正方横長行列)
		//! @param[out]	X	解行列
		template<size_t MB, size_t NB, typename TB = double, size_t MX, size_t NX, typename TX = double>
		static constexpr void linsolve_mat_nsqh(const ArcsMat<M,N,T>& A, const ArcsMat<MB,NB,TB>& B, ArcsMat<MX,NX,TX>& X){
			static_assert(M < N,   "ArcsMat: Size Error");		// 非正方横長行列チェック
			static_assert(MB == M,  "ArcsMat: Size Error");		// 行列のサイズチェック
			static_assert(MX == N,  "ArcsMat: Size Error");		// 行列のサイズチェック
			static_assert(NX == NB, "ArcsMat: Size Error");		// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TB>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TX>, "ArcsMat: Type Error");	// 対応可能型チェック

			// 列ごとに線形方程式を解く
			ArcsMat<MX,1,TX> x;
			for(size_t i = 1; i <= NB; ++i){
				ArcsMat<M,N,T>::linsolve_vec_nsqh(A, ArcsMat<MB,NB,TB>::getcolumn(B, i), x);
				ArcsMat<MX,NX,TX>::setcolumn(X, x, i);
			}
		}

		//! @brief AX = Bの形の線形方程式をXについて解く関数(引数渡し版)
		//! @tparam	MB, NB, TB, MX, NX, TX	BとXの高さ, 幅, 要素の型
		//! @param[in]	A	係数行列(正方行列・非正方行列)
		//! @param[in]	B	係数ベクトル・行列
		//! @param[out]	X	解ベクトル・行列
		template<size_t MB, size_t NB, typename TB = double, size_t MX, size_t NX, typename TX = double>
		static constexpr void linsolve(const ArcsMat<M,N,T>& A, const ArcsMat<MB,NB,TB>& B, ArcsMat<MX,NX,TX>& X){
			// 行列のサイズでアルゴリズムを変える
			if constexpr(M == 1 && N == 1){
				// スカラーの場合
				X[1] = B[1]/A(1,1);	// スカラーのときはそのまま単純に除算
			}else if constexpr((M == N) && (NB == 1)){
				// Aが正方行列で、Bが縦ベクトルの場合
				ArcsMat<M,N,T>::linsolve_vec(A, B, X);		// LU分解を使った線形方程式のベクトル版の解法を使用
			}else if constexpr((M == N) && (NB != 1)){
				// Aが正方行列で、Bも行列の場合
				ArcsMat<M,N,T>::linsolve_mat(A, B, X);		// LU分解を使った線形方程式の行列版の解法を使用
			}else if constexpr((N < M) && (NB == 1)){
				// Aが縦長行列で、Bが縦ベクトルの場合
				ArcsMat<M,N,T>::linsolve_vec_nsqv(A, B ,X);	// QR分解を使った線形方程式のベクトル版の解法を使用
			}else if constexpr((N < M) && (NB != 1)){
				// Aが縦長行列で、Bが行列の場合
				ArcsMat<M,N,T>::linsolve_mat_nsqv(A, B ,X);	// QR分解を使った線形方程式の行列版の解法を使用
			}else if constexpr((M < N) && (NB == 1)){
				// Aが横長行列で、Bが縦ベクトルの場合
				ArcsMat<M,N,T>::linsolve_vec_nsqh(A, B ,X);	// QR分解を使った線形方程式のベクトル版の解法を使用(この場合はMATLABとは異なる解を出力する、ただしもちろん、Ax = b は成立)
			}else if constexpr((M < N) && (NB != 1)){
				// Aが横長行列で、Bが行列の場合
				ArcsMat<M,N,T>::linsolve_mat_nsqh(A, B ,X);	// QR分解を使った線形方程式の行列版の解法を使用(この場合はMATLABとは異なる解を出力する、ただしもちろん、AX = B は成立)
			}else{
				arcs_assert(false);			// ここには来ない
			}
		}
		
		//! @brief AX = Bの形の線形方程式をXについて解く関数(戻り値返し版)
		//! @tparam	MB, NB, TB	B行列の高さ, 幅, 要素の型
		//! @param[in]	A	係数行列(正方行列・非正方行列)
		//! @param[in]	B	係数ベクトル・行列
		//! @return	解ベクトル・行列
		template<size_t MB, size_t NB, typename TB = double>
		static constexpr ArcsMat<N,NB,T> linsolve(const ArcsMat<M,N,T>& A, const ArcsMat<MB,NB,TB>& B){
			ArcsMat<N,NB,T> X;
			ArcsMat<M,N,T>::linsolve(A, B, X);
			return X;
		}

		//! @brief XA = Bの形の線形方程式をXについて解く関数(引数渡し版)
		//! @tparam	MB, NB, TB, MX, NX, TX	BとXの高さ, 幅, 要素の型
		//! @param[in]	A	係数行列(正方行列・非正方行列)
		//! @param[in]	B	係数ベクトル・行列
		//! @param[out]	X	解ベクトル・行列
		template<size_t MB, size_t NB, typename TB = double, size_t MX, size_t NX, typename TX = double>
		static constexpr void linsolveXAB(const ArcsMat<M,N,T>& A, const ArcsMat<MB,NB,TB>& B, ArcsMat<MX,NX,TX>& X){
			// 下記のように転置で変形
			//    XA = B
			// (XA)' = B'
			//  A'X' = B'
			ArcsMat<NX,MX,TX> Xt;					// X'
			ArcsMat<N,M,T>::linsolve(~A, ~B, Xt);	// A'X' = B' を解いて、
			X = ~Xt;	// (X')' として元に戻す
		}
		
		//! @brief XA = Bの形の線形方程式をXについて解く関数(戻り値返し版)
		//! @tparam	MB, NB, TB	B行列の高さ, 幅, 要素の型
		//! @param[in]	A	係数行列(正方行列・非正方行列)
		//! @param[in]	B	係数ベクトル・行列
		//! @return	解ベクトル・行列
		template<size_t MB, size_t NB, typename TB = double>
		static constexpr ArcsMat<MB,M,T> linsolveXAB(const ArcsMat<M,N,T>& A, const ArcsMat<MB,NB,TB>& B){
			ArcsMat<MB,M,T> X;
			ArcsMat<M,N,T>::linsolveXAB(A, B, X);
			return X;
		}

		//! @brief 逆行列を返す関数(引数渡し版)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	A	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void inv(const ArcsMat<M,N,T>& A, ArcsMat<P,Q,R>& Y){
			static_assert(M == N, "ArcsMat: Size Error");	// 正方行列チェック
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			const auto I = ArcsMat<M,N,T>::eye();			// 単位行列の生成 [C++20移行時にconstexprに改修]
			ArcsMat<M,N,T>::linsolve(A, I, Y);				// AX = I となる行列Xを linsolve で見つける
		}

		//! @brief 逆行列を返す関数(戻り値返し版)
		//! @param[in]	A	入力行列
		//! @return	出力行列
		static constexpr ArcsMat<M,N,T> inv(const ArcsMat<M,N,T>& A){
			ArcsMat<M,N,T> Y;
			ArcsMat<M,N,T>::inv(A, Y);
			return Y;
		}
		
		//! @brief Moore-Penroseの擬似逆行列を返す関数(引数渡し版)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	A	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void pinv(const ArcsMat<M,N,T>& A, ArcsMat<P,Q,R>& Y){
			static_assert(M != N, "ArcsMat: Size Error");	// 非正方行列チェック
			static_assert(M == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			if constexpr(N < M){
				// 縦長行列の場合
				Y = ArcsMat<N,N,T>::inv(~A*A)*~A;
			}else if constexpr(M < N){
				// 横長行列の場合
				Y = ~A*ArcsMat<M,M,T>::inv(A*~A);
			}
		}
		
		//! @brief Moore-Penroseの疑似逆行列を返す関数(戻り値返し版)
		//! @param[in]	A	入力行列
		//! @return	出力行列
		static constexpr ArcsMat<N,M,T> pinv(const ArcsMat<M,N,T>& A){
			ArcsMat<N,M,T> Y;
			ArcsMat<M,N,T>::pinv(A, Y);
			return Y;
		}
		
		//! @brief Hessenberg分解(引数渡し版)
		//!        複素数の場合、この関数はMATLABとは異なる解を出力する。
		//! @tparam	MP, NP, TP, MH, NH, TH	出力行列の高さ, 幅, 要素の型
		//! @param[in]	A	入力行列
		//! @param[out]	P	ユニタリ行列
		//! @param[out]	H	ヘッセンベルグ行列
		template<size_t MP, size_t NP, typename TP = double, size_t MH, size_t NH, typename TH = double>
		static constexpr void Hessenberg(const ArcsMat<M,N,T>& A, ArcsMat<MP,NP,TP>& P, ArcsMat<MH,NH,TH>& H){
			static_assert(M == N, "ArcsMat: Size Error");	// 正方行列チェック
			static_assert(MP == M, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(NP == N, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(MH == M, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(NH == N, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TP>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TH>, "ArcsMat: Type Error");	// 対応可能型チェック

			ArcsMat<M,N,T> W;
			ArcsMat<M,1,T> h, u, v;
			P = ArcsMat<M,N,T>::eye();
			H = A;
			for(size_t k = 1; k <= M - 2; ++k){
				h = ArcsMat<M,1,T>::shiftup( ArcsMat<M,N,T>::getcolumn(H, k), k);
				if constexpr(ArcsMatrix::IsComplex<T>){
					u[1] = -std::exp( std::complex( 0.0, std::arg(h[1])) )*std::sqrt( (~h*h)[1] );
				}else{
					if(std::abs(h[1]) < ArcsMatrix::EPSILON){
						u[1] = -std::sqrt( (~h*h)[1] );
					}else{
						u[1] = -ArcsMatrix::sgn(h[1])*std::sqrt( (~h*h)[1] );
					}
				}
				v = h - u;
				ArcsMat<M,1,T>::Householder(v, W, k);
				H = W*H*~W;
				P = P*~W;
			}
		}
		
		//! @brief Hessenberg分解(タプル返し版)
		//!        複素数の場合、この関数はMATLABとは異なる解を出力する。
		//! @param[in]	A	入力行列
		//! @return	(ユニタリ行列P, ヘッセンベルグ行列H) のタプル
		static constexpr std::tuple<ArcsMat<M,N,T>, ArcsMat<M,N,T>> Hessenberg(const ArcsMat<M,N,T>& A){
			ArcsMat<M,N,T> P, H;
			ArcsMat<M,N,T>::Hessenberg(A, P, H);
			return {P, H};
		}

		//! @brief 複素Schur分解(引数渡し版)
		//!        この関数はMATLABとは異なる解を出力する、ただしもちろん、A = USU' は成立
		//! @tparam	MU, NU, TU, MS, NS, TS 	入出力行列の高さ, 幅, 要素の型
		//! @param[in]	A	入力行列
		//! @param[out]	U	ユニタリ行列
		//! @param[out]	S	上三角行列または疑似上三角行列
		template<size_t MU, size_t NU, typename TU = double, size_t MS, size_t NS, typename TS = double>
		static constexpr void Schur(const ArcsMat<M,N,T>& A, ArcsMat<MU,NU,TU>& U, ArcsMat<MS,NS,TS>& S){
			static_assert(M == N, "ArcsMat: Size Error");	// 正方行列チェック
			static_assert(MU == M, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(NU == N, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(MS == M, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(NS == N, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TU>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TS>, "ArcsMat: Type Error");	// 対応可能型チェック

			// ヘッセンベルグ分解を用いた複素シュール分解
			const auto I = ArcsMat<M,N,T>::eye();		// [C++20移行時にconstexprに改修]
			auto [P, H] = ArcsMat<M,N,T>::Hessenberg(A);
			size_t k = M;
			U = P;
			S = H;
			T a = 0, b = 0, c = 0, d= 0, alpha_p = 0, alpha_m = 0, alpha = 0;
			ArcsMat<M,N,T> W, Q, R, V;
			while(1 < k){
				k = (ArcsMat<M,N,T>::getdiag(S, -1)).GetNumOfNonZero() + 1;
				if(1 < k){
					// Wilkinsonシフト量αの計算
					a = S(k,k);
					b = S(k-1,k-1);
					c = S(k-1,k);
					d = S(k,k-1);
					alpha_p = ( a + b + std::sqrt( std::pow(a - b, 2) + 4.0*c*d ) )/2.0;
					alpha_m = ( a + b - std::sqrt( std::pow(a - b, 2) + 4.0*c*d ) )/2.0;
					if( std::abs(alpha_p - a) < std::abs(alpha_m - a) ){
						alpha = alpha_p;
					}else{
						alpha = alpha_m;
					}

					// 原点シフトを用いたQR分解法によりシュール分解
					W.FillAllZero();
					ArcsMat<M,N,T>::copymatrix(S - alpha*I,1,k,1,k, W,1,1);
					ArcsMat<M,N,T>::QR(W, Q, R);
					V = I;
					ArcsMat<M,N,T>::copymatrix(Q,1,k,1,k, V,1,1);
					U = U*V;
					S = ~V*S*V;
					S.ZeroingTriLo();
				}
			}
		}

		//! @brief 複素Schur分解(タプル返し版)
		//!        この関数はMATLABとは異なる解を出力する、ただしもちろん、A = USU' は成立
		//! @param[in]	A	入力行列
		//! @return	(ユニタリ行列U, 上三角行列または疑似上三角行列S)のタプル
		static constexpr std::tuple<ArcsMat<M,N,T>, ArcsMat<M,N,T>> Schur(const ArcsMat<M,N,T>& A){
			ArcsMat<M,N,T> U, S;
			ArcsMat<M,N,T>::Schur(A, U, S);
			return {U, S};
		}
		
		//! @brief 固有値を返す関数(引数渡し版)
		//! @tparam	MV, NV, TV	出力行列の高さ, 幅, 要素の型
		//! @param[in]	A	入力行列
		//! @param[out]	v	固有値の縦ベクトル
		template<size_t MV, size_t NV, typename TV = std::complex<double>>
		static constexpr void eig(const ArcsMat<M,N,T>& A, ArcsMat<MV,NV,TV>& v){
			static_assert(M == N, "ArcsMat: Size Error");	// 正方行列チェック
			static_assert(MV == M, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(NV == 1, "ArcsMat: Vector Error");// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>, "ArcsMat: Type Error");				// 対応可能型チェック
			static_assert(ArcsMatrix::IsComplex<TV>, "ArcsMat: Type Error (Not Complex)");	// 出力は複素数型のみ対応

			const ArcsMat<M,M,TV> U = A;
			const auto [W, S] = ArcsMat<M,N,TV>::Schur(U);	// シュール分解を実行
			v = ArcsMat<M,M,TV>::getdiag(S);				// シュール分解の対角要素が固有値
		}

		//! @brief 固有値を返す関数(戻り値返し版)
		//! @tparam	TV	出力行列の要素の型
		//! @param[in]	A	入力行列
		//! @param[out]	v	固有値の縦ベクトル
		template<typename TV = std::complex<double>>
		static constexpr ArcsMat<M,1,TV> eig(const ArcsMat<M,N,T>& A){
			ArcsMat<M,1,TV> v;
			ArcsMat<M,N,T>::eig(A, v);
			return v;
		}

		//! @brief 固有値を持つ対角行列Dと、各々の固有値に対応する固有ベクトルを持つ行列Vを返す関数(引数渡し版)
		//! @tparam	MV, NV, TV, MD, ND, TD	出力行列の高さ, 幅, 要素の型
		//! @param[in]	A	入力行列
		//! @param[out]	V	各々の固有ベクトルを縦ベクトルとして持つ行列V
		//! @param[out]	D	固有値を対角に持つ対角行列D
		template<size_t MV, size_t NV, typename TV = std::complex<double>, size_t MD, size_t ND, typename TD = std::complex<double>>
		static constexpr void eigvec(const ArcsMat<M,N,T>& A, ArcsMat<MV,NV,TV>& V, ArcsMat<MD,ND,TD>& D){
			static_assert(M == N, "ArcsMat: Size Error");	// 正方行列チェック
			static_assert(MV == NV, "ArcsMat: Size Error");	// 正方行列チェック
			static_assert(MD == ND, "ArcsMat: Size Error");	// 正方行列チェック
			static_assert(MV == M, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(MD == M, "ArcsMat: Vector Error");// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>, "ArcsMat: Type Error");				// 対応可能型チェック
			static_assert(ArcsMatrix::IsComplex<TV>, "ArcsMat: Type Error (Not Complex)");	// 出力は複素数型のみ対応
			static_assert(ArcsMatrix::IsComplex<TD>, "ArcsMat: Type Error (Not Complex)");	// 出力は複素数型のみ対応

			// シュール分解による固有値と固有ベクトルの計算
			const ArcsMat<M,M,TV> Ax = A;			// 予め複素行列として読み込み
			const auto l = ArcsMat<M,M,TV>::eig(A);	// シュール分解で固有値を先に求める
			ArcsMat<M,1,TV>::diag(l, D);			// 固有値を対角に持つ対角行列を生成
			const auto I = ArcsMat<M,M,TV>::eye();	// [C++20移行時にconstexprに改修]
			ArcsMat<M,M,TV> Q, R;
			for(size_t k = 1; k <= M; ++k){
					//(~(Ax - l[k]*I)).Disp("% 8.4f");
				ArcsMat<M,M,TV>::QR( ~(Ax - l[k]*I), Q, R );
					//Q.Disp("% 8.4f");
				ArcsMat<M,M,TV>::setvvector(V, ArcsMat<M,M,TV>::getvvector(Q, 1, M), 1, k);
			}
		}

		//! @brief 固有値を持つ対角行列Dと、各々の固有値に対応する固有ベクトルを持つ行列Vを返す関数(タプル返し版)
		//! @tparam	TV	出力行列の要素の型
		//! @param[in]	A	入力行列
		//! @return	(V, D)	(各々の固有ベクトルを縦ベクトルとして持つ行列V, 固有値を対角に持つ対角行列D) のタプル
		template<typename TV = std::complex<double>>
		static constexpr std::tuple< ArcsMat<M,N,TV>, ArcsMat<M,N,TV> > eigvec(const ArcsMat<M,N,T>& A){
			ArcsMat<M,N,TV> V, D;
			ArcsMat<M,N,T>::eigvec(A, V, D);
			return {V, D};
		}
		
		//! @brief クロネッカー積(引数渡し版)
		//! @tparam	MR, NR, TR, MY, NY, TY	入出力行列の高さ, 幅, 要素の型
		//! @param[in] L 演算子の左側
		//! @param[in] R 演算子の右側
		//! @param[out]	Y 結果
		template<size_t MR, size_t NR, typename TR = double, size_t MY, size_t NY, typename TY = double>
		static constexpr void Kron(const ArcsMat<M,N,T>& L, const ArcsMat<MR,NR,TR>& R, ArcsMat<MY,NY,TY>& Y){
			static_assert(MY == M*MR, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(NY == N*NR, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TR>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TY>, "ArcsMat: Type Error");	// 対応可能型チェック
			ArcsMat<MR,NR,T> A;
			
			// 横方向に小行列で埋めていく
			for(size_t i = 1; i <= N; ++i){
				// 縦方向に小行列で埋めていく
				for(size_t j = 1; j <= M; ++j){
					A = L(j,i)*R;
					ArcsMat<MY,NY,TY>::setsubmatrix(Y, A, (j - 1)*MR + 1, (i - 1)*NR + 1);
				}
			}
		}

		//! @brief クロネッカー積(戻り値返し版)
		//! @tparam	MR, NR, TR, MY, NY, TY	入力行列の高さ, 幅, 要素の型
		//! @param[in] L 演算子の左側
		//! @param[in] R 演算子の右側
		//! @return 結果
		template<size_t MR, size_t NR, typename TR = double>
		static constexpr ArcsMat<M*MR,N*NR,T> Kron(const ArcsMat<M,N,T>& L, const ArcsMat<MR,NR,TR>& R){
			ArcsMat<M*MR,N*NR,T> Y;
			ArcsMat<M,N,T>::Kron(L, R, Y);
			return Y;
		}

		//! @brief クロス積 ベクトル版 (内部用引数渡し版のみ)
		//! @tparam	MR, NR, TR, MY, NY, TY	入出力行列の高さ, 幅, 要素の型
		//! @param[in] l 演算子の左側縦ベクトル
		//! @param[in] r 演算子の右側縦ベクトル
		//! @param[out]	y 結果縦ベクトル
		template<size_t MR, size_t NR, typename TR = double, size_t MY, size_t NY, typename TY = double>
		static constexpr void cross_vec(const ArcsMat<M,N,T>& l, const ArcsMat<MR,NR,TR>& r, ArcsMat<MY,NY,TY>& y){
			static_assert(M  == 3, "ArcsMat: Size Error");	// 行列のサイズチェック(3次元のみ)
			static_assert(MR == 3, "ArcsMat: Size Error");	// 行列のサイズチェック(3次元のみ)
			static_assert(MY == 3, "ArcsMat: Size Error");	// 行列のサイズチェック(3次元のみ)
			static_assert(N  == 1, "ArcsMat: Vector Error");	// 行列のサイズチェック
			static_assert(NR == 1, "ArcsMat: Vector Error");	// 行列のサイズチェック
			static_assert(NY == 1, "ArcsMat: Vector Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TR>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TY>, "ArcsMat: Type Error");	// 対応可能型チェック

			y[1] = l[2]*r[3] - l[3]*r[2];
			y[2] = l[3]*r[1] - l[1]*r[3];
			y[3] = l[1]*r[2] - l[2]*r[1];
		}

		//! @brief クロス積 行列版 (内部用引数渡し版のみ)
		//! @tparam	MR, NR, TR, MY, NY, TY	入出力行列の高さ, 幅, 要素の型
		//! @param[in] L 演算子の左側
		//! @param[in] R 演算子の右側
		//! @param[out]	Y 結果
		template<size_t MR, size_t NR, typename TR = double, size_t MY, size_t NY, typename TY = double>
		static constexpr void cross_mat(const ArcsMat<M,N,T>& L, const ArcsMat<MR,NR,TR>& R, ArcsMat<MY,NY,TY>& Y){
			static_assert(M  == 3, "ArcsMat: Size Error");	// 行列のサイズチェック(3次元のみ)
			static_assert(MR == 3, "ArcsMat: Size Error");	// 行列のサイズチェック(3次元のみ)
			static_assert(MY == 3, "ArcsMat: Size Error");	// 行列のサイズチェック(3次元のみ)
			static_assert(NR == N, "ArcsMat: Vector Error");	// 行列のサイズチェック
			static_assert(NY == N, "ArcsMat: Vector Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TR>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TY>, "ArcsMat: Type Error");	// 対応可能型チェック

			// 列ごとにクロス積を計算
			ArcsMat<MY,1,TY> y;
			for(size_t i = 1; i <= N; ++i){
				ArcsMat<M,1,T>::cross_vec(ArcsMat<M,N,T>::getcolumn(L, i), ArcsMat<MR,NR,TR>::getcolumn(R, i), y);
				ArcsMat<MY,NY,TY>::setcolumn(Y, y, i);
			}
		}

		//! @brief クロス積 (引数渡し版)
		//! @tparam	MR, NR, TR, MY, NY, TY	入出力行列の高さ, 幅, 要素の型
		//! @param[in] L 演算子の左側
		//! @param[in] R 演算子の右側
		//! @param[out]	Y 結果
		template<size_t MR, size_t NR, typename TR = double, size_t MY, size_t NY, typename TY = double>
		static constexpr void cross(const ArcsMat<M,N,T>& L, const ArcsMat<MR,NR,TR>& R, ArcsMat<MY,NY,TY>& Y){
			// ベクトル入力か行列入力かによってアルゴリズムを変える
			if constexpr(N == 1){
				// ベクトルの場合
				ArcsMat<M,N,T>::cross_vec(L, R, Y);	// ベクトル版を使用
			}else{
				// 行列の場合
				ArcsMat<M,N,T>::cross_mat(L, R, Y);	// 行列版を使用
			}
		}

		//! @brief クロス積 (戻り値返し版)
		//! @tparam	MR, NR, TR 入力行列の高さ, 幅, 要素の型
		//! @param[in] L 演算子の左側
		//! @param[in] R 演算子の右側
		//! @return	結果
		template<size_t MR, size_t NR, typename TR = double>
		static constexpr ArcsMat<M,N,T> cross(const ArcsMat<M,N,T>& L, const ArcsMat<MR,NR,TR>& R){
			ArcsMat<M,N,T> Y;
			ArcsMat<M,N,T>::cross(L, R, Y);
			return Y;
		}
		
		//! @brief vec作用素(行列→縦ベクトル) (引数渡し版)
		//! @tparam	MY, NY, TY 出力ベクトルの高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	y	出力ベクトル
		template<size_t MY, size_t NY, typename TY = double>
		static constexpr void vec(const ArcsMat<M,N,T>& U, ArcsMat<MY,NY,TY>& y){
			static_assert(MY == M*N, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(NY == 1, "ArcsMat: Vector Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TY>, "ArcsMat: Type Error");	// 対応可能型チェック
			
			// 横方向に走査
			size_t k = 0;
			for(size_t i = 1; i <= N; ++i){
				// 縦方向に走査
				for(size_t j = 1; j <= M; ++j){
					k++;
					y[k] = U(j,i);
				}
			}
		}

		//! @brief vec作用素(行列→縦ベクトル) (戻り値返し版)
		//! @param[in]	U	入力行列
		//! @param[out]	y	出力ベクトル
		static constexpr ArcsMat<M*N,1,T> vec(const ArcsMat<M,N,T>& U){
			ArcsMat<M*N,1,T> y;
			ArcsMat<M,N,T>::vec(U, y);
			return y;
		}
		
		//! @brief vec作用素の逆(縦ベクトル→行列) (引数渡し版)
		//! @tparam	MY, NY, TY 出力行列の高さ, 幅, 要素の型
		//! @param[in]	u	入力縦ベクトル
		//! @param[out]	Y	再構成後の行列
		template<size_t MY, size_t NY, typename TY = double>
		static constexpr void vecinv(const ArcsMat<M,N,T>& u, ArcsMat<MY,NY,TY>& Y){
			static_assert(M == MY*NY, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == 1, "ArcsMat: Vector Error");		// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TY>, "ArcsMat: Type Error");	// 対応可能型チェック
			
			// 横方向に走査
			size_t k = 0;
			for(size_t i = 1; i <= NY; ++i){
				// 縦方向に走査
				for(size_t j = 1; j <= MY; ++j){
					k++;
					Y(j,i) = u[k];
				}
			}
		}
		
		//! @brief vec作用素の逆(縦ベクトル→行列) (戻り値返し版)
		//! @tparam	MY, NY	出力行列の高さ, 幅
		//! @param[in]	u	入力縦ベクトル
		//! @return	再構成後の行列
		template<size_t MY, size_t NY>
		static constexpr ArcsMat<MY,NY,T> vecinv(const ArcsMat<M,N,T>& u){
			ArcsMat<MY,NY,T> Y;
			ArcsMat<M,N,T>::vecinv(u, Y);
			return Y;
		}
		
		//! @brief 行列指数関数 (引数渡し版)
		//! @tparam	K	パデ近似の次数 (デフォルト値 = 13)
		//! @tparam	MY, NY, TY 出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<size_t K = 13, size_t MY, size_t NY, typename TY = double>
		static constexpr void expm(const ArcsMat<M,N,T>& U, ArcsMat<MY,NY,TY>& Y){
			static_assert(M == N, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(MY == M, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(NY == N, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TY>, "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(0 < K, "ArcsMat: Setting Error");	// パデ近似の次数は1次以上

			// L∞ノルムでスケーリング
			int e = 0;
			std::frexp(ArcsMat<M,N,T>::norm<ArcsMatrix::NormType::AMT_LINF>(U), &e);
			ArcsMat<M,N,T> A;
			if(0 < e){
				A = std::pow(0.5, e + 1)*U;
			}else{
				e = 0;
				A = 0.5*U;
			}

			// 指数行列のパデ近似の計算
			const auto I = ArcsMat<M,N,T>::eye();	// 単位行列の生成
			ArcsMat<M,N,T> L = I, R = I, X = I, cX;
			T c = 1;
			bool signflag = false;					// 係数の符号生成用フラグ
			int coefsign = -1;						// 係数の符号
			for(size_t i = 1; i <= K; ++i){
				c = c*static_cast<T>(K - i + 1) / static_cast<T>(i*(2*K - i + 1));	// 指数行列のパデ近似係数の計算
				X = A*X;									// A^Mの計算
				cX = c*X;									// cM*A^Mの計算
				R += cX;									// R = I + c1*A + c2*A*A + c3*A*A*A + ... + cM*A^M ←パデ近似の分子
				coefsign = 2*static_cast<int>(signflag) - 1;// 係数の符号
				L += static_cast<T>(coefsign)*cX;			// L = I + c1*A + c2*A*A + c3*A*A*A + ... + cM*A^M ←パデ近似の分母
				signflag = !signflag;						// 係数の符号生成用フラグ
			}
			ArcsMat<M,N,T>::linsolve(L, R, Y);		// パデ近似 Y = L/R → 行列パデ近似 Y = inv(L)*R → 線形方程式表示 LY = R

			// スケールを元に戻す
			for(size_t i = 0; i < static_cast<size_t>(e) + 1; ++i) Y *= Y;
		}
		
		//! @brief 行列指数関数 (戻り値返し版)
		//! @tparam	K	パデ近似の次数 (デフォルト値 = 13)
		//! @param[in]	U	入力行列
		//! @return	出力行列
		template<size_t K = 13>
		static constexpr ArcsMat<M,N,T> expm(const ArcsMat<M,N,T>& U){
			ArcsMat<M,N,T> Y;
			ArcsMat<M,N,T>::expm<K>(U, Y);
			return Y;
		}
		
		//! @brief 縦方向の平均を計算する関数 行列入力-横ベクトル出力版 (引数渡し版)
		//! @tparam	MY, NY, TY 出力ベクトルの高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	y	出力ベクトル
		template<size_t MY, size_t NY, typename TY = double>
		static constexpr void meancolumn(const ArcsMat<M,N,T>& U, ArcsMat<MY,NY,TY>& y){
			static_assert(MY == 1, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(NY == N, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TY>, "ArcsMat: Type Error");	// 対応可能型チェック
			y = ArcsMat<M,N,T>::sumcolumn(U)/static_cast<T>(M);	// 縦方向に加算して行列の高さで割って平均を算出
		}

		//! @brief 縦方向の平均を計算する関数 行列入力-横ベクトル出力版 (戻り値返し版)
		//! @tparam	MY, NY, TY 出力ベクトルの高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @return	出力ベクトル
		static constexpr ArcsMat<1,N,T> meancolumn(const ArcsMat<M,N,T>& U){
			ArcsMat<1,N,T> y;
			ArcsMat<M,N,T>::meancolumn(U, y);
			return y;
		}
		
		//! @brief 横方向の平均を計算する関数 行列入力-縦ベクトル出力版 (引数渡し版)
		//! @tparam	MY, NY, TY 出力ベクトルの高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	y	出力ベクトル
		template<size_t MY, size_t NY, typename TY = double>
		static constexpr void meanrow(const ArcsMat<M,N,T>& U, ArcsMat<MY,NY,TY>& y){
			static_assert(MY == M, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(NY == 1, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TY>, "ArcsMat: Type Error");	// 対応可能型チェック
			y = ArcsMat<M,N,T>::sumrow(U)/static_cast<T>(N);	// 横方向に加算して行列の幅で割って平均を算出
		}

		//! @brief 横方向の平均を計算する関数 行列入力-縦ベクトル出力版 (戻り値返し版)
		//! @tparam	MY, NY, TY 出力ベクトルの高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @return	出力ベクトル
		static constexpr ArcsMat<M,1,T> meanrow(const ArcsMat<M,N,T>& U){
			ArcsMat<M,1,T> y;
			ArcsMat<M,N,T>::meanrow(U, y);
			return y;
		}

		//! @brief ベクトルの平均を計算する関数 ベクトル入力-スカラ出力版 (戻り値返し版のみ)
		//! @param[in]	u	入力ベクトル
		//! @return	平均値(スカラー)
		static constexpr T meanvec(const ArcsMat<M,N,T>& u){
			static_assert(M == 1 || N == 1, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック

			// 縦ベクトルか横ベクトルかでアルゴリズムを変える
			ArcsMat<1,1,T> y;
			if constexpr(N == 1){
				// 入力が縦ベクトルの場合
				y = ArcsMat<M,N,T>::meancolumn(u);
			}else if constexpr(M == 1){
				// 入力が横ベクトルの場合
				y = ArcsMat<M,N,T>::meanrow(u);
			}else{
				arcs_assert(false);	// ここには来ない
			}
			
			return y[1];
		}
		
		//! @brief 行列全体の平均を計算する関数 (戻り値返し版のみ)
		//! @param[in]	U	入力行列
		//! @return	平均値(スカラー)
		static constexpr T mean(const ArcsMat<M,N,T>& U){
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			return ArcsMat<1,N,T>::meanvec(ArcsMat<M,N,T>::meancolumn(U));		// 列優先で平均を計算
		}

		//! @brief 縦方向の中央値を計算する関数 行列入力-横ベクトル出力版 (引数渡し版)
		//! @tparam	MY, NY, TY 出力ベクトルの高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	y	出力ベクトル
		template<size_t MY, size_t NY, typename TY = double>
		static constexpr void mediancolumn(const ArcsMat<M,N,T>& U, ArcsMat<MY,NY,TY>& y){
			static_assert(MY == 1, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(NY == N, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TY>, "ArcsMat: Type Error");	// 対応可能型チェック
			ArcsMat<M,N,T> Usrt;
			ArcsMat<M,N,T>::sortcolumn(U, Usrt);		// 縦方向にソート
			if constexpr((M % 2) == 1){
				// 奇数行のとき
				ArcsMat<M,N,T>::getrow(Usrt, y, M/2);	// 中央値の行を抽出 
			}else{
				// 偶数行のとき
				ArcsMat<M,N,T>::getrow(Usrt, y, M/2);	// 中央値の行を抽出
			}
		}

		//! @brief 縦方向の中央値を計算する関数 行列入力-横ベクトル出力版 (戻り値返し版)
		//! @tparam	MY, NY, TY 出力ベクトルの高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @return	出力ベクトル
		static constexpr ArcsMat<1,N,T> mediancolumn(const ArcsMat<M,N,T>& U){
			ArcsMat<1,N,T> y;
			ArcsMat<M,N,T>::mediancolumn(U, y);
			return y;
		}

		//! @brief 横方向の中央値を計算する関数 行列入力-横ベクトル出力版 (引数渡し版)
		//! @tparam	MY, NY, TY 出力ベクトルの高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	y	出力ベクトル
		template<size_t MY, size_t NY, typename TY = double>
		static constexpr void medianrow(const ArcsMat<M,N,T>& U, ArcsMat<MY,NY,TY>& y){
			static_assert(MY == M, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(NY == 1, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TY>, "ArcsMat: Type Error");	// 対応可能型チェック

		}

		//! @brief 行列全体の中央値を計算する関数 (戻り値返し版のみ)
		//! @param[in]	U	入力行列
		//! @return	中央値(スカラー)
		static constexpr T median(const ArcsMat<M,N,T>& U){
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック

		}
		
		//! @brief 横方向の中央値を計算する関数 行列入力-横ベクトル出力版 (戻り値返し版)
		//! @tparam	MY, NY, TY 出力ベクトルの高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @return	出力ベクトル
		static constexpr ArcsMat<M,1,T> medianrow(const ArcsMat<M,N,T>& U){
			ArcsMat<M,1,T> y;
			ArcsMat<M,N,T>::medianrow(U, y);
			return y;
		}
		
		//! @brief 縦方向の分散を計算する関数 行列入力-横ベクトル出力版 (引数渡し版)
		//! @tparam	MY, NY, TY 出力ベクトルの高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	y	出力ベクトル
		template<size_t MY, size_t NY, typename TY = double>
		static constexpr void varcolumn(const ArcsMat<M,N,T>& U, ArcsMat<MY,NY,TY>& y){
			static_assert(MY == 1, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(NY == N, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TY>, "ArcsMat: Type Error");	// 対応可能型チェック
			const ArcsMat<1,N,T> ubar = ArcsMat<M,N,T>::meancolumn(U);	// 列ごとの平均を計算
			ArcsMat<M,N,T> W;
			for(size_t m = 1; m <= M; ++m){
				ArcsMat<M,N,T>::setrow(W, ArcsMat<M,N,T>::getrow(U, m) - ubar, m);	// 列の平均で引く
			}
			y = ArcsMat<M,N,T>::sumcolumn(W & W)/static_cast<T>(M - 1);	// 要素ごとに2乗して列の総和、不偏分散を計算
		}

		//! @brief 縦方向の分散を計算する関数 行列入力-横ベクトル出力版 (戻り値返し版)
		//! @param[in]	U	入力行列
		//! @return	出力ベクトル
		static constexpr ArcsMat<1,N,T> varcolumn(const ArcsMat<M,N,T>& U){
			ArcsMat<1,N,T> y;
			ArcsMat<M,N,T>::varcolumn(U, y);
			return y;
		}

		//! @brief 横方向の分散を計算する関数 行列入力-縦ベクトル出力版 (引数渡し版)
		//! @tparam	MY, NY, TY 出力ベクトルの高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	y	出力ベクトル
		template<size_t MY, size_t NY, typename TY = double>
		static constexpr void varrow(const ArcsMat<M,N,T>& U, ArcsMat<MY,NY,TY>& y){
			static_assert(MY == M, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(NY == 1, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TY>, "ArcsMat: Type Error");	// 対応可能型チェック
			const ArcsMat<M,1,T> ubar = ArcsMat<M,N,T>::meanrow(U);	// 行ごとの平均を計算
			ArcsMat<M,N,T> W;
			for(size_t n = 1; n <= N; ++n){
				ArcsMat<M,N,T>::setcolumn(W, ArcsMat<M,N,T>::getcolumn(U, n) - ubar, n);	// 行の平均で引く
			}
			y = ArcsMat<M,N,T>::sumrow(W & W)/static_cast<T>(N - 1);	// 要素ごとに2乗して行の総和、不偏分散を計算
		}

		//! @brief 横方向の分散を計算する関数 行列入力-縦ベクトル出力版 (戻り値返し版)
		//! @param[in]	U	入力行列
		//! @return	出力ベクトル
		static constexpr ArcsMat<M,1,T> varrow(const ArcsMat<M,N,T>& U){
			ArcsMat<M,1,T> y;
			ArcsMat<M,N,T>::varrow(U, y);
			return y;
		}

		//! @brief ベクトルの分散を計算する関数 ベクトル入力-スカラ出力版 (戻り値返し版のみ)
		//! @param[in]	u	入力ベクトル
		//! @return	分散(スカラー)
		static constexpr T varvec(const ArcsMat<M,N,T>& u){
			static_assert(M == 1 || N == 1, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック

			// 縦ベクトルか横ベクトルかでアルゴリズムを変える
			ArcsMat<1,1,T> y;
			if constexpr(N == 1){
				// 入力が縦ベクトルの場合
				y = ArcsMat<M,N,T>::varcolumn(u);
			}else if constexpr(M == 1){
				// 入力が横ベクトルの場合
				y = ArcsMat<M,N,T>::varrow(u);
			}else{
				arcs_assert(false);	// ここには来ない
			}
			
			return y[1];
		}

		//! @brief 行列全体の分散を計算する関数 (戻り値返し版のみ)
		//! @param[in]	U	入力行列
		//! @return	分散
		static constexpr T var(const ArcsMat<M,N,T>& U){
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			const T Ubar = ArcsMat<M,N,T>::mean(U);	// 要素全体の平均を計算
			const ArcsMat<M,N,T> W = U - Ubar;		// 平均で引く
			return ArcsMat<M,N,T>::sum(W & W)/static_cast<T>(M*N - 1);	// 要素ごとに2乗して総和、不偏分散を計算
		}

		//! @brief 縦方向の標準偏差を計算する関数 行列入力-横ベクトル出力版 (引数渡し版)
		//! @tparam	MY, NY, TY 出力ベクトルの高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	y	出力ベクトル
		template<size_t MY, size_t NY, typename TY = double>
		static constexpr void stdevcolumn(const ArcsMat<M,N,T>& U, ArcsMat<MY,NY,TY>& y){
			static_assert(MY == 1, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(NY == N, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TY>, "ArcsMat: Type Error");	// 対応可能型チェック
			y = ArcsMat<1,N,T>::sqrt( ArcsMat<M,N,T>::varcolumn(U) );	// 分散を計算して平方根を通して出力
		}

		//! @brief 縦方向の標準偏差を計算する関数 行列入力-横ベクトル出力版 (戻り値返し版)
		//! @param[in]	U	入力行列
		//! @return	出力ベクトル
		static constexpr ArcsMat<1,N,T> stdevcolumn(const ArcsMat<M,N,T>& U){
			ArcsMat<1,N,T> y;
			ArcsMat<M,N,T>::stdevcolumn(U, y);
			return y;
		}

		//! @brief 横方向の標準偏差を計算する関数 行列入力-縦ベクトル出力版 (引数渡し版)
		//! @tparam	MY, NY, TY 出力ベクトルの高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	y	出力ベクトル
		template<size_t MY, size_t NY, typename TY = double>
		static constexpr void stdevrow(const ArcsMat<M,N,T>& U, ArcsMat<MY,NY,TY>& y){
			static_assert(MY == M, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(NY == 1, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TY>, "ArcsMat: Type Error");	// 対応可能型チェック
			y = ArcsMat<M,1,T>::sqrt( ArcsMat<M,N,T>::varrow(U) );	// 分散を計算して平方根を通して出力
		}

		//! @brief 横方向の標準偏差を計算する関数 行列入力-縦ベクトル出力版 (戻り値返し版)
		//! @param[in]	U	入力行列
		//! @return	出力ベクトル
		static constexpr ArcsMat<M,1,T> stdevrow(const ArcsMat<M,N,T>& U){
			ArcsMat<M,1,T> y;
			ArcsMat<M,N,T>::stdevrow(U, y);
			return y;
		}

		//! @brief 行列全体の標準偏差を計算する関数 (戻り値返し版のみ)
		//! @param[in]	U	入力行列
		//! @return	標準偏差
		static constexpr T stdev(const ArcsMat<M,N,T>& U){
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			return std::sqrt( ArcsMat<M,N,T>::var(U) );	// 全体の分散を計算して平方根を通して出力
		}
		
		//! @brief 対称ハンケル行列を作成する関数 (引数渡し版)
		//! @tparam	MY, NY, TY 出力行列の高さ, 幅, 要素の型
		//! @param	u	入力ベクトル
		//! @param	Y	出力行列
		template<size_t MY, size_t NY, typename TY = double>
		static constexpr void hankel(const ArcsMat<M,N,T>& u, ArcsMat<MY,NY,TY>& Y){
			static_assert(N == 1, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(MY == M, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(NY == M, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TY>, "ArcsMat: Type Error");	// 対応可能型チェック

			size_t i_start = 1;
			for(size_t m = 1; m <= M; ++m){
				for(size_t i = i_start; i <= M; ++i){
					Y(m, i - i_start + 1) = u[i];
				}
				i_start++;
			}
		}

		//! @brief 対称ハンケル行列を作成する関数 (戻り値返し版)
		//! @param	u	入力ベクトル
		//! @return	出力行列
		static constexpr ArcsMat<M,M,T> hankel(const ArcsMat<M,N,T>& u){
			ArcsMat<M,M,T> Y;
			ArcsMat<M,N,T>::hankel(u, Y);
			return Y;
		}

		//! @brief 根(極)から多項式の係数を計算する関数 (引数渡し版)
		//! MATLABでいうところのpoly関数
		//! 参考文献： Calif. Inst. of Technology, "Compute Polynomial Coeffcients from Roots," Jul. 2015.
		//!  (z - z1)*(z - z2)*...*(z - zn)
		//! 上記の根z1～znから、下記の係数c(1)～c(n+1) を求める。
		//!  c(1)*z^n + ... + c(n)*z + c(n+1)
		//! 使い方の例：
		//!  z[1] = 2 + 1i;
		//!  z[2] = 3 + 2i;
		//!  polycoeff(z, c)
		//!  計算結果： c = {1.0000 + 0.0000i, -5.0000 - 3.0000i, 4.0000 + 7.0000i}
		//! @tparam	MY, NY, TY 出力行列の高さ, 幅, 要素の型
		//!	@param[in]	u	根(極)の値が羅列された入力縦ベクトル
		//! @param[out]	y	多項式(M次方程式)の係数が羅列された出力縦ベクトル
		//! @param[in]	Tol	許容誤差(デフォルト値 1e-10)
		template<size_t MY, size_t NY, typename TY = double, typename R = double>
		static constexpr void polycoeff(const ArcsMat<M,N,T>& u, ArcsMat<MY,NY,TY>& y, const R Tol = 1e-10){
			static_assert(M == MY - 1, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == 1, "ArcsMat: Size Error");		// 行列のサイズチェック
			static_assert(NY == 1, "ArcsMat: Size Error");		// 行列のサイズチェック
			static_assert(ArcsMatrix::IsApplicable<T>,  "ArcsMat: Type Error");	// 対応可能型チェック
			static_assert(ArcsMatrix::IsApplicable<TY>, "ArcsMat: Type Error");	// 対応可能型チェック
			
			// 多項式の係数を求めるラムダ式(関数オブジェクト)
			auto getcoeff = [](const auto& in, auto& out){
				out[1] = 1;
				out[2] = -in[1];
				for(size_t i = 2; i <= M; ++i){
					out[i+1] = -out[i]*in[i];
					for(ssize_t j = i; 2 <= j; --j){
						out[j] = out[j] - out[j-1]*in[i];
					}
				}
				return;
			};

			// 入力引数と出力引数の型によって宣言を変える
			if constexpr (ArcsMatrix::IsComplex<T> && ArcsMatrix::IsComplex<TY>){
				// 入力引数も出力引数も複素数なら、そのまま多項式の係数を求める
				getcoeff(u, y);		// ラムダ式を呼び出す
			}else if constexpr (ArcsMatrix::IsComplex<T> && ArcsMatrix::IsIntFloat<TY>){
				// 入力引数が複素数で、出力引数が実数なら、
				ArcsMat<MY,1,std::complex<TY>> yx;	// 複素数に拡張した係数ベクトルを一旦定義
				getcoeff(u, yx);	// ラムダ式を呼び出す

				// 実数係数に限定する
				yx.ZeroingImag(Tol);						// 浮動小数点演算で出たゼロに近い虚数部を、ゼロに置き換え
				arcs_assert(yx.GetNumOfNonZeroImag() == 0);	// 複素係数になっていないかチェック
				ArcsMat<MY,1,std::complex<TY>>::real(yx, y);// この時点で虚数部がほぼゼロとみなせるので、実数係数として出力
				
			}else if constexpr (ArcsMatrix::IsIntFloat<T> && ArcsMatrix::IsIntFloat<TY>){
				// 入力引数も出力引数も実数なら、そのまま多項式の係数を求める
				getcoeff(u, y);		// ラムダ式を呼び出す
			}else{
				// 入力が実数で出力が複素数はあり得ないので、
				arcs_assert(false);	// ここには来ない
			}
		}
		
		//! @brief 根(極)から多項式の係数を計算する関数 (戻り値返し版)
		//! MATLABでいうところのpoly関数
		//!	@param[in]	u	根(極)の値が羅列された入力縦ベクトル
		//! @param[in]	Tol	許容誤差(デフォルト値 1e-10)
		//! @return	多項式(M次方程式)の係数が羅列された出力縦ベクトル
		static constexpr ArcsMat<M+1,1,T> polycoeff(const ArcsMat<M,N,T>& u, const T Tol = 1e-10){
			ArcsMat<M+1,1,T> y;
			ArcsMat<M,N,T>::polycoeff(u, y, Tol);
			return y;
		}
		
	private:
		// 非公開版基本定数
		static constexpr size_t ITERATION_MAX = 10000;	//!< 反復計算の最大値
		static constexpr char CMPLX_UNIT = 'i';			//!< 虚数単位記号
		
		// 内部処理用
		size_t Nindex;	//!< 横方向カウンタ
		size_t Mindex;	//!< 縦方向カウンタ
		ArcsMatrix::MatStatus Status = ArcsMatrix::MatStatus::AMT_NA;	// 行列の状態

		// 行列データの実体
		// ArcsMatは行列の縦方向にメモリアドレスが連続しているので、縦ベクトル優先。
		std::array<std::array<T, M>, N> Data;//!< データ格納用変数 配列要素の順番は Data[N列(横方向)][M行(縦方向)]
};

namespace ArcsMatrix {
	// グローバル版関数の定義 (ADL問題に当たる可能性があるが便利さ優先)
	
	//! @brief 行列のサイズの表示 (この関数はマクロを介して呼ばれることを想定している)
	//! @tparam	M, N, T	行列の高さ, 幅, 要素の型
	//! @param[in] U		表示する行列
	//! @param[in] varname	変数名
	template <size_t M, size_t N, typename T = double>
	static void dispsize_macro(const ArcsMat<M,N,T>& U, const std::string& varname){
		printf("%s = [ Height: %zu  x  Width: %zu ]\n\n", varname.c_str(), M, N);
	}
	
	//! @brief 行列の要素を表示 (書式指定あり版，この関数はマクロを介して呼ばれることを想定している)
	//! @tparam	M, N, T	行列の高さ, 幅, 要素の型
	//! @param[in] U		表示する行列
	//! @param[in] format	表示形式 (%1.3e とか %5.3f とか printfと同じ)
	//! @param[in] varname	変数名
	template <size_t M, size_t N, typename T = double>
	static void dispf_macro(const ArcsMat<M,N,T>& U, const std::string& format, const std::string& varname){
		printf("%s = \n", varname.c_str());
		U.Disp(format);
	}
	
	//! @brief 行列の要素を表示 (書式指定なし版，この関数はマクロを介して呼ばれることを想定している)
	//! @tparam	M, N, T	行列の高さ, 幅, 要素の型
	//! @param[in] U		表示する行列
	//! @param[in] varname	変数名
	template <size_t M, size_t N, typename T = double>
	static void disp_macro(const ArcsMat<M,N,T>& U, const std::string& varname){
		// データ型によって書式指定子を変える
		// 浮動小数点型の場合
		if constexpr(std::is_floating_point_v<T>){
			dispf_macro(U, "% 4g", varname);
			return;
		}
		
		// int型の場合
		if constexpr(std::is_same_v<T, int>){
			dispf_macro(U, "% d", varname);
			return;
		}
		
		// long型の場合
		if constexpr(std::is_same_v<T, long>){
			dispf_macro(U, "% ld", varname);
			return;
		}
		
		// size_t型の場合
		if constexpr(std::is_same_v<T, size_t>){
			dispf_macro(U, "% zu", varname);
			return;
		}
		
		// 複素数型の場合
		if constexpr(std::is_same_v<T, std::complex<double>> || std::is_same_v<T, std::complex<float>>){
			dispf_macro(U, "% 4g", varname);
			return;
		}
	}

	//! @brief M行N列の単位行列を返す関数
	//! @tparam	M, N, T	行列の高さ, 幅, 要素の型
	//! @return 単位行列
	template <size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> eye(void){
		return ArcsMat<M,N,T>::eye();
	}
	
	//! @brief m行n列の零行列を返す関数
	//! @tparam	M, N, T	行列の高さ, 幅, 要素の型
	//! @return 零行列
	template <size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> zeros(void){
		return ArcsMat<M,N,T>::zeros();
	}
	
	//! @brief m行n列の要素がすべて1の行列を返す関数
	//! @tparam	M, N, T	行列の高さ, 幅, 要素の型
	//! @return 1行列
	template <size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> ones(void){
		return ArcsMat<M,N,T>::ones();
	}
	
	//! @brief 単調増加の縦ベクトルを返す関数
	//! @tparam	M, N, T	行列の高さ, 幅, 要素の型
	//! @return 縦ベクトル
	template <size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> ramp(void){
		return ArcsMat<M,N,T>::ramp();
	}
	
	//! @brief 指定した列から縦ベクトルとして抽出する関数 (引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列の高さ, 幅, 要素の型, 出力ベクトルの高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	y	出力ベクトル
	//! @param[in]	n	抽出したい列
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void getcolumn(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& y, const size_t n){
		ArcsMat<M,N,T>::getcolumn(U, y, n);
	}
	
	//! @brief 指定した列から縦ベクトルとして抽出する関数 (戻り値渡し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[in]	n	抽出したい列
	//! @return	出力ベクトル
	template <size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,1,T> getcolumn(const ArcsMat<M,N,T>& U, const size_t n){
		return ArcsMat<M,N,T>::getcolumn(U, n);
	}
	
	//! @brief 指定した列を縦ベクトルで上書きする関数 (引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入出力行列の高さ, 幅, 要素の型, 入力ベクトルの高さ, 幅, 要素の型
	//! @param[in,out]	UY	入出力行列
	//! @param[in]	u	入力ベクトル
	//! @param[in]	n	上書きしたい列
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void setcolumn(ArcsMat<M,N,T>& UY, const ArcsMat<P,Q,R>& u, const size_t n){
		ArcsMat<M,N,T>::setcolumn(UY, u, n);
	}
	
	//! @brief 指定した列を縦ベクトルで上書きする関数 (戻り値渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列の高さ, 幅, 要素の型, 入力ベクトルの高さ, 幅, 要素の型
	//! @param[in]	u	入力ベクトル
	//! @param[in]	n	上書きしたい列
	//! @param[in]	U	入力行列
	//! @return	出力ベクトル
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr ArcsMat<M,N,T> setcolumn(const ArcsMat<P,Q,R>& u, const size_t n, const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::setcolumn(u, n, U);
	}
	
	//! @brief 指定した列と列を入れ替える関数 (引数渡し版)
	//! @tparam	M, N, T	入出力行列の高さ, 幅, 要素の型
	//! @param[in,out]	UY	入出力行列
	//! @param[in]	n1	指定列1
	//! @param[in]	n2	指定列2
	template <size_t M, size_t N, typename T = double>
	constexpr void swapcolumn(ArcsMat<M,N,T>& UY, const size_t n1, const size_t n2){
		ArcsMat<M,N,T>::swapcolumn(UY, n1, n2);
	}
	
	//! @brief 指定した列と列を入れ替える関数 (戻り値渡し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	n1	指定列1
	//! @param[in]	n2	指定列2
	//! @param[in]	U	入力行列
	//! @return	Y	出力行列
	template <size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> swapcolumn(const size_t n1, const size_t n2, const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::swapcolumn(n1, n2, U);
	}

	//! @brief n列目のm1行目からm2行目までを数値aで埋める関数 (m1 <= m2 であること) (引数渡し版)
	//! @tparam	M, N, T, R	入出力行列の高さ, 幅, 要素の型, 埋める値の型
	//! @param[in,out]	UY	入出力行列
	//! @param[in]	a	埋める値
	//! @param[in]	n	指定列
	//! @param[in]	m1	開始行
	//! @param[in]	m2	終了行
	template<size_t M, size_t N, typename T = double, typename R = double>
	constexpr void fillcolumn(ArcsMat<M,N,T>& UY, const R a, const size_t n, const size_t m1, const size_t m2){
		ArcsMat<M,N,T>::fillcolumn(UY, a, n, m1, m2);
	}
	
	//! @brief n列目のm1行目からm2行目までを数値aで埋める関数 (m1 <= m2 であること) (引数渡し版)
	//! @tparam	M, N, T, R	入出力行列の高さ, 幅, 要素の型, 埋める値の型
	//! @param[in]	a	埋める値
	//! @param[in]	n	指定列
	//! @param[in]	m1	開始行
	//! @param[in]	m2	終了行
	//! @param[in]	U	入力行列
	//! @return	出力行列
	template<size_t M, size_t N, typename T = double, typename R = double>
	constexpr ArcsMat<M,N,T> fillcolumn(const R a, const size_t n, const size_t m1, const size_t m2, const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::fillcolumn(a, n, m1, m2, U);
	}
	
	//! @brief 並び替え指定横ベクトルuが昇順になるように，行列Uの列を並び替える関数 (戻り値渡し版のみ)
	//! @tparam	M, N, T, P, Q, R	入力行列の高さ, 幅, 要素の型, 並び替え指定横ベクトルの高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[in]	u	並び替え指定横ベクトル
	//! @return	出力行列
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = size_t>
	constexpr ArcsMat<M,N,T> ordercolumn(const ArcsMat<M,N,T>& U, const ArcsMat<P,Q,R>& u){
		return ArcsMat<M,N,T>::ordercolumn(U, u);
	}
	
	//! @brief 並び替え指定横ベクトルuが昇順になるように，行列Uの列と指定横ベクトルの両方を並び替える関数 (引数渡し版のみ)
	//! @tparam	M, N, T, P, Q, R	入出力行列の高さ, 幅, 要素の型, 並び替え指定横ベクトルの高さ, 幅, 要素の型
	//! @param[in,out]	UY	入出力行列
	//! @param[in]	uy	並び替え指定横ベクトル
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = size_t>
	constexpr void ordercolumn_and_vec(ArcsMat<M,N,T>& UY, ArcsMat<P,Q,R>& uy){
		ArcsMat<M,N,T>::ordercolumn_and_vec(UY, uy);
	}
	
	//! @brief 各列の総和を計算して横ベクトルを出力する関数 (引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列の高さ, 幅, 要素の型, 出力ベクトルの高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	y	出力ベクトル
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void sumcolumn(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& y){
		ArcsMat<M,N,T>::sumcolumn(U, y);
	}
	
	//! @brief 各列の総和を計算して横ベクトルを出力する関数 (戻り値渡し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	出力ベクトル
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<1,N,T> sumcolumn(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::sumcolumn(U);
	}

	//! @brief 指定した行から横ベクトルとして抽出する関数 (引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列の高さ, 幅, 要素の型, 出力ベクトルの高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	y	出力ベクトル
	//! @param[in]	m	抽出したい行
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void getrow(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& y, const size_t m){
		ArcsMat<M,N,T>::getrow(U, y, m);
	}
	
	//! @brief 指定した行から横ベクトルとして抽出する関数 (戻り値渡し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[in]	m	抽出したい行
	//! @return	出力ベクトル
	template <size_t M, size_t N, typename T = double>
	constexpr ArcsMat<1,N,T> getrow(const ArcsMat<M,N,T>& U, const size_t m){
		return ArcsMat<M,N,T>::getrow(U, m);
	}
	
	//! @brief 指定した行を横ベクトルで上書きする関数 (引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入出力行列の高さ, 幅, 要素の型, 入力ベクトルの高さ, 幅, 要素の型
	//! @param[in,out]	UY	入出力行列
	//! @param[in]	u	入力ベクトル
	//! @param[in]	m	上書きしたい行
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void setrow(ArcsMat<M,N,T>& UY, const ArcsMat<P,Q,R>& u, const size_t m){
		ArcsMat<M,N,T>::setrow(UY, u, m);
	}
	
	//! @brief 指定した行を横ベクトルで上書きする関数 (戻り値渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列の高さ, 幅, 要素の型, 入力ベクトルの高さ, 幅, 要素の型
	//! @param[in]	u	入力ベクトル
	//! @param[in]	m	上書きしたい行
	//! @param[in]	U	入力行列
	//! @return	出力ベクトル
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr ArcsMat<M,N,T> setrow(const ArcsMat<P,Q,R>& u, const size_t m, const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::setrow(u, m, U);
	}
	
	//! @brief 指定した行と行を入れ替える関数 (引数渡し版)
	//! @tparam	M, N, T	入出力行列の高さ, 幅, 要素の型
	//! @param[in,out]	UY	入出力行列
	//! @param[in]	m1	指定行1
	//! @param[in]	m2	指定行2
	template <size_t M, size_t N, typename T = double>
	constexpr void swaprow(ArcsMat<M,N,T>& UY, const size_t m1, const size_t m2){
		ArcsMat<M,N,T>::swaprow(UY, m1, m2);
	}
	
	//! @brief 指定した行と行を入れ替える関数 (戻り値渡し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	m1	指定行1
	//! @param[in]	m2	指定行2
	//! @param[in]	U	入力行列
	//! @return	Y	出力行列
	template <size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> swaprow(const size_t m1, const size_t m2, const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::swaprow(m1, m2, U);
	}

	//! @brief m行目のn1列目からn2列目までを数値aで埋める関数 (n1 <= n2 であること) (引数渡し版)
	//! @tparam	M, N, T, R	入出力行列の高さ, 幅, 要素の型, 埋める値の型
	//! @param[in,out]	UY	入出力行列
	//! @param[in]	a	埋める値
	//! @param[in]	m	指定行
	//! @param[in]	n1	開始列
	//! @param[in]	n2	終了列
	template<size_t M, size_t N, typename T = double, typename R = double>
	constexpr void fillrow(ArcsMat<M,N,T>& UY, const R a, const size_t m, const size_t n1, const size_t n2){
		ArcsMat<M,N,T>::fillrow(UY, a, m, n1, n2);
	}
	
	//! @brief m行目のn1列目からn2列目までを数値aで埋める関数 (n1 <= n2 であること) (戻り値渡し版)
	//! @tparam	M, N, T, R	入出力行列の高さ, 幅, 要素の型, 埋める値の型
	//! @param[in]	a	埋める値
	//! @param[in]	m	指定行
	//! @param[in]	n1	開始列
	//! @param[in]	n2	終了列
	//! @param[in]	U	入力行列
	//! @return	出力行列
	template<size_t M, size_t N, typename T = double, typename R = double>
	constexpr ArcsMat<M,N,T> fillrow(const R a, const size_t m, const size_t n1, const size_t n2, const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::fillrow(a, m, n1, n2, U);
	}
	
	//! @brief 並び替え指定縦ベクトルuが昇順になるように，行列Uの行を並び替える関数 (戻り値渡し版のみ)
	//! @tparam	M, N, T, P, Q, R	入力行列の高さ, 幅, 要素の型, 並び替え指定縦ベクトルの高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[in]	u	並び替え指定縦ベクトル
	//! @return	出力行列
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = size_t>
	constexpr ArcsMat<M,N,T> orderrow(const ArcsMat<M,N,T>& U, const ArcsMat<P,Q,R>& u){
		return ArcsMat<M,N,T>::orderrow(U, u);
	}
	
	//! @brief 並び替え指定縦ベクトルuが昇順になるように，行列Uの行と指定縦ベクトルの両方を並び替える関数 (引数渡し版のみ)
	//! @tparam	M, N, T, P, Q, R	入出力行列の高さ, 幅, 要素の型, 並び替え指定縦ベクトルの高さ, 幅, 要素の型
	//! @param[in,out]	UY	入出力行列
	//! @param[in]	uy	並び替え指定縦ベクトル
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = size_t>
	constexpr void orderrow_and_vec(ArcsMat<M,N,T>& UY, ArcsMat<P,Q,R>& uy){
		ArcsMat<M,N,T>::orderrow_and_vec(UY, uy);
	}
	
	//! @brief 各行の総和を計算して縦ベクトルを出力する関数 (引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列の高さ, 幅, 要素の型, 出力ベクトルの高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	y	出力ベクトル
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void sumrow(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& y){
		ArcsMat<M,N,T>::sumrow(U, y);
	}
	
	//! @brief 各行の総和を計算して縦ベクトルを出力する関数 (戻り値渡し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	出力ベクトル
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,1,T> sumrow(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::sumrow(U);
	}
	
	//! @brief 指定位置から縦ベクトルを抽出する関数 (引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列の高さ, 幅, 要素の型, 出力ベクトルの高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[in]	y	出力ベクトル
	//! @param[in]	m	先頭位置 m行目
	//! @param[in]	n	先頭位置 n列目
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void getvvector(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& y, const size_t m, const size_t n){
		ArcsMat<M,N,T>::getvvector(U, y, m, n);
	}
	
	//! @brief 指定位置から縦ベクトルを抽出する関数 (戻り値渡し版)
	//! @tparam	P, M, N, T	出力ベクトルの幅, 入力行列の高さ, 幅, 要素の型, 出力ベクトルの幅
	//! @param[in]	U	入力行列
	//! @param[in]	m	先頭位置 m行目
	//! @param[in]	n	先頭位置 n列目
	//! @return	縦ベクトル
	template<size_t P, size_t M, size_t N, typename T = double>
	constexpr ArcsMat<P,1,T> getvvector(const ArcsMat<M,N,T>& U, const size_t m, const size_t n){
		return ArcsMat<M,N,T>::template getvvector<P>(U, m, n);
	}
	
	//! @brief 指定位置に縦ベクトルで上書きする関数 (引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列の高さ, 幅, 要素の型, 入力ベクトルの高さ, 幅, 要素の型
	//! @param[in,out]	UY	入出力行列
	//! @param[in]	u	入力ベクトル
	//! @param[in]	m	先頭位置 m行目
	//! @param[in]	n	先頭位置 n列目
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void setvvector(ArcsMat<M,N,T>& UY, const ArcsMat<P,Q,R>& u, const size_t m, const size_t n){
		ArcsMat<M,N,T>::setvvector(UY, u, m, n);
	}
	
	//! @brief 指定位置に縦ベクトルで上書きする関数 (戻り値渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列の高さ, 幅, 要素の型, 入力ベクトルの高さ, 幅, 要素の型
	//! @param[in]	u	入力ベクトル
	//! @param[in]	m	先頭位置 m行目
	//! @param[in]	n	先頭位置 n列目
	//! @param[in]	U	入力行列
	//! @return	出力行列
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr ArcsMat<M,N,T> setvvector(const ArcsMat<P,Q,R>& u, const size_t m, const size_t n, const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::setvvector(u, m, n, U);
	}
	
	//! @brief 指定位置から横ベクトルを抽出する関数 (引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列の高さ, 幅, 要素の型, 出力ベクトルの高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[in]	y	出力ベクトル
	//! @param[in]	m	先頭位置 m行目
	//! @param[in]	n	先頭位置 n列目
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void gethvector(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& y, const size_t m, const size_t n){
		ArcsMat<M,N,T>::gethvector(U, y, m, n);
	}
	
	//! @brief 指定位置から横ベクトルを抽出する関数 (戻り値渡し版)
	//! @tparam	Q, M, N, T	出力ベクトルの幅, 入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[in]	m	先頭位置 m行目
	//! @param[in]	n	先頭位置 n列目
	//! @return	横ベクトル
	template<size_t Q, size_t M, size_t N, typename T = double>
	constexpr ArcsMat<1,Q,T> gethvector(const ArcsMat<M,N,T>& U, const size_t m, const size_t n){
		return ArcsMat<M,N,T>::template gethvector<Q>(U, m, n);
	}
	
	//! @brief 指定位置に横ベクトルで上書きする関数 (引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列の高さ, 幅, 要素の型, 入力ベクトルの高さ, 幅, 要素の型
	//! @param[in,out]	UY	入出力行列
	//! @param[in]	u	入力ベクトル
	//! @param[in]	m	先頭位置 m行目
	//! @param[in]	n	先頭位置 n列目
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void sethvector(ArcsMat<M,N,T>& UY, const ArcsMat<P,Q,R>& u, const size_t m, const size_t n){
		ArcsMat<M,N,T>::sethvector(UY, u, m, n);
	}
	
	//! @brief 指定位置に横ベクトルで上書きする関数 (戻り値渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列の高さ, 幅, 要素の型, 入力ベクトルの高さ, 幅, 要素の型
	//! @param[in]	u	入力ベクトル
	//! @param[in]	m	先頭位置 m行目
	//! @param[in]	n	先頭位置 n列目
	//! @param[in]	U	入力行列
	//! @return	出力行列
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr ArcsMat<M,N,T> sethvector(const ArcsMat<P,Q,R>& u, const size_t m, const size_t n, const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::sethvector(u, m, n, U);
	}
		
	//! @brief 行列から指定位置の小行列を抽出する関数 (引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列の高さ, 幅, 要素の型, 小行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	抽出した小行列
	//! @param[in]	m	抽出する縦位置（小行列の上）
	//! @param[in]	n	抽出する横位置（小行列の左）
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void getsubmatrix(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y, const size_t m, const size_t n){
		ArcsMat<M,N,T>::getsubmatrix(U, Y, m, n);
	}
	
	//! @brief 行列から指定位置の小行列を抽出する関数 (戻り値渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列の高さ, 幅, 要素の型, 小行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[in]	m	抽出する縦位置（小行列の上）
	//! @param[in]	n	抽出する横位置（小行列の左）
	//! @return	抽出した小行列
	template<size_t P, size_t Q, size_t M, size_t N, typename T = double>
	constexpr ArcsMat<P,Q,T> getsubmatrix(const ArcsMat<M,N,T>& U, const size_t m, const size_t n){
		return ArcsMat<M,N,T>::template getsubmatrix<P,Q>(U, m, n);
	}
	
	//! @brief 小行列を行列の指定位置に上書きする関数 (引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入出力行列の高さ, 幅, 要素の型, 入力行列の高さ, 幅, 要素の型
	//! @param[in,out]	UY	入出力行列
	//! @param[in]	U	書き込む小行列
	//! @param[in]	m	上書きしたい縦位置（小行列の上）
	//! @param[in]	n	上書きしたい横位置（小行列の左）
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void setsubmatrix(ArcsMat<M,N,T>& UY, const ArcsMat<P,Q,R>& U, const size_t m, const size_t n){
		ArcsMat<M,N,T>::setsubmatrix(UY, U, m, n);
	}
	
	//! @brief 小行列を行列の指定位置に上書きする関数 (戻り値渡し版)
	//! @tparam	M, N, T, P, Q, R	入出力行列の高さ, 幅, 要素の型, 入力行列の高さ, 幅, 要素の型
	//! @param[in]	Us	書き込む小行列
	//! @param[in]	m	上書きしたい縦位置（小行列の上）
	//! @param[in]	n	上書きしたい横位置（小行列の左）
	//! @param[in]	U	入力行列
	//! @return 出力行列	
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr ArcsMat<M,N,T> setsubmatrix(const ArcsMat<P,Q,R>& Us, const size_t m, const size_t n, const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::setsubmatrix(Us, m, n, U);
	}

	//! @brief 行列Uから別の行列Yへ位置とサイズを指定してコピーする関数(引数渡し版のみ)
	//! 等価なMATLABコード: Y(my:my+(m2-m1), ny:ny+(n2-n1)) = U(m1:m2, n1:n2)
	//! @tparam	M, N, T, P, Q, R	入出力行列の高さ, 幅, 要素の型, 入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	コピー元の行列
	//! @param[in]	m1	コピー元の縦方向の抽出開始位置
	//! @param[in]	m2	コピー元の縦方向の抽出終了位置
	//! @param[in]	n1	コピー元の横方向の抽出開始位置
	//! @param[in]	n2	コピー元の横方向の抽出終了位置
	//! @param[in,out]	Y	コピー先の行列
	//! @param[in]	my	コピー先の縦方向の書き込み開始位置
	//! @param[in]	ny	コピー先の横方向の書き込み開始位置
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void copymatrix(
		const ArcsMat<M,N,T>& U, const size_t m1, const size_t m2, const size_t n1, const size_t n2,
		ArcsMat<P,Q,R>& Y, const size_t my, const size_t ny
	){
		ArcsMat<M,N,T>::copymatrix(U, m1, m2, n1, n2, Y, my, ny);
	}

	//! @brief 行列の各要素を上に1行分シフトする関数(下段の行はゼロになる)(引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列の高さ, 幅, 要素の型, 出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	出力行列
	//! @param[in]	m	シフトする行数(デフォルト値 = 1)
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void shiftup(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y, const size_t m = 1){
		ArcsMat<M,N,T>::shiftup(U, Y, m);
	}
	
	//! @brief 行列の各要素を上に1行分シフトする関数(下段の行はゼロになる)(戻り値渡し版)
	//! @tparam	M, N, T	入出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[in]	m	シフトする行数(デフォルト値 = 1)
	//! @return	Y	出力行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> shiftup(const ArcsMat<M,N,T>& U, const size_t m = 1){
		return ArcsMat<M,N,T>::shiftup(U, m);
	}
	
	//! @brief 行列の各要素を下にm行分シフトする関数(上段の行はゼロになる)(引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列の高さ, 幅, 要素の型, 出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	出力行列
	//! @param[in]	m	シフトする行数(デフォルト値 = 1)
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void shiftdown(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y, const size_t m = 1){
		ArcsMat<M,N,T>::shiftdown(U, Y, m);
	}
	
	//! @brief 行列の各要素を下にm行分シフトする関数(上段の行はゼロになる)(戻り値渡し版)
	//! @tparam	M, N, T	入出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[in]	m	シフトする行数(デフォルト値 = 1)
	//! @return	出力行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> shiftdown(const ArcsMat<M,N,T>& U, const size_t m = 1){
		return ArcsMat<M,N,T>::shiftdown(U, m);
	}
	
	//! @brief 行列の各要素を左にn列分シフトする関数(右段の列はゼロになる)(引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列の高さ, 幅, 要素の型, 出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	出力行列
	//! @param[in]	n	シフトする列数(デフォルト値 = 1)
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void shiftleft(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y, const size_t n = 1){
		ArcsMat<M,N,T>::shiftleft(U, Y, n);
	}
	
	//! @brief 行列の各要素を左にn列分シフトする関数(右段の列はゼロになる)(戻り値渡し版)
	//! @tparam	M, N, T	入出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[in]	n	シフトする列数(デフォルト値 = 1)
	//! @return	Y	出力行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> shiftleft(const ArcsMat<M,N,T>& U, const size_t n = 1){
		return ArcsMat<M,N,T>::shiftleft(U, n);
	}
	
	//! @brief 行列の各要素を右にn列分シフトする関数(左段の列はゼロになる)(引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列の高さ, 幅, 要素の型, 出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	出力行列
	//! @param[in]	n	シフトする列数(デフォルト値 = 1)
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void shiftright(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y, const size_t n = 1){
		ArcsMat<M,N,T>::shiftright(U, Y, n);
	}

	//! @brief 行列の各要素を右にn列分シフトする関数(左段の列はゼロになる)(戻り値渡し版)
	//! @tparam	M, N, T	入出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[in]	n	シフトする列数(デフォルト値 = 1)
	//! @return	出力行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> shiftright(const ArcsMat<M,N,T>& U, const size_t n = 1){
		return ArcsMat<M,N,T>::shiftright(U, n);
	}

	//! @brief 行列を縦方向にソートする関数 (引数渡し版)
	//! @tparam	ST	ソート方法： AMT_ASCENT = 昇順(デフォルト), AMT_DESCENT = 降順
	//! @tparam	M, N, T, P, Q, R	入力行列の高さ, 幅, 要素の型, 出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	出力行列
	template<ArcsMatrix::SortType ST = ArcsMatrix::SortType::AMT_ASCENT, size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void sortcolumn(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
		ArcsMat<M,N,T>::template sortcolumn<ST>(U, Y);
	}

	//! @brief 行列を縦方向にソートする関数 (戻り値返し版)
	//! @tparam	ST	ソート方法： AMT_ASCENT = 昇順(デフォルト), AMT_DESCENT = 降順
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	出力行列
	template<ArcsMatrix::SortType ST = ArcsMatrix::SortType::AMT_ASCENT, size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> sortcolumn(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::sortcolumn(U);
	}

	//! @brief 行列を横方向にソートする関数 (引数渡し版)
	//! @tparam	ST	ソート方法： AMT_ASCENT = 昇順(デフォルト), AMT_DESCENT = 降順
	//! @tparam	M, N, T, P, Q, R	入力行列の高さ, 幅, 要素の型, 出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	出力行列
	template<ArcsMatrix::SortType ST = ArcsMatrix::SortType::AMT_ASCENT, size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void sortrow(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
		ArcsMat<M,N,T>::template sortrow<ST>(U, Y);
	}

	//! @brief 行列を横方向にソートする関数 (戻り値返し版)
	//! @tparam	ST	ソート方法： AMT_ASCENT = 昇順(デフォルト), AMT_DESCENT = 降順
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	出力行列
	template<ArcsMatrix::SortType ST = ArcsMatrix::SortType::AMT_ASCENT, size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> sortrow(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::sortrow(U);
	}

	//! @brief 行列を縦に連結する関数(引数渡し版)
	//! @tparam	M, N, T, P, Q, R, D, E, F	入力行列1と2と出力行列の高さ, 幅, 要素の型
	//! @param[in]	U1	入力行列1
	//! @param[in]	U2	入力行列2
	//! @param[out]	Y	出力行列
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double, size_t D, size_t E, typename F = double>
	constexpr void concatv(const ArcsMat<M,N,T>& U1, const ArcsMat<P,Q,R>& U2, ArcsMat<D,E,F>& Y){
		ArcsMat<M,N,T>::concatv(U1, U2, Y);
	}
	
	//! @brief 行列を縦に連結する関数(戻り値渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列1と2の高さ, 幅, 要素の型
	//! @param[in]	U1	入力行列1
	//! @param[in]	U2	入力行列2
	//! @param[out]	Y	出力行列
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	static constexpr ArcsMat<M+P,N,T> concatv(const ArcsMat<M,N,T>& U1, const ArcsMat<P,Q,R>& U2){
		return ArcsMat<M,N,T>::concatv(U1, U2);
	}
	
	//! @brief 行列を横に連結する関数(引数渡し版)
	//! @tparam	M, N, T, P, Q, R, D, E, F	入力行列1と2と出力行列の高さ, 幅, 要素の型
	//! @param[in]	U1	入力行列1
	//! @param[in]	U2	入力行列2
	//! @param[out]	Y	出力行列
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double, size_t D, size_t E, typename F = double>
	constexpr void concath(const ArcsMat<M,N,T>& U1, const ArcsMat<P,Q,R>& U2, ArcsMat<D,E,F>& Y){
		ArcsMat<M,N,T>::concath(U1, U2, Y);
	}
	
	//! @brief 行列を横に連結する関数(戻り値渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列1と2の高さ, 幅, 要素の型
	//! @param[in]	U1	入力行列1
	//! @param[in]	U2	入力行列2
	//! @param[out]	Y	出力行列
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr ArcsMat<M,N+Q,T> concath(const ArcsMat<M,N,T>& U1, const ArcsMat<P,Q,R>& U2){
		return ArcsMat<M,N,T>::concath(U1, U2);
	}
	
	//! @brief 4つの行列を1つに連結する関数(戻り値渡し版)
	//! @tparam	M, N, T, P, Q, R, D, E, F, G, H, L, V, W, X	入力行列11-22と出力行列の高さ, 幅, 要素の型
	//! @param[in]	U11	入力行列11(左上)
	//! @param[in]	U12	入力行列12(右上)
	//! @param[in]	U21	入力行列11(左下)
	//! @param[in]	U22	入力行列12(右下)
	//! @param[out]	Y	出力行列
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double, size_t D, size_t E, typename F = double, size_t G, size_t H, typename L = double, size_t V, size_t W, typename X = double>
	constexpr void concat4(const ArcsMat<M,N,T>& U11, const ArcsMat<P,Q,R>& U12, const ArcsMat<D,E,F>& U21, const ArcsMat<G,H,L>& U22, ArcsMat<V,W,X>& Y){
		ArcsMat<M,N,T>::concat4(U11, U12, U21, U22, Y);
	}
	
	//! @brief 4つの行列を1つに連結する関数(引数渡し版)
	//! @tparam	M, N, T, P, Q, R, D, E, F, G, H, L	入力行列11-22と出力行列の高さ, 幅, 要素の型
	//! @param[in]	U11	入力行列11(左上)
	//! @param[in]	U12	入力行列12(右上)
	//! @param[in]	U21	入力行列11(左下)
	//! @param[in]	U22	入力行列12(右下)
	//! @return	出力行列
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double, size_t D, size_t E, typename F = double, size_t G, size_t H, typename L = double>
	constexpr ArcsMat<M+D,N+Q,T> concat4(const ArcsMat<M,N,T>& U11, const ArcsMat<P,Q,R>& U12, const ArcsMat<D,E,F>& U21, const ArcsMat<G,H,L>& U22){
		return ArcsMat<M,N,T>::concat4(U11, U12, U21, U22);
	}
	
	//! @brief 縦ベクトルの各要素を対角要素に持つ正方行列を生成する関数(引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力ベクトルと出力行列の高さ, 幅, 要素の型
	//! @param[in]	u	入力ベクトル
	//! @param[in]	Y	出力行列
	//! @return	出力行列
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void diag(const ArcsMat<M,N,T>& u, ArcsMat<P,Q,R>& Y){
		ArcsMat<M,N,T>::diag(u, Y);
	}
	
	//! @brief 縦ベクトルの各要素を対角要素に持つ正方行列を生成する関数(戻り値渡し版)
	//! @tparam	M, N, T 正方行列のサイズ, 要素の型
	//! @param[in]	u	入力ベクトル
	//! @return	出力行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,M,T> diag(const ArcsMat<M,N,T>& u){
		return ArcsMat<M,N,T>::diag(u);
	}
	
	//! @brief 行列の対角要素を縦ベクトルとして取得する関数(引数渡し版)
	//! @tparam	M, N, T, P, Q, R, L	入力行列と出力ベクトルの高さ, 幅, 要素の型, 対角要素の数
	//! @param[in]	U	入力行列
	//! @param[in]	k	k番目の対角, k=0で主対角、K<0で主対角より下、0<kで主対角より上 (デフォルト値 = 0)
	//! @param[out]	y	出力ベクトル
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double, size_t L = std::min(M,N)>
	constexpr void getdiag(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& y, const ssize_t k = 0){
		ArcsMat<M,N,T>::getdiag(U, y, k);
	}
	
	//! @brief 行列の対角要素を縦ベクトルとして取得する関数(戻り値渡し版)
	//! @tparam	M, N, T, L	入力行列の高さ, 幅, 要素の型, 対角要素の数
	//! @param[in]	U	入力行列
	//! @param[in]	k	k番目の対角, k=0で主対角、K<0で主対角より下、0<kで主対角より上 (デフォルト値 = 0)
	//! @return	出力ベクトル
	template<size_t M, size_t N, typename T = double, size_t L = std::min(M,N)>
	constexpr ArcsMat<L,1,T> getdiag(const ArcsMat<M,N,T>& U, const ssize_t k = 0){
		return ArcsMat<M,N,T>::getdiag(U, k);
	}
	
	//! @brief 行列のトレースを返す関数(戻り値渡し版のみ)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	結果
	template<size_t M, size_t N, typename T = double>
	constexpr T trace(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::trace(U);
	}
	
	//! @brief 行列の対角要素の総積を返す関数(戻り値渡し版のみ)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	結果
	template<size_t M, size_t N, typename T = double>
	constexpr T multdiag(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::multdiag(U);
	}
	
	//! @brief 行列要素の最大値の要素番号を返す関数(戻り値渡し版のみ)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	結果
	template<size_t M, size_t N, typename T = double>
	constexpr std::tuple<size_t,size_t> maxidx(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::maxidx(U);
	}
	
	//! @brief 行列要素の最大値を返す関数(戻り値渡し版のみ)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	結果
	template<size_t M, size_t N, typename T = double>
	constexpr T max(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::max(U);
	}
	
	//! @brief 行列要素の最小値の要素番号を返す関数(戻り値渡し版のみ)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	結果
	template<size_t M, size_t N, typename T = double>
	constexpr std::tuple<size_t,size_t> minidx(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::minidx(U);
	}
	
	//! @brief 行列要素の最小値を返す関数(戻り値渡し版のみ)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	結果
	template<size_t M, size_t N, typename T = double>
	constexpr T min(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::min(U);
	}

	//! @brief 行列要素の総和を返す関数(戻り値渡し版のみ)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	結果
	template<size_t M, size_t N, typename T = double>
	constexpr T sum(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::sum(U);
	}

	//! @brief 行列要素の指数関数を計算する関数(引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力ベクトルと出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	出力行列
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void exp(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
		ArcsMat<M,N,T>::exp(U, Y);
	}

	//! @brief 行列要素の指数関数を計算する関数(戻り値渡し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	Y	出力行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> exp(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::exp(U);
	}

	//! @brief 行列要素の対数関数(底e版)を計算する関数(引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力ベクトルと出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	出力行列
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void log(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
		ArcsMat<M,N,T>::log(U, Y);
	}

	//! @brief 行列要素の対数関数(底e版)を計算する関数(戻り値渡し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	Y	出力行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> log(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::log(U);
	}

	//! @brief 行列要素の対数関数(底10版)を計算する関数(引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力ベクトルと出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	出力行列
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void log10(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
		ArcsMat<M,N,T>::log10(U, Y);
	}

	//! @brief 行列要素の対数関数(底10版)を計算する関数(戻り値渡し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	Y	出力行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> log10(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::log10(U);
	}

	//! @brief 行列要素の正弦関数を計算する関数(引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力ベクトルと出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	出力行列
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void sin(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
		ArcsMat<M,N,T>::sin(U, Y);
	}

	//! @brief 行列要素の正弦関数を計算する関数(戻り値渡し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	Y	出力行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> sin(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::sin(U);
	}

	//! @brief 行列要素の余弦関数を計算する関数(引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力ベクトルと出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	出力行列
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void cos(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
		ArcsMat<M,N,T>::cos(U, Y);
	}

	//! @brief 行列要素の余弦関数を計算する関数(戻り値渡し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	Y	出力行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> cos(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::cos(U);
	}

	//! @brief 行列要素の正接関数を計算する関数(引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力ベクトルと出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	出力行列
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void tan(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
		ArcsMat<M,N,T>::tan(U, Y);
	}

	//! @brief 行列要素の正接関数を計算する関数(戻り値渡し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	Y	出力行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> tan(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::tan(U);
	}

	//! @brief 行列要素の双曲線正接関数を計算する関数(引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力ベクトルと出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	出力行列
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void tanh(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
		ArcsMat<M,N,T>::tanh(U, Y);
	}

	//! @brief 行列要素の双曲線正接関数を計算する関数(戻り値渡し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	Y	出力行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> tanh(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::tanh(U);
	}

	//! @brief 行列要素の平方根を計算する関数(引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力ベクトルと出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	出力行列
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void sqrt(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
		ArcsMat<M,N,T>::sqrt(U, Y);
	}

	//! @brief 行列要素の平方根を計算する関数(戻り値渡し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	Y	出力行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> sqrt(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::sqrt(U);
	}

	//! @brief 行列要素の符号関数を計算する関数(引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力ベクトルと出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	出力行列
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void sign(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
		ArcsMat<M,N,T>::sign(U, Y);
	}

	//! @brief 行列要素の符号関数を計算する関数(戻り値渡し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	Y	出力行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> sign(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::sign(U);
	}

	//! @brief 行列要素の絶対値を計算する関数(引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列と出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	出力行列
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void abs(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
		ArcsMat<M,N,T>::abs(U, Y);
	}

	//! @brief 行列要素の絶対値を計算する関数(戻り値渡し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	Y	出力行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> abs(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::abs(U);
	}

	//! @brief 行列要素の絶対値を計算する関数(戻り値渡し版, 複素数版特殊化)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	Y	出力行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> abs(const ArcsMat<M,N,std::complex<T>>& U){
		ArcsMat<M,N,T> Y;							// 実数返し用
		ArcsMat<M,N,std::complex<T>>::abs(U, Y);	// 入力は複素数、出力は実数
		return Y;	// 実数で返す
	}

	//! @brief 複素数行列要素の偏角を計算する関数(引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列と出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	出力行列
	template<size_t M, size_t N, typename T = std::complex<double>, size_t P, size_t Q, typename R = double>
	constexpr void arg(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
		ArcsMat<M,N,T>::arg(U, Y);
	}

	//! @brief 複素数行列要素の偏角を計算する関数(戻り値渡し版, 複素数版特殊化)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	Y	出力行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> arg(const ArcsMat<M,N,std::complex<T>>& U){
		ArcsMat<M,N,T> Y;							// 実数返し用
		ArcsMat<M,N,std::complex<T>>::arg(U, Y);	// 入力は複素数、出力は実数
		return Y;	// 実数で返す
	}

	//! @brief 複素数行列要素の実数部を取得する関数(引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列と出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	出力行列
	template<size_t M, size_t N, typename T = std::complex<double>, size_t P, size_t Q, typename R = double>
	constexpr void real(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
		ArcsMat<M,N,T>::real(U, Y);
	}

	//! @brief 複素数行列要素の実数部を計算する関数(戻り値渡し版, 複素数版特殊化)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	Y	出力行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> real(const ArcsMat<M,N,std::complex<T>>& U){
		ArcsMat<M,N,T> Y;							// 実数返し用
		ArcsMat<M,N,std::complex<T>>::real(U, Y);	// 入力は複素数、出力は実数
		return Y;	// 実数で返す
	}

	//! @brief 複素数行列要素の虚数部を取得する関数(引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列と出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	出力行列
	template<size_t M, size_t N, typename T = std::complex<double>, size_t P, size_t Q, typename R = double>
	constexpr void imag(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
		ArcsMat<M,N,T>::imag(U, Y);
	}

	//! @brief 複素数行列要素の虚数部を計算する関数(戻り値渡し版, 複素数版特殊化)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	Y	出力行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> imag(const ArcsMat<M,N,std::complex<T>>& U){
		ArcsMat<M,N,T> Y;							// 実数返し用
		ArcsMat<M,N,std::complex<T>>::imag(U, Y);	// 入力は複素数、出力は実数
		return Y;	// 実数で返す
	}

	//! @brief 複素数行列要素の複素共役を取得する関数(引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列と出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	出力行列
	template<size_t M, size_t N, typename T = std::complex<double>, size_t P, size_t Q, typename R = std::complex<double>>
	constexpr void conj(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
		ArcsMat<M,N,T>::conj(U, Y);
	}

	//! @brief 複素数行列要素の複素共役を計算する関数(戻り値渡し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	Y	出力行列
	template<size_t M, size_t N, typename T = std::complex<double>>
	constexpr ArcsMat<M,N,T> conj(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::conj(U);
	}
	
	//! @brief 転置行列を返す関数 (引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列の高さ, 幅, 要素の型, 結果の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	結果
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void tp(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
		ArcsMat<M,N,T>::tp(U, Y);
	}
	
	//! @brief 転置行列を返す関数 (戻り値渡し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	Y	出力行列
	template <size_t M, size_t N, typename T = double>
	constexpr ArcsMat<N,M,T> tp(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::tp(U);
	}
		
	//! @brief エルミート転置行列を返す関数 (引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	出力行列
	template<size_t M, size_t N, typename T = std::complex<double>, size_t P, size_t Q, typename R = std::complex<double>>
	constexpr void Htp(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
		ArcsMat<M,N,T>::Htp(U, Y);
	}
	
	//! @brief エルミート転置行列を返す関数 (戻り値渡し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	出力行列
	template<size_t M, size_t N, typename T = std::complex<double>>
	constexpr ArcsMat<N,M,T> Htp(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::Htp(U);
	}

	//! @brief 行列のノルムを返す関数(戻り値渡し版のみ)
	//! @tparam	NRM, M, N, T	ノルムのタイプ, 入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	結果
	template<NormType NRM = ArcsMatrix::NormType::AMT_L2, typename T = double, size_t M, size_t N>
	constexpr T norm(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::template norm<NRM>(U);
	}
	
	//! @brief 行列のノルムを返す関数(戻り値渡し版のみ, 複素数版特殊化)
	//! @tparam	NRM, M, N, T	ノルムのタイプ, 入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	結果
	template<NormType NRM = ArcsMatrix::NormType::AMT_L2, typename R = double, typename T = double, size_t M, size_t N>
	constexpr R norm(const ArcsMat<M,N,std::complex<T>>& U){
		return ArcsMat<M,N,std::complex<T>>::template norm<NRM, R>(U);
	}
	
	//! @brief n列目を左端として右上の上三角部分のみを返す関数(下三角部分はゼロ)(引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列の高さ, 幅, 要素の型, 出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	出力行列
	//! @param[in]	n	切り出す左端位置 n列目(デフォルト値 = 1)
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void gettriup(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y, const size_t n = 1){
		ArcsMat<M,N,T>::gettriup(U, Y, n);
	}

	//! @brief n列目を左端として右上の上三角部分のみを返す関数(下三角部分はゼロ)(戻り値渡し版)
	//! @tparam	M, N, T	入出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[in]	n	切り出す左端位置 n列目(デフォルト値 = 1)
	//! @return	出力行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> gettriup(const ArcsMat<M,N,T>& U, const size_t n = 1){
		return ArcsMat<M,N,T>::gettriup(U, n);
	}

	//! @brief m行目を上端として左下の下三角部分のみを返す関数(上三角部分はゼロ)(引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入力行列の高さ, 幅, 要素の型, 出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	出力行列
	//! @param[in]	m	切り出す上端位置 m行目(デフォルト値 = 1)
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void gettrilo(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y, const size_t m = 1){
		ArcsMat<M,N,T>::gettrilo(U, Y, m);
	}

	//! @brief m行目を上端として左下の下三角部分のみを返す関数(上三角部分はゼロ)(戻り値渡し版)
	//! @tparam	M, N, T	入出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[in]	m	切り出す上端位置 m行目(デフォルト値 = 1)
	//! @return	出力行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> gettrilo(const ArcsMat<M,N,T>& U, const size_t m = 1){
		return ArcsMat<M,N,T>::gettrilo(U, m);
	}

	//! @brief LU分解の結果と置換行列を返す関数(引数渡し版)
	//! @tparam	M, N, T, ML, NL, TL, MU, NU, TU, MP, NP, TP	入出力行列,L,U,P行列の高さ, 幅, 要素の型
	//! @param[in]	A	入力行列
	//! @param[out]	L	下三角行列
	//! @param[out]	U	上三角行列
	//! @param[out]	P	置換行列
	template<
		size_t M, size_t N, typename T = double, size_t ML, size_t NL, typename TL = double,
		size_t MU, size_t NU, typename TU = double, size_t MP, size_t NP, typename TP = double
	>
	constexpr void LUP(const ArcsMat<M,N,T>& A, ArcsMat<ML,NL,TL>& L, ArcsMat<MU,NU,TU>& U, ArcsMat<MP,NP,TP>& P){
		ArcsMat<M,N,T>::LUP(A, L, U, P);
	}

	//! @brief LU分解の結果と置換行列を返す関数(タプル返し版)
	//! @tparam	M, N, T	入出力行列の高さ, 幅, 要素の型
	//! @param[in]	A	入力行列
	//! @return	(L, U, P)	(下三角行列, 上三角行列, 置換行列)のタプル
	template<size_t M, size_t N, typename T = double>
	constexpr std::tuple<ArcsMat<M,N,T>, ArcsMat<M,N,T>, ArcsMat<M,N,T>> LUP(const ArcsMat<M,N,T>& A){
		return ArcsMat<M,N,T>::LUP(A);
	}
	
	//! @brief LU分解の結果のみ返す関数(引数渡し版)
	//! @tparam	M, N, T, ML, NL, TL, MU, NU, TU	入出力行列,L,U行列の高さ, 幅, 要素の型
	//! @param[in]	A	入力行列
	//! @param[out]	L	下三角行列
	//! @param[out]	U	上三角行列
	template<
		size_t M, size_t N, typename T = double,
		size_t ML, size_t NL, typename TL = double, size_t MU, size_t NU, typename TU = double
	>
	constexpr void LU(const ArcsMat<M,N,T>& A, ArcsMat<ML,NL,TL>& L, ArcsMat<MU,NU,TU>& U){
		ArcsMat<M,N,T>::LU(A, L, U);
	}

	//! @brief LU分解の結果のみ返す関数(タプル返し版)
	//! @tparam	M, N, T	入出力行列の高さ, 幅, 要素の型
	//! @param[in]	A	入力行列
	//! @return	(L, U)	(下三角行列, 上三角行列)のタプル
	template<size_t M, size_t N, typename T = double>
	constexpr std::tuple<ArcsMat<M,N,T>, ArcsMat<M,N,T>> LU(const ArcsMat<M,N,T>& A){
		return ArcsMat<M,N,T>::LU(A);
	}

	//! @brief 行列式の値を返す関数(戻り値返し版のみ)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	A	入力行列
	//! @return	結果
	template<size_t M, size_t N, typename T = double>
	constexpr T det(const ArcsMat<M,N,T>& A){
		return ArcsMat<M,N,T>::det(A);
	}
	
	//! @brief Householder行列を生成する関数(引数渡し版)
	//! @tparam	M, N, T, MH, NH, TH	入力ベクトルとH行列の高さ, 幅, 要素の型
	//! @param[in]	v	入力縦ベクトル
	//! @param[out]	H	ハウスホルダー行列
	//! @param[in]	k	次元 (デフォルト値 = 1)
	template<size_t M, size_t N, typename T = double, size_t MH, size_t NH, typename TH = double>
	constexpr void Householder(const ArcsMat<M,N,T>& v, ArcsMat<MH,NH,TH>& H, const size_t k = 1){
		ArcsMat<M,N,T>::Householder(v, H, k);
	}
	
	//! @brief Householder行列を生成する関数(戻り値返し版)
	//! @tparam	M, N, T	入力ベクトルの高さ, 幅, 要素の型
	//! @param[in]	v	入力縦ベクトル
	//! @param[in]	k	次元 (デフォルト値 = 1)
	//! @return	ハウスホルダー行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,M,T> Householder(const ArcsMat<M,N,T>& v, const size_t k = 1){
		return ArcsMat<M,N,T>::Householder(v, k);
	}
	
	//! @brief QR分解(引数渡し版)
	//! @tparam	M, N, T, MQ, NQ, TQ, MQ, NQ, TQ	入力行列,Q,R行列の高さ, 幅, 要素の型
	//! @param[in]	A	入力行列
	//! @param[out]	Q	ユニタリ行列 Q行列
	//! @param[out] R	上三角行列 R行列
	template<
		size_t M, size_t N, typename T = double,
		size_t MQ, size_t NQ, typename TQ = double, size_t MR, size_t NR, typename TR = double
	>
	constexpr void QR(const ArcsMat<M,N,T>& A, ArcsMat<MQ,NQ,TQ>& Q, ArcsMat<MR,NR,TR>& R){
		ArcsMat<M,N,T>::QR(A, Q, R);
	}

	//! @brief QR分解(タプル返し版)
	//! @tparam	M, N, T	入出力行列の高さ, 幅, 要素の型
	//! @param[in]	A	入力行列
	//! @return	(Q, R)	(ユニタリ行列, 上三角行列)のタプル
	template<size_t M, size_t N, typename T = double>
	constexpr std::tuple<ArcsMat<M,M,T>, ArcsMat<M,N,T>> QR(const ArcsMat<M,N,T>& A){
		return ArcsMat<M,N,T>::QR(A);
	}

	//! @brief SVD特異値分解(引数渡し版)
	//! @tparam	M, N, T, MU, NU, TU, MS, NS, TS, MV, NV, TV	入力行列,U,Σ,V行列の高さ, 幅, 要素の型
	//! @param[in]	A	入力行列
	//! @param[out]	U	U行列
	//! @param[out]	S	Σ行列
	//! @param[out]	V	V行列
	template<
		size_t M, size_t N, typename T = double, size_t MU, size_t NU, typename TU = double,
		size_t MS, size_t NS, typename TS = double, size_t MV, size_t NV, typename TV = double
	>
	constexpr void SVD(const ArcsMat<M,N,T>& A, ArcsMat<MU,NU,TU>& U, ArcsMat<MS,NS,TS>& S, ArcsMat<MV,NV,TV>& V){
		ArcsMat<M,N,T>::SVD(A, U, S, V);
	}

	//! @brief SVD特異値分解(タプル返し版)
	//! @tparam	M, N, T	入出力行列の高さ, 幅, 要素の型
	//! @param[in]	A	入力行列
	//! @return	(U, S, V)	(U行列, Σ行列, V行列)のタプル
	template<size_t M, size_t N, typename T = double>
	constexpr std::tuple<ArcsMat<M,M,T>, ArcsMat<M,N,T>, ArcsMat<N,N,T>> SVD(const ArcsMat<M,N,T>& A){
		return ArcsMat<M,N,T>::SVD(A);
	}
	
	//! @brief 行列の階数を返す関数(戻り値返し版のみ)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	A	入力行列
	//! @return	結果
	template<size_t M, size_t N, typename T = double>
	constexpr size_t rank(const ArcsMat<M,N,T>& A, const T eps = ArcsMatrix::EPSILON){
		return ArcsMat<M,N,T>::rank(A, eps);
	}

	//! @brief 修正コレスキー分解(LDL分解) (引数渡し版)
	//! @tparam	M, N, T, ML, NL, TL, MD, ND, TD	A, L, D行列の高さ, 幅, 要素の型
	//! @param[in]	A	入力行列
	//! @param[out]	L	下三角行列
	//! @param[out]	D	対角行列
	template<size_t M, size_t N, typename T = double, size_t ML, size_t NL, typename TL = double, size_t MD, size_t ND, typename TD = double>
	constexpr void LDL(const ArcsMat<M,N,T>& A, ArcsMat<ML,NL,TL>& L, ArcsMat<MD,ND,TD>& D){
		ArcsMat<M,N,T>::LDL(A, L, D);
	}

	//! @brief 修正コレスキー分解(LDL分解) (タプル返し版)
	//! @tparam	M, N, T	入出力行列の高さ, 幅, 要素の型
	//! @param[in]	A	入力行列
	//! @return	(L, D)	(L行列, D行列)のタプル
	template<size_t M, size_t N, typename T = double>
	constexpr std::tuple<ArcsMat<M,N,T>, ArcsMat<M,N,T>> LDL(const ArcsMat<M,N,T>& A){
		return ArcsMat<M,N,T>::LDL(A);
	}

	//! @brief 修正コレスキー分解(引数渡し版)
	//! @tparam	M, N, T, ML, NL, TL	入力出力行列の高さ, 幅, 要素の型
	//! @param[in]	A	入力行列
	//! @param[out]	R	上三角行列
	template<size_t M, size_t N, typename T = double, size_t ML, size_t NL, typename TL = double>
	constexpr void Cholesky(const ArcsMat<M,N,T>& A, ArcsMat<ML,NL,TL>& R){
		ArcsMat<M,N,T>::Cholesky(A, R);
	}

	//! @brief 修正コレスキー分解(戻り値返し版)
	//! @tparam	M, N, T	入出力行列の高さ, 幅, 要素の型
	//! @param[in]	A	入力行列
	//! @return	上三角行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> Cholesky(const ArcsMat<M,N,T>& A){
		return ArcsMat<M,N,T>::Cholesky(A);
	}

	//! @brief AX = Bの形の線形方程式をXについて解く関数(引数渡し版)
	//! @tparam	M, N, T, MB, NB, TB, MX, NX, TX	A、BとXの高さ, 幅, 要素の型
	//! @param[in]	A	係数行列(正方行列・非正方行列)
	//! @param[in]	B	係数ベクトル・行列
	//! @param[out]	X	解ベクトル・行列
	template<size_t M, size_t N, typename T = double, size_t MB, size_t NB, typename TB = double, size_t MX, size_t NX, typename TX = double>
	constexpr void linsolve(const ArcsMat<M,N,T>& A, const ArcsMat<MB,NB,TB>& B, ArcsMat<MX,NX,TX>& X){
		ArcsMat<M,N,T>::linsolve(A, B, X);
	}
	
	//! @brief AX = Bの形の線形方程式をXについて解く関数(戻り値返し版)
	//! @tparam	M, N, T, MB, NB, TB	A、Bの高さ, 幅, 要素の型
	//! @param[in]	A	係数行列(正方行列・非正方行列)
	//! @param[in]	B	係数ベクトル・行列
	//! @return	X	解ベクトル・行列
	template<size_t M, size_t N, typename T = double, size_t MB, size_t NB, typename TB = double>
	constexpr ArcsMat<N,NB,T> linsolve(const ArcsMat<M,N,T>& A, const ArcsMat<MB,NB,TB>& B){
		return ArcsMat<M,N,T>::linsolve(A, B);
	}

	//! @brief XA = Bの形の線形方程式をXについて解く関数(引数渡し版)
	//! @tparam	M, N, T, MB, NB, TB, MX, NX, TX	A、BとXの高さ, 幅, 要素の型
	//! @param[in]	A	係数行列(正方行列・非正方行列)
	//! @param[in]	B	係数ベクトル・行列
	//! @param[out]	X	解ベクトル・行列
	template<size_t M, size_t N, typename T = double, size_t MB, size_t NB, typename TB = double, size_t MX, size_t NX, typename TX = double>
	constexpr void linsolveXAB(const ArcsMat<M,N,T>& A, const ArcsMat<MB,NB,TB>& B, ArcsMat<MX,NX,TX>& X){
		ArcsMat<M,N,T>::linsolveXAB(A, B, X);
	}
	
	//! @brief AX = Bの形の線形方程式をXについて解く関数(戻り値返し版)
	//! @tparam	M, N, T, MB, NB, TB	A、Bの高さ, 幅, 要素の型
	//! @param[in]	A	係数行列(正方行列・非正方行列)
	//! @param[in]	B	係数ベクトル・行列
	//! @return	X	解ベクトル・行列
	template<size_t M, size_t N, typename T = double, size_t MB, size_t NB, typename TB = double>
	constexpr ArcsMat<MB,M,T> linsolveXAB(const ArcsMat<M,N,T>& A, const ArcsMat<MB,NB,TB>& B){
		return ArcsMat<M,N,T>::linsolveXAB(A, B);
	}

	//! @brief 逆行列を返す関数(引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入出力行列の高さ, 幅, 要素の型
	//! @param[in]	A	入力行列
	//! @param[out]	Y	出力行列
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void inv(const ArcsMat<M,N,T>& A, ArcsMat<P,Q,R>& Y){
		ArcsMat<M,N,T>::inv(A, Y);
	}

	//! @brief 逆行列を返す関数(戻り値返し版)
	//! @tparam	M, N, T	入出力行列の高さ, 幅, 要素の型
	//! @param[in]	A	入力行列
	//! @return	出力行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> inv(const ArcsMat<M,N,T>& A){
		return ArcsMat<M,N,T>::inv(A);
	}

	//! @brief Moore-Penroseの擬似逆行列を返す関数(引数渡し版)
	//! @tparam	M, N, T, P, Q, R	入出力行列の高さ, 幅, 要素の型
	//! @param[in]	A	入力行列
	//! @param[out]	Y	出力行列
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double>
	constexpr void pinv(const ArcsMat<M,N,T>& A, ArcsMat<P,Q,R>& Y){
		ArcsMat<M,N,T>::pinv(A, Y);
	}

	//! @brief Moore-Penroseの疑似逆行列を返す関数(戻り値返し版)
	//! @tparam	M, N, T	入出力行列の高さ, 幅, 要素の型
	//! @param[in]	A	入力行列
	//! @return	出力行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<N,M,T> pinv(const ArcsMat<M,N,T>& A){
		return ArcsMat<M,N,T>::pinv(A);
	}
	
	//! @brief Hessenberg分解(引数渡し版)
	//!        複素数の場合、この関数はMATLABとは異なる解を出力する。
	//! @tparam	M, N, T, MP, NP, TP, MH, NH, TH	出力行列の高さ, 幅, 要素の型
	//! @param[in]	A	入力行列
	//! @param[out]	P	ユニタリ行列
	//! @param[out]	H	ヘッセンベルグ行列
	template<size_t M, size_t N, typename T = double, size_t MP, size_t NP, typename TP = double, size_t MH, size_t NH, typename TH = double>
	constexpr void Hessenberg(const ArcsMat<M,N,T>& A, ArcsMat<MP,NP,TP>& P, ArcsMat<MH,NH,TH>& H){
		ArcsMat<M,N,T>::Hessenberg(A, P, H);
	}

	//! @brief Hessenberg分解(タプル返し版)
	//!        複素数の場合、この関数はMATLABとは異なる解を出力する。
	//! @tparam	M, N, T	入出力行列の高さ, 幅, 要素の型
	//! @param[in]	A	入力行列
	//! @return	(ユニタリ行列P, ヘッセンベルグ行列H) のタプル
	template<size_t M, size_t N, typename T = double>
	constexpr std::tuple<ArcsMat<M,N,T>, ArcsMat<M,N,T>> Hessenberg(const ArcsMat<M,N,T>& A){
		return ArcsMat<M,N,T>::Hessenberg(A);
	}

	//! @brief 複素Schur分解(引数渡し版)
	//!        この関数はMATLABとは異なる解を出力する、ただしもちろん、A = USU' は成立
	//! @tparam	M, N, T, MU, NU, TU, MS, NS, TS 	入出力行列の高さ, 幅, 要素の型
	//! @param[in]	A	入力行列
	//! @param[out]	U	ユニタリ行列
	//! @param[out]	S	上三角行列または疑似上三角行列
	template<size_t M, size_t N, typename T = double, size_t MU, size_t NU, typename TU = double, size_t MS, size_t NS, typename TS = double>
	constexpr void Schur(const ArcsMat<M,N,T>& A, ArcsMat<MU,NU,TU>& U, ArcsMat<MS,NS,TS>& S){
		ArcsMat<M,N,T>::Schur(A, U, S);
	}

	//! @brief Schur分解(タプル返し版)
	//!        この関数はMATLABとは異なる解を出力する、ただしもちろん、A = USU' は成立
	//! @tparam	M, N, T	入出力行列の高さ, 幅, 要素の型
	//! @param[in]	A	入力行列
	//! @return	(ユニタリ行列U, 上三角行列または疑似上三角行列S)のタプル
	template<size_t M, size_t N, typename T = double>
	constexpr std::tuple<ArcsMat<M,N,T>, ArcsMat<M,N,T>> Schur(const ArcsMat<M,N,T>& A){
		return ArcsMat<M,N,T>::Schur(A);
	}

	//! @brief 固有値を返す関数
	//! @tparam	M, N, T, MV, NV, TV	入出力行列の高さ, 幅, 要素の型
	//! @param[in]	A	入力行列
	//! @param[out]	v	固有値の縦ベクトル
	template<size_t M, size_t N, typename T = double, size_t MV, size_t NV, typename TV = std::complex<double>>
	constexpr void eig(const ArcsMat<M,N,T>& A, ArcsMat<MV,NV,TV>& v){
		ArcsMat<M,N,T>::eig(A, v);
	}

	//! @brief 固有値を返す関数(戻り値返し版)
	//! @tparam	M, N, T, TV	入力行列の高さ, 幅, 要素の型, 出力行列の要素の型
	//! @param[in]	A	入力行列
	//! @return	固有値の縦ベクトル
	template<size_t M, size_t N, typename T = double, typename TV = std::complex<double>>
	constexpr ArcsMat<M,1,TV> eig(const ArcsMat<M,N,T>& A){
		return ArcsMat<M,N,T>::eig(A);
	}

	//! @brief 固有値を持つ対角行列Dと、各々の固有値に対応する固有ベクトルを持つ行列Vを返す関数(引数渡し版)
	//! @tparam	M, N, T, MV, NV, TV, MD, ND, TD	入出力行列の高さ, 幅, 要素の型
	//! @param[in]	A	入力行列
	//! @param[out]	V	各々の固有ベクトルを縦ベクトルとして持つ行列V
	//! @param[out]	D	固有値を対角に持つ対角行列D
	template<size_t M, size_t N, typename T = double, size_t MV, size_t NV, typename TV = std::complex<double>, size_t MD, size_t ND, typename TD = std::complex<double>>
	constexpr void eigvec(const ArcsMat<M,N,T>& A, ArcsMat<MV,NV,TV>& V, ArcsMat<MD,ND,TD>& D){
		ArcsMat<M,N,T>::eigvec(A, V, D);
	}

	//! @brief 固有値を持つ対角行列Dと、各々の固有値に対応する固有ベクトルを持つ行列Vを返す関数(タプル返し版)
	//! @tparam	M, N, T, TV	入出力行列の高さ, 幅, 要素の型
	//! @param[in]	A	入力行列
	//! @return	(V, D)	(各々の固有ベクトルを縦ベクトルとして持つ行列V, 固有値を対角に持つ対角行列D) のタプル
	template<size_t M, size_t N, typename T = double, typename TV = std::complex<double>>
	constexpr std::tuple< ArcsMat<M,N,TV>, ArcsMat<M,N,TV> > eigvec(const ArcsMat<M,N,T>& A){
		return ArcsMat<M,N,T>::eigvec(A);
	}

	//! @brief クロネッカー積(引数渡し版)
	//! @tparam	M, N, T, MR, NR, TR, MY, NY, TY	入出力行列の高さ, 幅, 要素の型
	//! @param[in] L 演算子の左側
	//! @param[in] R 演算子の右側
	//! @param[out]	Y 結果
	template<size_t M, size_t N, typename T = double, size_t MR, size_t NR, typename TR = double, size_t MY, size_t NY, typename TY = double>
	constexpr void Kron(const ArcsMat<M,N,T>& L, const ArcsMat<MR,NR,TR>& R, ArcsMat<MY,NY,TY>& Y){
		ArcsMat<M,N,T>::Kron(L, R, Y);
	}

	//! @brief クロネッカー積(戻り値返し版)
	//! @tparam	M, N, T, MR, NR, TR	入力行列の高さ, 幅, 要素の型
	//! @param[in] L 演算子の左側
	//! @param[in] R 演算子の右側
	//! @return 結果
	template<size_t M, size_t N, typename T = double, size_t MR, size_t NR, typename TR = double>
	constexpr ArcsMat<M*MR,N*NR,T> Kron(const ArcsMat<M,N,T>& L, const ArcsMat<MR,NR,TR>& R){
		return ArcsMat<M,N,T>::Kron(L, R);
	}

	//! @brief クロス積 (引数渡し版)
	//! @tparam	M, N, T, MR, NR, TR, MY, NY, TY	入出力行列の高さ, 幅, 要素の型
	//! @param[in] L 演算子の左側
	//! @param[in] R 演算子の右側
	//! @param[out]	Y 結果
	template<size_t M, size_t N, typename T = double, size_t MR, size_t NR, typename TR = double, size_t MY, size_t NY, typename TY = double>
	constexpr void cross(const ArcsMat<M,N,T>& L, const ArcsMat<MR,NR,TR>& R, ArcsMat<MY,NY,TY>& Y){
		ArcsMat<M,N,T>::cross(L, R, Y);
	}

	//! @brief クロス積 (戻り値返し版)
	//! @tparam	M, N, T, MR, NR, TR 入力行列の高さ, 幅, 要素の型
	//! @param[in] L 演算子の左側
	//! @param[in] R 演算子の右側
	//! @return 結果
	template<size_t M, size_t N, typename T = double, size_t MR, size_t NR, typename TR = double>
	constexpr ArcsMat<M,N,T> cross(const ArcsMat<M,N,T>& L, const ArcsMat<MR,NR,TR>& R){
		return ArcsMat<M,N,T>::cross(L, R);
	}

	//! @brief vec作用素(行列→縦ベクトル) (引数渡し版)
	//! @tparam	M, N, T, MY, NY, TY 入力行列と出力ベクトルの高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	y	出力ベクトル
	template<size_t M, size_t N, typename T = double, size_t MY, size_t NY, typename TY = double>
	constexpr void vec(const ArcsMat<M,N,T>& U, ArcsMat<MY,NY,TY>& y){
		ArcsMat<M,N,T>::vec(U, y);
	}

	//! @brief vec作用素(行列→縦ベクトル) (戻り値返し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	y	出力ベクトル
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M*N,1,T> vec(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::vec(U);
	}

	//! @brief vec作用素の逆(縦ベクトル→行列) (引数渡し版)
	//! @tparam	M, N, T, MY, NY, TY 入力縦ベクトルと出力行列の高さ, 幅, 要素の型
	//! @param[in]	u	入力縦ベクトル
	//! @param[out]	Y	再構成後の行列
	template<size_t M, size_t N, typename T = double, size_t MY, size_t NY, typename TY = double>
	constexpr void vecinv(const ArcsMat<M,N,T>& u, ArcsMat<MY,NY,TY>& Y){
		ArcsMat<M,N,T>::vecinv(u, Y);
	}

	//! @brief vec作用素の逆(縦ベクトル→行列) (戻り値返し版)
	//! @tparam	MY, NY	出力行列の高さ, 幅
	//! @tparam M, N, T	入力縦ベクトルの高さ, 幅, 要素の型
	//! @param[in]	u	入力縦ベクトル
	//! @return	再構成後の行列
	template<size_t MY, size_t NY, size_t M, size_t N, typename T = double>
	constexpr ArcsMat<MY,NY,T> vecinv(const ArcsMat<M,N,T>& u){
		return ArcsMat<M,N,T>::template vecinv<MY,NY>(u);
	}

	//! @brief 行列指数関数 (引数渡し版)
	//! @tparam	K	パデ近似の次数 (デフォルト値 = 13)
	//! @tparam	M, N, T, MY, NY, TY 入出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	Y	出力行列
	template<size_t K = 13, size_t M, size_t N, typename T = double, size_t MY, size_t NY, typename TY = double>
	constexpr void expm(const ArcsMat<M,N,T>& U, ArcsMat<MY,NY,TY>& Y){
		ArcsMat<M,N,T>::template expm<K>(U, Y);
	}
		
	//! @brief 行列指数関数 (戻り値返し版)
	//! @tparam	K	パデ近似の次数 (デフォルト値 = 13)
	//! @tparam	M, N, T	入出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	出力行列
	template<size_t K = 13, size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,N,T> expm(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::template expm<K>(U);
	}

	//! @brief 縦方向の平均を計算する関数 行列入力-横ベクトル出力版 (引数渡し版)
	//! @tparam	M, N, T, MY, NY, TY 入出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	y	出力ベクトル
	template<size_t M, size_t N, typename T = double, size_t MY, size_t NY, typename TY = double>
	constexpr void meancolumn(const ArcsMat<M,N,T>& U, ArcsMat<MY,NY,TY>& y){
		ArcsMat<M,N,T>::meancolumn(U, y);
	}

	//! @brief 縦方向の平均を計算する関数 行列入力-横ベクトル出力版 (戻り値返し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	出力ベクトル
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<1,N,T> meancolumn(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::meancolumn(U);
	}
	
	//! @brief 横方向の平均を計算する関数 行列入力-縦ベクトル出力版 (引数渡し版)
	//! @tparam	M, N, T, MY, NY, TY 入出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	y	出力ベクトル
	template<size_t M, size_t N, typename T = double, size_t MY, size_t NY, typename TY = double>
	constexpr void meanrow(const ArcsMat<M,N,T>& U, ArcsMat<MY,NY,TY>& y){
		ArcsMat<M,N,T>::meanrow(U, y);
	}

	//! @brief 横方向の平均を計算する関数 行列入力-縦ベクトル出力版 (戻り値返し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	出力ベクトル
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,1,T> meanrow(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::meanrow(U);
	}

	//! @brief ベクトルの平均を計算する関数 ベクトル入力-スカラ出力版 (戻り値返し版のみ)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	u	入力ベクトル
	//! @return	平均値(スカラー)
	template<size_t M, size_t N, typename T = double>
	constexpr T meanvec(const ArcsMat<M,N,T>& u){
		return ArcsMat<M,N,T>::meanvec(u);
	}

	//! @brief 行列全体の平均を計算する関数 (戻り値返し版のみ)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	平均値(スカラー)
	template<size_t M, size_t N, typename T = double>
	constexpr T mean(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::mean(U);
	}

	//! @brief 縦方向の中央値を計算する関数 行列入力-横ベクトル出力版 (引数渡し版)
	//! @tparam	M, N, T, MY, NY, TY 入出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	y	出力ベクトル
	template<size_t M, size_t N, typename T = double, size_t MY, size_t NY, typename TY = double>
	constexpr void mediancolumn(const ArcsMat<M,N,T>& U, ArcsMat<MY,NY,TY>& y){
		ArcsMat<M,N,T>::mediancolumn(U, y);
	}

	//! @brief 縦方向の中央値を計算する関数 行列入力-横ベクトル出力版 (戻り値返し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	出力ベクトル
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<1,N,T> mediancolumn(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::mediancolumn(U);
	}

	//! @brief 横方向の中央値を計算する関数 行列入力-横ベクトル出力版 (引数渡し版)
	//! @tparam	M, N, T, MY, NY, TY 入出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	y	出力ベクトル
	template<size_t M, size_t N, typename T = double, size_t MY, size_t NY, typename TY = double>
	constexpr void medianrow(const ArcsMat<M,N,T>& U, ArcsMat<MY,NY,TY>& y){
		ArcsMat<M,N,T>::medianrow(U, y);
	}

	//! @brief 横方向の中央値を計算する関数 行列入力-横ベクトル出力版 (戻り値返し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	出力ベクトル
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,1,T> medianrow(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::medianrow(U);
	}
	
	//! @brief 行列全体の中央値を計算する関数 (戻り値返し版のみ)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	中央値(スカラー)
	template<size_t M, size_t N, typename T = double>
	constexpr T median(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::median(U);
	}

	//! @brief 縦方向の分散を計算する関数 行列入力-横ベクトル出力版 (引数渡し版)
	//! @tparam	M, N, T, MY, NY, TY 入出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	y	出力ベクトル
	template<size_t M, size_t N, typename T = double, size_t MY, size_t NY, typename TY = double>
	constexpr void varcolumn(const ArcsMat<M,N,T>& U, ArcsMat<MY,NY,TY>& y){
		ArcsMat<M,N,T>::varcolumn(U, y);
	}

	//! @brief 縦方向の分散を計算する関数 行列入力-横ベクトル出力版 (戻り値返し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	出力ベクトル
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<1,N,T> varcolumn(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::varcolumn(U);
	}

	//! @brief 横方向の分散を計算する関数 行列入力-縦ベクトル出力版 (引数渡し版)
	//! @tparam	M, N, T, MY, NY, TY 入出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	y	出力ベクトル
	template<size_t M, size_t N, typename T = double, size_t MY, size_t NY, typename TY = double>
	constexpr void varrow(const ArcsMat<M,N,T>& U, ArcsMat<MY,NY,TY>& y){
		ArcsMat<M,N,T>::varrow(U, y);
	}
	
	//! @brief 横方向の分散を計算する関数 行列入力-縦ベクトル出力版 (戻り値返し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	出力ベクトル
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,1,T> varrow(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::varrow(U);
	}

	//! @brief ベクトルの分散を計算する関数 ベクトル入力-スカラ出力版 (戻り値返し版のみ)
	//! @param[in]	u	入力ベクトル
	//! @return	分散(スカラー)
	template<size_t M, size_t N, typename T = double>
	constexpr T varvec(const ArcsMat<M,N,T>& u){
		return ArcsMat<M,N,T>::varvec(u);
	}

	//! @brief 行列全体の分散を計算する関数 (戻り値返し版のみ)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	分散(スカラー)
	template<size_t M, size_t N, typename T = double>
	constexpr T var(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::var(U);
	}

	//! @brief 縦方向の標準偏差を計算する関数 行列入力-横ベクトル出力版 (引数渡し版)
	//! @tparam	M, N, T, MY, NY, TY 入出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	y	出力ベクトル
	template<size_t M, size_t N, typename T = double, size_t MY, size_t NY, typename TY = double>
	constexpr void stdevcolumn(const ArcsMat<M,N,T>& U, ArcsMat<MY,NY,TY>& y){
		ArcsMat<M,N,T>::stdevcolumn(U, y);
	}

	//! @brief 縦方向の標準偏差を計算する関数 行列入力-横ベクトル出力版 (戻り値返し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	出力ベクトル
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<1,N,T> stdevcolumn(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::stdevcolumn(U);
	}

	//! @brief 横方向の標準偏差を計算する関数 行列入力-縦ベクトル出力版 (引数渡し版)
	//! @tparam	M, N, T, MY, NY, TY 入出力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @param[out]	y	出力ベクトル
	template<size_t M, size_t N, typename T = double, size_t MY, size_t NY, typename TY = double>
	constexpr void stdevrow(const ArcsMat<M,N,T>& U, ArcsMat<MY,NY,TY>& y){
		ArcsMat<M,N,T>::stdevrow(U, y);
	}

	//! @brief 横方向の標準偏差を計算する関数 行列入力-縦ベクトル出力版 (戻り値返し版)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	出力ベクトル
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,1,T> stdevrow(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::stdevrow(U);
	}

	//! @brief 行列全体の標準偏差を計算する関数 (戻り値返し版のみ)
	//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
	//! @param[in]	U	入力行列
	//! @return	分散(スカラー)
	template<size_t M, size_t N, typename T = double>
	constexpr T stdev(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::stdev(U);
	}

	//! @brief 対称ハンケル行列を作成する関数 (引数渡し版)
	//! @tparam	M, N, T	入力ベクトルの高さ, 幅, 要素の型
	//! @tparam	MY, NY, TY 出力行列の高さ, 幅, 要素の型
	//! @param	u	入力ベクトル
	//! @param	Y	出力行列
	template<size_t M, size_t N, typename T = double, size_t MY, size_t NY, typename TY = double>
	constexpr void hankel(const ArcsMat<M,N,T>& u, ArcsMat<MY,NY,TY>& Y){
		ArcsMat<M,N,T>::hankel(u, Y);
	}

	//! @brief 対称ハンケル行列を作成する関数 (戻り値返し版)
	//! @tparam	M, N, T	入力ベクトルの高さ, 幅, 要素の型
	//! @param	u	入力ベクトル
	//! @return	出力行列
	template<size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M,M,T> hankel(const ArcsMat<M,N,T>& u){
		return ArcsMat<M,N,T>::hankel(u);
	}

	//! @brief 根(極)から多項式の係数を計算する関数 (引数渡し版)
	//! MATLABでいうところのpoly関数
	//! 参考文献： Calif. Inst. of Technology, "Compute Polynomial Coeffcients from Roots," Jul. 2015.
	//!  (z - z1)*(z - z2)*...*(z - zn)
	//! 上記の根z1～znから、下記の係数c(1)～c(n+1) を求める。
	//!  c(1)*z^n + ... + c(n)*z + c(n+1)
	//! 使い方の例：
	//!  z[1] = 2 + 1i;
	//!  z[2] = 3 + 2i;
	//!  polycoeff(z, c)
	//!  計算結果： c = {1.0000 + 0.0000i, -5.0000 - 3.0000i, 4.0000 + 7.0000i}
	//! @tparam	M, N, T	入力ベクトルの高さ, 幅, 要素の型
	//! @tparam	MY, NY, TY 出力ベクトルの高さ, 幅, 要素の型
	//!	@param	u	根(極)の値が羅列された入力縦ベクトル
	//! @param	y	多項式(M次方程式)の係数が羅列された出力縦ベクトル
	//! @param	Tol	許容誤差(デフォルト値 1e-10)
	template<size_t M, size_t N, typename T = double, size_t MY, size_t NY, typename TY = double>
	constexpr void polycoeff(const ArcsMat<M,N,T>& u, ArcsMat<MY,NY,TY>& y, const T Tol = 1e-10){
		ArcsMat<M,N,T>::template polycoeff<MY,NY,TY>(u, y, Tol);
	}

	//! @brief 根(極)から多項式の係数を計算する関数 (戻り値返し版)
	//! MATLABでいうところのpoly関数
	//! @tparam	M, N, T	入力ベクトルの高さ, 幅, 要素の型
	//!	@param	u	根(極)の値が羅列された入力縦ベクトル
	//! @param	Tol	許容誤差(デフォルト値 1e-10)
	//! @return	多項式(M次方程式)の係数が羅列された出力縦ベクトル
	template<typename TY = double, size_t M, size_t N, typename T = double>
	constexpr ArcsMat<M+1,1,TY> polycoeff(const ArcsMat<M,N,T>& u, const T Tol = 1e-10){
		ArcsMat<M+1,1,TY> y;
		ArcsMat<M,N,T>::template polycoeff<M+1,1,TY>(u, y, Tol);
		return y;
	}

}

// ArcsMatrixに付随するクラス
namespace ArcsMatrix {
	//! @brief MATファイル保存クラス(MATLAB Level 4対応)
	//! @tparam 
	//template <>
	class MatExport {
		public:
			//! @brief コンストラクタ
			//! @param[in]	FileName	MATファイル名 (拡張子matも含むこと)
			MatExport(const std::string& FileName)
				: fp(nullptr), MatFileName(FileName), MatHeader()
			{
				fp = std::fopen(FileName.c_str(), "wb");// バイナリ書き込みモードでMATファイルを新規作成
				arcs_assert(fp != nullptr);				// ファイル作成に失敗した場合
			}

			//! @brief ムーブコンストラクタ
			//! @param[in]	r	演算子右側
			MatExport(MatExport&& r)
				: fp(r.fp), MatFileName(r.MatFileName), MatHeader(r.MatHeader)
			{
				r.fp = nullptr;	// ムーブ元の所有権を解放
			}

			//! @brief ムーブ代入演算子
			//! @param[in]	r	演算子右側
			MatExport& operator=(MatExport&& r) noexcept {
				fp = r.fp;		// ムーブ元からムーブ先への所有権の移動
				r.fp = nullptr;	// ムーブ元の所有権を解放
				return *this;
			}

			//! @brief デストラクタ
			~MatExport(){
				std::fclose(fp);	// MATファイルを閉じる
			}

			//! @brief MATファイルへの行列データの書き出し
			//! @tparam	M, N, T	入力行列の高さ, 幅, 要素の型
			//! @param[in]	MatName		変数名
			//! @param[in]	U			入力行列
			template<size_t M, size_t N, typename T = double>
			void Save(const std::string& MatName, const ArcsMat<M,N,T> U){
				// ヘッダの初期化
				MatHeader.Type = 0000;		// データタイプの初期化
				MatHeader.NumOfRows = M;	// 行列の縦の長さ
				MatHeader.NumOfColumn = N;	// 行列の横の長さ
				MatHeader.LenOfName = MatName.length() + 1;	// 変数名の長さ＋１

				// 処理系のバイトオーダー(エンディアン)によってヘッダの設定を変える
				#if   __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
					//                MOPT
					MatHeader.Type += 0000;	// リトルエンディアンの場合(x86系)
				#elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
					//                MOPT
					MatHeader.Type += 1000;	// ビッグエンディアンの場合
				#else
					static_assert(false);	// エンディアンが不明の場合
				#endif

				// データ型によってヘッダの設定を変える
				if constexpr(std::is_same_v<double, T> || std::is_same_v<std::complex<double>, T>){
					//                MOPT
					MatHeader.Type += 0000;
				}else if constexpr(std::is_same_v<float, T> || std::is_same_v<std::complex<float>, T>){
					//                MOPT
					MatHeader.Type += 0010;
				}else if constexpr(std::is_same_v<int32_t, T>){
					//                MOPT
					MatHeader.Type += 0020;
				}else if constexpr(std::is_same_v<int16_t, T>){
					//                MOPT
					MatHeader.Type += 0030;
				}else if constexpr(std::is_same_v<uint16_t, T>){
					//                MOPT
					MatHeader.Type += 0040;
				}else if constexpr(std::is_same_v<uint8_t, T>){
					//                MOPT
					MatHeader.Type += 0050;
				}else{
					arcs_assert(false);	// ここには来ない
				}

				// 複素数型か実数型でヘッダの設定を変える
				if constexpr(ArcsMatrix::IsComplex<T>){
					MatHeader.HasImag = 1;
				}else{
					MatHeader.HasImag = 0;
				}

				// MATファイルにヘッダ部分を書き出す
				std::fwrite(&MatHeader, sizeof(MatHeader), 1, fp);						// ヘッダ書き出し
				std::fwrite(MatName.c_str(), sizeof(char), MatHeader.LenOfName, fp);	// 変数名書き出し

				// 複素数型か実数型で書き出し動作を変える
				if constexpr(ArcsMatrix::IsComplex<T>){
					// 複素数型の場合
					const ArcsMat<M*N,1,T> vecU  = ArcsMat<M,N,T>::vec(U);	// 縦ベクトル化
					const auto RealU = ArcsMat<M*N,1,T>::real(vecU);		// 実数部を抽出
					const auto ImagU = ArcsMat<M*N,1,T>::imag(vecU);		// 虚数部を抽出
					const std::array<std::array<double, M*N>, 1> MatReal = RealU.GetData();	// 実数部を1次元配列へ格納
					const std::array<std::array<double, M*N>, 1> MatImag = ImagU.GetData();	// 虚数部を1次元配列へ格納
					std::fwrite(MatReal.at(0).data(), sizeof(double), M*N, fp);				// 行列の実数部要素書き出し
					std::fwrite(MatImag.at(0).data(), sizeof(double), M*N, fp);				// 行列の虚数部要素書き出し
				}else{
					// 実数型の場合
					const std::array<std::array<T, M*N>, 1> MatData = (ArcsMat<M,N,T>::vec(U)).GetData();	// 縦ベクトル化して1次元配列へ格納
					std::fwrite(MatData.at(0).data(), sizeof(T), M*N, fp);									// 行列要素書き出し
				}
			}

		private:
			MatExport(const MatExport&) = delete;					//!< コピーコンストラクタ使用禁止
			const MatExport& operator=(const MatExport&) = delete;	//!< コピー代入演算子使用禁止

			FILE* fp;					//!< ファイルポインタ
			std::string MatFileName;	//!< MATファイル名
			
			//!@brief ヘッダデータの定義
			struct {
				// M = 0 リトルエンディアン, M = 1 ビッグエンディアン
				// O = 0 予約、常にゼロ
				// P = 0 double, P = 1 single, P = 2 int32_t, P = 3 int16_t, P = 4 uint16_t, P = 5 uint8_t 
				// T = 0 Numeric, T = 1 Text, T = 2 Sparse
				//                     MOPT
				uint32_t Type		 = 0000;	// データタイプの設定
				uint32_t NumOfRows	 = 1;		// 行数M
				uint32_t NumOfColumn = 1;		// 列数N
				uint32_t HasImag	 = 0;		// = 0 実数データ, = 1 複素数データ
				uint32_t LenOfName	 = 0;		// 変数名の長さ＋１
			} MatHeader;

	};
}

//! @brief ARCS-Variable-Matrix 可変行列演算クラス
//! @tparam T	データ型(デフォルトはdouble型)
template <typename T = double>
class ArcsVarMat {
	private:
		// 非公開版基本定数
		static constexpr size_t ITERATION_MAX = 10000;	//!< 反復計算の最大値
		static constexpr char CMPLX_UNIT = 'i';			//!< 虚数単位記号

		// 内部処理用
		size_t M;		//!< 高さ(行数)
		size_t N;		//!< 縦幅(列数)
		size_t Mindex;	//!< 縦方向カウンタ
		size_t Nindex;	//!< 横方向カウンタ

		// 行列データの実体
		// ArcsVarMatは行列の縦方向にメモリアドレスが連続しているので、縦ベクトル優先。
		std::vector<std::vector<T>> Data;	//!< データ格納用変数 配列要素の順番は Data[N列(横方向)][M行(縦方向)]
		
	public:
		//! @brief コンストラクタ
		//! @param	m	行列の高さ
		//! @param	n	行列の幅
		constexpr ArcsVarMat(const size_t m = 1, const size_t n = 1) noexcept
			: M(m), N(n), Mindex(0), Nindex(0), Data({0})
		{
			static_assert(ArcsMatrix::IsApplicable<T>, "ArcsMat: Type Error");	// 対応可能型チェック
			arcs_assert(m != 0);	// サイズゼロの行列は禁止
			arcs_assert(n != 0);	// サイズゼロの行列は禁止

			Resize(M, N);	// サイズを指定
			FillAll(0);		// すべての要素を零で初期化
		}

		//! @brief コンストラクタ(初期化リスト版)
		//! @param[in]	InitList	初期化リスト
		constexpr ArcsVarMat(const size_t m, const size_t n, const std::initializer_list<T> InitList) noexcept
			: M(m), N(n), Mindex(0), Nindex(0), Data({0})
		{
			static_assert(ArcsMatrix::IsApplicable<T>, "ArcsMat: Type Error");	// 対応可能型チェック
			arcs_assert(N != 0);	// サイズゼロの行列は禁止
			arcs_assert(M != 0);	// サイズゼロの行列は禁止

			Resize(M, N);			// サイズを指定

			const T* ListVal = InitList.begin();		// 初期化リストの最初のポインタ位置
			size_t Ni = 0;				// 横方向カウンタ
			size_t Mi = 0;				// 縦方向カウンタ
			for(size_t i = 0; i < InitList.size(); ++i){
				// 初期化リストを順番に読み込んでいく
				arcs_assert(Ni < N);	// 横方向カウンタが行の長さ以内かチェック
				arcs_assert(Mi < M);	// 縦方向カウンタが列の高さ以内かチェック
				Data[Ni][Mi] = static_cast<T>(ListVal[i]);	// キャストしてから行列の要素を埋める
				Ni++;					// 横方向カウンタをカウントアップ
				if(Ni == N){			// 横方向カウンタが最後まで行き着いたら，
					Ni = 0;				// 横方向カウンタを零に戻して，
					Mi++;				// その代わりに，縦方向カウンタをカウントアップ
				}
			}
		}

		//! @brief コピーコンストラクタ(型が同じ行列の場合の定義)
		//! @param[in]	right	演算子右側
		constexpr ArcsVarMat(const ArcsVarMat<T>& right) noexcept
			: M(1), N(1), Mindex(0), Nindex(0), Data({0})
		{
			static_assert(ArcsMatrix::IsApplicable<T>, "ArcsMat: Type Error");	// 対応可能型チェック

			// メンバを取り込む処理
			Resize(right.M, right.N);
			Data = right.Data;
		}

		//! @brief コピーコンストラクタ(型が違う行列の場合の定義)
		//! @tparam	R	演算子右側の行列の要素の型
		//! @param[in]	right	演算子右側
		template<typename R = double>
		constexpr ArcsVarMat(const ArcsVarMat<R>& right) noexcept
			: M(1), N(1), Mindex(0), Nindex(0), Data({0})
		{
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			arcs_assert(M == right.M);	// 行列のサイズチェック
			arcs_assert(N == right.N);	// 行列のサイズチェック

			// メンバを取り込む処理
			Resize(right.M, right.N);
			
			// 型の種類によってキャスト処理を変更
			if constexpr(ArcsMatrix::IsIntFloat<R> && ArcsMatrix::IsIntFloat<T>){
				// 「整数・浮動小数点型 → 整数・浮動小数点型」のキャスト処理
				Data = static_cast<T>(right.GetData());
			}else if constexpr(ArcsMatrix::IsIntFloat<R> && ArcsMatrix::IsComplex<T>){
				// 「整数・浮動小数点型 → 複素数型」のキャスト処理
				arcs_assert(false);	// 未実装
				const std::vector<std::vector<R>>& RightData = right.ReadOnlyRef();
				for(size_t i = 0; i < N; ++i){
					//for(size_t j = 0; j < M; ++j) Data[i][j].real(RightData[i][j]);	// 未実装
				}
			}else if constexpr(ArcsMatrix::IsComplex<R> && ArcsMatrix::IsComplex<T>){
				// 「複素数型 → 複素数型」のキャスト処理
				Data = static_cast<T>(right.GetData());
			}else{
				// キャスト不能の場合の処理
				arcs_assert(false);	// 変換不能、もしここに来たら問答無用でAssertion Failed
			}
		}
		
		//! @brief ムーブコンストラクタ
		//! @param[in]	right	演算子右側
		constexpr ArcsVarMat(ArcsVarMat<T>&& right) noexcept
			: M(right.M), N(right.N), Mindex(0), Nindex(0), Data(right.Data)
		{
			// メンバを取り込む以外の処理は無し
		}
		
		//! @brief ムーブコンストラクタ(型が違う行列の場合の定義)
		//! @tparam	R	演算子右側の行列要素の型
		//! @param[in]	right	演算子右側
		template<typename R = double>
		constexpr ArcsVarMat(ArcsVarMat<R>&& right) noexcept
			: M(right.M), N(right.N), Mindex(0), Nindex(0), Data(right.Data)
		{
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			arcs_assert(M == right.M);	// 行列のサイズチェック
			arcs_assert(N == right.N);	// 行列のサイズチェック

			// メンバを取り込む処理
			Resize(right.M, right.N);
			
			// 型の種類によってキャスト処理を変更
			if constexpr(ArcsMatrix::IsIntFloat<R> && ArcsMatrix::IsIntFloat<T>){
				// 「整数・浮動小数点型 → 整数・浮動小数点型」のキャスト処理
				Data = static_cast<T>(right.GetData());
			}else if constexpr(ArcsMatrix::IsIntFloat<R> && ArcsMatrix::IsComplex<T>){
				// 「整数・浮動小数点型 → 複素数型」のキャスト処理
				arcs_assert(false);	// 未実装
				const std::vector<std::vector<R>>& RightData = right.ReadOnlyRef();
				for(size_t i = 0; i < N; ++i){
					//for(size_t j = 0; j < M; ++j) Data[i][j].real(RightData[i][j]);	// 未実装
				}
			}else if constexpr(ArcsMatrix::IsComplex<R> && ArcsMatrix::IsComplex<T>){
				// 「複素数型 → 複素数型」のキャスト処理
				Data = static_cast<T>(right.GetData());
			}else{
				// キャスト不能の場合の処理
				arcs_assert(false);	// 変換不能、もしここに来たら問答無用でAssertion Failed
			}
		}

		//! @brief 行列コピー代入演算子(サイズと型が同じ同士の行列の場合)
		//! @param[in] right 演算子の右側
		//! @return 結果
		constexpr ArcsVarMat<T>& operator=(const ArcsVarMat<T>& right) noexcept {
			static_assert(ArcsMatrix::IsApplicable<T>, "ArcsMat: Type Error");	// 対応可能型チェック
			arcs_assert(M == right.M);	// 行列のサイズチェック
			arcs_assert(N == right.N);	// 行列のサイズチェック
			//for(size_t i = 1; i <= N; ++i){
			//	for(size_t j = 1; j <= M; ++j) (*this)(j,i) = right(j,i);
				// 上記は次とほぼ同等→ this->Data[i][j] = right.Data[i][j];
			//}
			Data = right.Data;
			return (*this);
		}
		
		//! @brief 行列コピー代入演算子(型が違う行列の場合の定義)
		//! @tparam	R	演算子右側の行列要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<typename R = double>
		constexpr ArcsVarMat<T>& operator=(const ArcsVarMat<R>& right) noexcept {
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			arcs_assert(M == right.GetHeight());	// 行列のサイズチェック
			arcs_assert(N == right.GetWidth());		// 行列のサイズチェック

			// 型の種類によってキャスト処理を変更
			if constexpr(ArcsMatrix::IsIntFloat<R> && ArcsMatrix::IsIntFloat<T>){
				// 「整数・浮動小数点型 → 整数・浮動小数点型」のキャスト処理
				const std::vector<std::vector<R>>& RightData = right.ReadOnlyRef();
				for(size_t i = 0; i < N; ++i){
					for(size_t j = 0; j < M; ++j) Data[i][j] = static_cast<T>(RightData[i][j]);
				}
			}else if constexpr(ArcsMatrix::IsIntFloat<R> && ArcsMatrix::IsComplex<T>){
				// 「整数・浮動小数点型 → 複素数型」のキャスト処理
				arcs_assert(false);	// 未実装
				const std::vector<std::vector<R>>& RightData = right.ReadOnlyRef();
				for(size_t i = 0; i < N; ++i){
					//for(size_t j = 0; j < M; ++j) Data[i][j].real(RightData[i][j]);	// 未実装
				}
			}else if constexpr(ArcsMatrix::IsComplex<R> && ArcsMatrix::IsComplex<T>){
				// 「複素数型 → 複素数型」のキャスト処理
				Data = static_cast<T>(right.GetData());
			}else{
				// キャスト不能の場合の処理
				arcs_assert(false);	// 変換不能、もしここに来たら問答無用でAssertion Failed
			}
			
			return (*this);
		}

		//! @brief 縦ベクトル添字演算子(縦ベクトルのm番目の要素の値を返す。x = A(m,1)と同じ意味)
		//!        備考：ArcsVarMatは縦ベクトル優先なので、横ベクトル添字演算子は無い。
		//! @param[in]	m	縦方向の要素番号( "1" 始まり)
		//! @return	要素の値
		constexpr T operator[](const size_t m) const{
			arcs_assert(N == 1);		// 縦ベクトルチェック
			return Data[0][m - 1];
		}
		
		//! @brief 縦ベクトル添字演算子(縦ベクトルのm番目の要素に値を設定する。A(m,1) = xと同じ意味)
		//!        備考：ArcsVarMatは縦ベクトル優先なので、横ベクトル添字演算子は無い。
		//! @param[in]	m	縦方向の要素番号( "1" 始まり)
		//! @return	設定後の縦ベクトル
		constexpr T& operator[](const size_t m){
			arcs_assert(N == 1);		// 縦ベクトルチェック
			return Data[0][m - 1];
		}

		//! @brief 行列括弧演算子(行列の(m,n)要素の値を返す。サイズチェック無し版)
		//! @param[in]	m	m行目(縦方向の位置)
		//! @param[in]	n	n列目(横方向の位置)
		//! @return	要素の値
		constexpr T operator()(const size_t m, const size_t n) const{
			return Data[n - 1][m - 1];
		}
		
		//! @brief 行列括弧演算子(行列の(m,n)要素に値を設定する。サイズチェック無し版)
		//! @param[in]	m	m行目(縦方向の位置)
		//! @param[in]	n	n列目(横方向の位置)
		//! @return	設定後の行列
		constexpr T& operator()(const size_t m, const size_t n){
			return Data[n - 1][m - 1];
		}
		
		//! @brief 行列括弧演算子(行列の(m,n)要素の値を返す。サイズチェック可能版)
		//! @param[in]	m	m行目(縦方向の位置)
		//! @param[in]	n	n列目(横方向の位置)
		//! @param[in]	chk	サイズチェックフラグ true = サイズチェックする, false = サイズチェックしない
		//! @return	要素の値
		constexpr T operator()(const size_t m, const size_t n, const bool chk) const{
			if(chk == true){
				arcs_assert(0 < m && m <= M);
				arcs_assert(0 < n && n <= N);
			}
			return Data[n - 1][m - 1];
		}
		
		//! @brief 行列括弧演算子(行列の(m,n)要素に値を設定する。サイズチェック可能版)
		//! @param[in]	m	m行目(縦方向の位置)
		//! @param[in]	n	n列目(横方向の位置)
		//! @param[in]	chk	サイズチェックフラグ true = サイズチェックする, false = サイズチェックしない
		//! @return	設定後の行列
		constexpr T& operator()(const size_t m, const size_t n, const bool chk){
			if(chk == true){
				arcs_assert(0 < m && m <= M);
				arcs_assert(0 < n && n <= N);
			}
			return Data[n - 1][m - 1];
		}

		//! @brief 単項プラス演算子
		//! @return 結果
		constexpr ArcsVarMat<T> operator+(void) const{
			ArcsVarMat<T> ret;
			//for(size_t i = 1; i <= N; ++i){
			//	for(size_t j = 1; j <= M; ++j) ret(j,i) = (*this)(j,i);
				// 上記は次とほぼ同等→ ret.Data[i][j] = Data[i][j];
			//}
			ret.Data = Data;
			return ret;
		}
		
		//! @brief 単項マイナス演算子
		//! @return 結果
		constexpr ArcsVarMat<T> operator-(void) const{
			ArcsVarMat<T> ret;
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) ret(j,i) = -(*this)(j,i);
				// 上記は次とほぼ同等→ ret.Data[i][j] = -Data[i][j];
			}
			return ret;
		}

		//! @brief 行列サイズを変更する関数
		//! @param	m	行列の高さ
		//! @param	n	行列の幅
		constexpr void Resize(const size_t m, const size_t n){
			M = m;
			N = n;
			Data.resize(n);			// 幅をN列に変更
			for(auto& column : Data){
				column.resize(m);	// 高さをM行に変更
			}
		}

		//! @brief 行列要素の各メモリアドレスを表示する関数
		constexpr void DispAddress(void) const{
			if(__builtin_constant_p(Data) == true) return;	// コンパイル時には処理を行わない

			printf("Mem addr of matrix:\n");
			for(size_t j = 0; j < M; ++j){
				printf("[ ");
				for(size_t i = 0; i < N; ++i){
					printf("%p ", &(Data[i][j]));
				}
				printf("]\n");
			}
			printf("\n");
		}

		//! @brief 行列の要素を表示
		//! @param[in] format	表示形式 (%1.3e とか %5.3f とか printfと同じ)
		constexpr void Disp(const std::string& format = "% g") const{
			if(__builtin_constant_p(Data) == true) return;	// コンパイル時には処理を行わない

			for(size_t j = 0; j < M; ++j){
				printf("[ ");
				for(size_t i = 0; i < N; ++i){
					// データ型によって表示方法を変える
					if constexpr(ArcsMatrix::IsComplex<T>){
						// 複素数型の場合
						// 実数部の表示
						printf(format.c_str(), Data[i][j].real());
						// 虚数部の表示
						if(0.0 <= Data[i][j].imag()){
							printf(" +");
						}else{
							printf(" -");
						}
						printf( format.c_str(), std::abs(Data[i][j].imag()) );
						printf("%c", CMPLX_UNIT);
					}else{
						// それ以外の場合
						printf(format.c_str(), Data[i][j]);
					}
					printf(" ");
				}
				printf("]\n");
			}
			printf("\n");
		}
		
		//! @brief 行列のサイズを表示
		constexpr void DispSize(void) const{
			if(__builtin_constant_p(Data) == true) return;	// コンパイル時には処理を行わない

			printf("[ Height: %zu  x  Width: %zu ]\n\n", M, N);
		}
		
		//! @brief 行列の高さ(行数)を返す関数
		//! @return 行列の幅
		constexpr size_t GetHeight(void) const{
			return M;
		}
		
		//! @brief 行列の幅(列数)を返す関数
		//! @return 行列の幅
		constexpr size_t GetWidth(void) const{
			return N;
		}
		
		//! @brief 行列要素に値を設定する関数
		//! @tparam	T1, T2	要素の型
		//! @param[in]	u1	要素1の値
		//! @param[in]	u2	要素2以降の値
		template<typename T1, typename... T2>				// 可変長引数テンプレート
		constexpr void Set(const T1& u1, const T2&... u2){	// 再帰で順番に可変長引数を読み込んでいく
			static_assert(ArcsMatrix::IsApplicable<T1>, "ArcsMat: Type Error");	// 対応可能型チェック
			arcs_assert(Mindex < M);		// 縦方向カウンタが高さ以内かチェック
			arcs_assert(Nindex < N);		// 横方向カウンタが幅以内かチェック
			Data[Nindex][Mindex] = static_cast<T>(u1);	// キャストしてから行列の要素を埋める
			Nindex++;			// 横方向カウンタをインクリメント
			if(Nindex == N){	// 横方向カウンタが最後まで行き着いたら，
				Nindex = 0;		// 横方向カウンタを零に戻して，
				Mindex++;		// その代わりに，縦方向カウンタをインクリメント
			}
			Set(u2...);			// 自分自身を呼び出す(再帰)
		}
		constexpr void Set(){
			// 再帰の最後に呼ばれる関数
			Nindex = 0;	// すべての作業が終わったので，
			Mindex = 0;	// 横方向と縦方向カウンタを零に戻しておく
		}
		
		//! @brief 行列要素から値を読み込む関数
		//! @tparam	T1, T2	要素の型
		//! @param[in]	u1	要素1の値
		//! @param[in]	u2	要素2以降の値
		template<typename T1, typename... T2>	// 可変長引数テンプレート
		constexpr void Get(T1& u1, T2&... u2){	// 再帰で順番に可変長引数を読み込んでいく
			static_assert(ArcsMatrix::IsApplicable<T1>, "ArcsMat: Type Error");	// 対応可能型チェック
			arcs_assert(Mindex < M);		// 縦方向カウンタが高さ以内かチェック
			arcs_assert(Nindex < N);		// 横方向カウンタが幅以内かチェック
			u1 = static_cast<T1>(Data[Nindex][Mindex]);	// 行列の要素からキャストして読み込み
			Nindex++;			// 横方向カウンタをインクリメント
			if(Nindex == N){	// 横方向カウンタが最後まで行き着いたら，
				Nindex = 0;		// 横方向カウンタを零に戻して，
				Mindex++;		// その代わりに，縦方向カウンタをインクリメント
			}
			Get(u2...);			// 自分自身を呼び出す(再帰)
		}
		constexpr void Get(){
			// 再帰の最後に呼ばれる関数
			Nindex = 0;	// すべての作業が終わったので，
			Mindex = 0;	// 横方向と縦方向カウンタを零に戻しておく
		}
		
		//! @brief すべての要素を指定した値で埋める関数
		//! @tparam	R	要素の型
		//! @param[in] u 埋める値
		template<typename R>
		constexpr void FillAll(const R& u){
			static_assert(ArcsMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j){
					// 型の種類によって値埋め処理を変更
					if constexpr(ArcsMatrix::IsIntFloat<R> && ArcsMatrix::IsIntFloat<T>){
						// 「整数・浮動小数点型 → 整数・浮動小数点型」の値埋め処理
						Data[i][j] = static_cast<T>(u);
					}else if constexpr(ArcsMatrix::IsIntFloat<R> && ArcsMatrix::IsComplex<T>){
						// 「整数・浮動小数点型 → 複素数型」の値埋め処理
						Data[i][j].real(u);
						Data[i][j].imag(0);
					}else if constexpr(ArcsMatrix::IsComplex<R> && ArcsMatrix::IsComplex<T>){
						// 「複素数型 → 複素数型」の値埋め処理
						Data[i][j] = static_cast<T>(u);
					}else{
						// 値埋め不能の場合の処理
						arcs_assert(false);	// 変換不能、もしここに来たら問答無用でAssertion Failed
					}
				}
			}
		}

		//! @brief std:arrayの2次元配列データをそのまま返す関数
		//! @return	2次元vectorコンテナ
		constexpr std::vector<std::vector<T>> GetData(void) const{
			return Data;
		}

		//! @brief std:arrayの2次元配列データの読み込み専用の参照を返す関数
		//! @return	2次元vectorコンテナの参照
		constexpr const std::vector<std::vector<T>>& ReadOnlyRef(void) const{
			return Data;
		}

};


}

#endif

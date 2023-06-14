//! @file ArcsMatrix.hh
//! @brief ARCS-Matrix 行列演算クラス
//!
//! 行列に関係する様々な演算を実行するクラス
//!
//! @date 2022/08/21
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2022 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the BSD License.
// For details, see the License.txt file.
//
// 以下，コメント。
// ・各々の関数における計算結果はMATLAB/Maximaと比較して合っていることを確認している。
// ・ただし，LU分解やコレスキー分解などの一見した表現が定まらない関数では，当然，MATLABとは異なる結果を出力するように見える。
// ・動的メモリ版に比べてかなり高速の行列演算が可能。

#ifndef ARCSMATRIX
#define ARCSMATRIX

#include <string>
#include <tuple>
#include <cmath>
#include <cassert>
#include <array>
#include <complex>

// ARCS組込み用マクロ
#ifdef ARCS_IN
	// ARCSに組み込まれる場合
	#include "ARCSassert.hh"
#else
	// ARCSに組み込まれない場合
	#define arcs_assert(a) (assert(a))
#endif

// 表示用マクロ
#define dispsize(a)  (dispsize_macro((a),#a))		//!< 行列サイズ表示マクロ
#define dispmatfmt(a,b) (dispmatfmt_macro((a),b,#a))//!< 行列要素表示マクロ(フォーマット指定あり版)
#define dispmat(a) (dispmat_macro((a),#a))			//!< 行列要素表示マクロ(フォーマット指定なし版)

namespace ARCS {	// ARCS名前空間

//! @brief ノルム計算方法の定義
enum class NormType {
	AMT_INFINITY,	//!< 無限大ノルム
	AMT_EUCLID		//!< ユークリッドノルム
};

//! @brief 行列の状態の定義
enum class MatStatus {
	AMT_NA,		//!< 状態定義該当なし
	AMT_LU_ODD,	//!< LU分解したときに並べ替えが奇数回発生
	AMT_LU_EVEN	//!< LU分解したときに並べ替えが偶数回発生
};

//! @brief ARCS-Matrix 行列演算クラス
//! @tparam M	行列の高さ
//! @tparam	N	行列の幅
//! @tparam T	データ型(デフォルトはdouble型)
template <size_t M, size_t N, typename T = double>
class ArcsMat {
	public:
		//! @brief コンストラクタ
		constexpr ArcsMat(void)
			: Nindex(0), Mindex(0), Data({0})
		{
			static_assert(N != 0);	// サイズゼロの行列は禁止
			static_assert(M != 0);	// サイズゼロの行列は禁止
			FillAll(0);				// すべての要素を零で初期化
		}
		
		//! @brief コンストラクタ(任意初期値版)
		//! @tparam	R	演算子右側の要素の型
		//! @param[in]	InitValue	行列要素の初期値
		template<typename R>
		constexpr explicit ArcsMat(const R InitValue)
			: Nindex(0), Mindex(0), Data({0})
		{
			static_assert(N != 0);	// サイズゼロの行列は禁止
			static_assert(M != 0);	// サイズゼロの行列は禁止
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			FillAll(InitValue);		// すべての要素を指定した値で初期化
		}
		
		//! @brief コンストラクタ(初期化リスト版)
		//! @tparam	R	演算子右側の要素の型
		//! @param[in]	InitList	初期化リスト
		template<typename R>
		constexpr ArcsMat(const std::initializer_list<R> InitList)
			: Nindex(0), Mindex(0), Data({0})
		{
			static_assert(N != 0);	// サイズゼロの行列は禁止
			static_assert(M != 0);	// サイズゼロの行列は禁止
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			const R* ListVal = InitList.begin();		// 初期化リストの最初のポインタ位置
			size_t Ni = 0;				// 横方向カウンタ
			size_t Mi = 0;				// 縦方向カウンタ
			for(size_t i = 0; i < InitList.size(); ++i){
				// 初期化リストを順番に読み込んでいく
				arcs_assert(Ni < N);	// 横方向カウンタが行の長さ以内かチェック
				arcs_assert(Mi < M);	// 縦方向カウンタが列の高さ以内かチェック
				Data[Ni][Mi] = (T)ListVal[i];	// キャストしてから行列の要素を埋める
				Ni++;					// 横方向カウンタをカウントアップ
				if(Ni == N){			// 横方向カウンタが最後まで行き着いたら，
					Ni = 0;				// 横方向カウンタを零に戻して，
					Mi++;				// その代わりに，縦方向カウンタをカウントアップ
				}
			}
		}
		
		//! @brief コピーコンストラクタ
		//! @param[in]	right	右辺値
		constexpr ArcsMat(const ArcsMat<M,N,T>& right)
			: Nindex(0), Mindex(0), Data(right.GetData())
		{
			
		}
		
		//! @brief コピーコンストラクタ(サイズと型が違う行列の場合, エラー検出用の定義)
		//! @tparam	P, Q, R	演算子右側の行列の高さ, 幅, 要素の型
		//! @param[in]	right	右辺値
		template<size_t P, size_t Q, typename R = double>
		constexpr ArcsMat(const ArcsMat<P,Q,R>& right)
			: Nindex(0), Mindex(0), Data(right.GetData())
		{
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			arcs_assert(false);	// もしここに来たら問答無用でAssertion Failed
		}
		
	
		//! @brief ムーブコンストラクタ
		//! @param[in]	right	右辺値
		constexpr ArcsMat(ArcsMat<M,N,T>&& right)
			: Nindex(0), Mindex(0), Data(right.GetData())
		{
			
		}
		
		//! @brief ムーブコンストラクタ(サイズと型が違う行列の場合, エラー検出用の定義)
		//! @tparam	P, Q, R	演算子右側の行列の高さ, 幅, 要素の型
		//! @param[in]	right	右辺値
		template<size_t P, size_t Q, typename R = double>
		constexpr ArcsMat(ArcsMat<M,N,T>&& right)
			: Nindex(0), Mindex(0), Data(right.GetData())
		{
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			arcs_assert(false);	// もしここに来たら問答無用でAssertion Failed
		}
		
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
		
		//! @brief 行列代入演算子(サイズと型が同じ同士の行列の場合)
		//! @param[in] right 演算子の右側
		//! @return 結果
		constexpr ArcsMat<M,N,T>& operator=(const ArcsMat<M,N,T>& right){
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) this->Data[i][j] = right.Data[i][j];
			}
			return (*this);
		}
		
		//! @brief 行列代入演算子(サイズと型が違う行列の場合, エラー検出用の定義)
		//! @tparam	P, Q, R	演算子右側の行列の高さ, 幅, 要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<size_t P, size_t Q, typename R = double>
		constexpr ArcsMat<M,N,T>& operator=(const ArcsMat<P,Q,R>& right){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			arcs_assert(false);	// もしここに来たら問答無用でAssertion Failed
			return (*this);
		}
		
		//! @brief 単項プラス演算子
		//! @return 結果
		constexpr ArcsMat<M,N,T> operator+(void) const{
			ArcsMat<M,N,T> ret;
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) ret.Data[i][j] = Data[i][j];
			}
			return ret;
		}
		
		//! @brief 単項マイナス演算子
		//! @return 結果
		constexpr ArcsMat<M,N,T> operator-(void) const{
			ArcsMat<M,N,T> ret;
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) ret.Data[i][j] = -Data[i][j];
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
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			ArcsMat<M,N,T> ret;
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) ret.Data[i][j] = Data[i][j] + right.Data[i][j];
			}
			return ret;
		}
		
		//! @brief 行列加算演算子(行列＝行列＋スカラーの場合)
		//! @tparam	R	演算子右側の要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<typename R>
		constexpr ArcsMat<M,N,T> operator+(const R& right) const{
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			ArcsMat<M,N,T> ret;
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) ret.Data[i][j] = Data[i][j] + right;
			}
			return ret;
		}
		
		//! @brief 行列減算演算子(行列＝行列ー行列の場合)
		//! @tparam	P, Q, R	演算子右側の行列の高さ, 幅, 要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<size_t P, size_t Q, typename R = double>
		constexpr ArcsMat<M,N,T> operator-(const ArcsMat<P,Q,R>& right) const{
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			ArcsMat<M,N,T> ret;
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) ret.Data[i][j] = Data[i][j] - right.Data[i][j];
			}
			return ret;
		}
		
		//! @brief 行列減算演算子(行列＝行列－スカラーの場合)
		//! @tparam	R	演算子右側の要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<typename R>
		constexpr ArcsMat<M,N,T> operator-(const R& right) const{
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			ArcsMat<M,N,T> ret;
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) ret.Data[i][j] = Data[i][j] - right;
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
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			ArcsMat<M,Q,T> ret;
			for(size_t k = 0; k < Q; ++k){
				for(size_t i = 0; i < N; ++i){
					for(size_t j = 0; j < M; ++j) ret.Data[k][j] += Data[i][j]*right.Data[k][i];
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
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			ArcsMat<M,N,T> ret;
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) ret.Data[i][j] = Data[i][j]*right;
			}
			return ret;
		}
		
		//! @brief 行列スカラー除算演算子(行列＝行列／スカラーの場合)
		//! @tparam	R	演算子右側の要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<typename R>
		constexpr ArcsMat<M,N,T> operator/(const R& right) const{
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			ArcsMat<M,N,T> ret;
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) ret.Data[i][j] = Data[i][j]/right;
			}
			return ret;
		}
		
		//! @brief 行列加算代入演算子(行列＝行列＋行列の場合)
		//! @tparam	P, Q, R	演算子右側の行列の高さ, 幅, 要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<size_t P, size_t Q, typename R = double>
		constexpr ArcsMat<M,N,T>& operator+=(const ArcsMat<P,Q,R>& right){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) this->Data[i][j] += right.Data[i][j];
			}
			return (*this);
		}
		
		//! @brief 行列加算代入演算子(行列＝行列＋スカラーの場合)
		//! @tparam	R	演算子右側の要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<typename R>
		constexpr ArcsMat<M,N,T>& operator+=(const R& right){
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) this->Data[i][j] += right;
			}
			return (*this);
		}
		
		//! @brief 行列減算代入演算子(行列＝行列－行列の場合)
		//! @tparam	P, Q, R	演算子右側の行列の高さ, 幅, 要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<size_t P, size_t Q, typename R = double>
		constexpr ArcsMat<M,N,T>& operator-=(const ArcsMat<P,Q,R>& right){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			for(size_t j = 0; j < M; ++j){
				for(size_t i = 0; i < N; ++i) this->Data[i][j] -= right.Data[i][j];
			}
			return (*this);
		}
		
		//! @brief 行列減算代入演算子(行列＝行列－スカラーの場合)
		//! @tparam	R	演算子右側の要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<typename R>
		constexpr ArcsMat<M,N,T>& operator-=(const R& right){
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) this->Data[i][j] -= right;
			}
			return (*this);
		}
		
		//! @brief 行列乗算代入演算子(行列＝行列＊行列の場合)
		//! @tparam	P, Q, R	演算子右側の行列の高さ, 幅, 要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<size_t P, size_t Q, typename R = double>
		constexpr ArcsMat<M,N,T>& operator*=(const ArcsMat<P,Q,R>& right){
			static_assert(M == N, "ArcsMat: Size Error");	// 正方行列チェック
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			(*this) = (*this)*right;
			return (*this);
		}
		
		//! @brief 行列乗算代入演算子(行列＝行列＊スカラーの場合)
		//! @tparam	R	演算子右側の要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<typename R>
		constexpr ArcsMat<M,N,T>& operator*=(const R& right){
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			(*this) = (*this)*right;
			return (*this);
		}
		
		//! @brief 行列べき乗演算子(正方行列のべき乗)
		//! @param[in] right 演算子の右側
		//! @return 結果
		constexpr ArcsMat<M,N,T> operator^(const size_t& right) const{
			static_assert(M == N, "ArcsMat: Size Error");					// 正方行列チェック
			ArcsMat<M,N,T> ret = ArcsMat<M,N,T>::eye();
			for(size_t k = 1; k <= right; ++k) ret *= (*this);
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
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			ArcsMat<M,N,T> ret;
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) ret.Data[i][j] = Data[i][j]*right.Data[i][j];
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
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			ArcsMat<M,N,T> ret;
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) ret.Data[i][j] = Data[i][j]/right.Data[i][j];
			}
			return ret;
		}
		
		//! @brief 行列加算演算子 (スカラー＋行列の場合)
		//! @param[in] left		左側のスカラー値
		//! @param[in] right	右側の行列
		constexpr friend ArcsMat<M,N,T> operator+(const T& left, const ArcsMat<M,N,T>& right){
			ArcsMat<M,N,T> ret;
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) ret.Data[i][j] = left + right.Data[i][j];
			}
			return ret;
		}
		
		//! @brief 行列減算演算子 (スカラー－行列の場合)
		//! @param[in] left		左側のスカラー値
		//! @param[in] right	右側の行列
		constexpr friend ArcsMat<M,N,T> operator-(const T& left, const ArcsMat<M,N,T>& right){
			ArcsMat<M,N,T> ret;
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) ret.Data[i][j] = left - right.Data[i][j];
			}
			return ret;
		}
		
		//! @brief 行列乗算演算子 (スカラー＊行列の場合)
		//! @param[in] left		左側のスカラー値
		//! @param[in] right	右側の行列
		constexpr friend ArcsMat<M,N,T> operator*(const T& left, const ArcsMat<M,N,T>& right){
			ArcsMat<N,M,T> ret;
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) ret.Data[i][j] = right.Data[i][j]*left;
			}
			return ret;
		}
		
		//! @brief 行列乗算演算子 (スカラー／行列の場合)
		//! @param[in] left		左側のスカラー値
		//! @param[in] right	右側の行列
		constexpr friend ArcsMat<M,N,T> operator/(const T& left, const ArcsMat<M,N,T>& right){
			arcs_assert(false);		// この演算子の使い方は使用禁止
			ArcsMat<N,M,T> ret;
			return ret;
		}
		
		//! @brief 行列要素の各メモリアドレスを表示する関数
		constexpr void DispAddress(void) const{
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
		//! @param[in] format	表示形式 (%1.3e とか %5.3f とか printfと同じ)
		constexpr void Disp(const std::string& format) const{
			for(size_t j = 0; j < M; ++j){
				printf("[ ");
				for(size_t i = 0; i < N; ++i){
					// データ型によって表示方法を変える
					if constexpr(std::is_same_v<T, std::complex<double>> || std::is_same_v<T, std::complex<float>>){
						// 複素数型の場合
						// 実数部の表示
						printf(format.c_str(), Data[i][j].real());
						// 虚数部の表示
						if(0.0 <= Data[i][j].imag()){
							printf(" + j");
						}else{
							printf(" - j");
						}
						printf( format.c_str(), std::abs(Data[i][j].imag()) );
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
			static_assert(std::is_convertible_v<T, T1>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			arcs_assert(Mindex < M);		// 縦方向カウンタが高さ以内かチェック
			arcs_assert(Nindex < N);		// 横方向カウンタが幅以内かチェック
			Data[Nindex][Mindex] = (T)u1;	// キャストしてから行列の要素を埋める
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
			static_assert(std::is_convertible_v<T, T1>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			arcs_assert(Mindex < M);		// 縦方向カウンタが高さ以内かチェック
			arcs_assert(Nindex < N);		// 横方向カウンタが幅以内かチェック
			u1 = (T1)Data[Nindex][Mindex];	// 行列の要素からキャストして読み込み
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
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) Data[i][j] = u;
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
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			for(size_t j = 0; j < M; ++j) Data[0][j] = Array[j];
		}
		
		//! @brief 縦ベクトルを1次元std::array配列に書き込む関数
		//! @tparam	P, R	配列の長さ, 要素の型
		//! @param[out]	Array	std::array配列(縦MM×横1)
		template<size_t P, typename R = double>
		constexpr void StoreArray(std::array<R, P>& Array) const{
			static_assert(N == 1, "ArcsMat: Vector Error");	// 縦ベクトルチェック
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			for(size_t j = 0; j < M; ++j) Array[j] = Data[0][j];
		}
		
		//! @brief 2次元std::array配列を行列として読み込む関数
		//! @tparam	P, Q, R	配列の高さ, 幅, 要素の型
		//! @param[in]	Array	std::array配列(縦MM×横NN)
		template<size_t P, size_t Q, typename R = double>
		constexpr void LoadArray(const std::array<std::array<R, P>, Q>& Array){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			Data = Array;
		}
		
		//! @brief 行列を2次元std::array配列に書き込む関数
		//! @tparam	P, Q, R	配列の高さ, 幅, 要素の型
		//! @param[in]	Array	std::array配列(縦MM×横NN)
		template<size_t P, size_t Q, typename R = double>
		constexpr void StoreArray(std::array<std::array<R, P>, Q>& Array) const{
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			Array = Data;
		}
		
		//! @brief std:arrayの2次元配列データをそのまま返す関数
		//! @return	2次元配列データ
		constexpr std::array<std::array<T, M>, N> GetData(void) const{
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
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
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
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
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
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
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
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			arcs_assert(0 < m && m <= M);			// サイズチェック
			arcs_assert(0 < n && Q + n - 1 <= N);	// はみ出しチェック
			for(size_t i = 1; i <= Q; ++i) (*this)(m, n + i - 1) = w(1,i);
		}
		
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
			if constexpr(std::is_same_v<T, std::complex<double>> || std::is_same_v<T, std::complex<float>>){
				// 複素数型の場合
				for(size_t i = 1; i <= N; ++i) ret(i,i) = std::complex(1.0, 0.0);	// 対角成分を 1 + j0 で埋める
			}else{
				// それ以外の場合
				for(size_t i = 1; i <= N; ++i) ret(i,i) = (T)1;	// 対角成分を 1 で埋める
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
		
		//! @brief 転置行列を返す関数 (引数渡し版)
		//! @tparam	P, Q, R	出力の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void tp(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(M == Q, "ArcsMat: Transpose Size Error");	// 転置サイズチェック
			static_assert(N == P, "ArcsMat: Transpose Size Error");	// 転置サイズチェック
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) Y(i,j) = U(j,i);
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
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
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
			static_assert(std::is_convertible_v<size_t, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			ArcsMat<M,N,T> Y;
			for(size_t i = 1; i <= N; ++i){
				setcolumn(Y,  getcolumn(U, (size_t)u(1,i)), i);
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
			static_assert(std::is_convertible_v<size_t, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			for(size_t i = 1; i <= N; ++i){
				swapcolumn(UY, i, uy(1,i));
				ArcsMat<P,Q,R>::swapcolumn(uy, i, uy(1,i));
			}
		}
		
		//! @brief 各列の総和を計算して横ベクトルを出力する関数 (引数渡し版)
		//! @tparam	P, Q, R	出力ベクトルの高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	y	出力ベクトル
		template<size_t P, size_t Q, typename R = size_t>
		static constexpr void sumcolumn(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& y){
			static_assert(P == 1, "ArcsMat: Vector Error");	// 横ベクトルチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			y.FillAllZero();
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) y(1,i) += U(j,i);
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
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			arcs_assert(0 < m && m <= M);		// 行が高さ以内かチェック
			arcs_assert(0 < n1 && n1 <= N);	// 列が幅以内かチェック
			arcs_assert(0 < n2 && n2 <= N);	// 列が幅以内かチェック
			arcs_assert(n1 <= n2);				// 開始列と終了列が入れ替わらないかチェック
			for(size_t i = n1; i <= n2; ++i) UY(m, i) = (T)a;
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
			static_assert(std::is_convertible_v<size_t, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			ArcsMat<M,N,T> Y;
			for(size_t j = 1; j <= M; ++j){
				setrow(Y,  getrow(U, (size_t)u(j,1)), j);
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
			static_assert(std::is_convertible_v<size_t, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			for(size_t j = 1; j <= M; ++j){
				swaprow(UY, j, uy(j,1));
				ArcsMat<P,Q,R>::swaprow(uy, j, uy(j,1));
			}
		}
		
		//! @brief 各行の総和を計算して縦ベクトルを出力する関数 (引数渡し版)
		//! @tparam	P, Q, R	出力ベクトルの高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	y	出力ベクトル
		template<size_t P, size_t Q, typename R = size_t>
		static constexpr void sumrow(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& y){
			static_assert(Q == 1, "ArcsMat: Vector Error");	// 縦ベクトルチェック
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			y.FillAllZero();
			for(size_t j = 1; j <= M; ++j){
				for(size_t i = 1; i <= N; ++i) y(j,1) += U(j,i);
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
		template<size_t P>
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
			static_assert(P <= M);		// 小行列の方が高さが小さいかチェック
			static_assert(Q <= N);		// 小行列の方が幅が小さいかチェック
			static_assert(std::is_convertible_v<size_t, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
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
			static_assert(P <= M);		// 小行列の方が高さが小さいかチェック
			static_assert(Q <= N);		// 小行列の方が幅が小さいかチェック
			static_assert(std::is_convertible_v<size_t, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
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
		
		//! @brief 行列の各要素を上にm行分シフトする関数(下段の行はゼロになる)(引数渡し版)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		//! @param[in]	m	シフトする行数(デフォルト値 = 1)
		template<size_t P, size_t Q, typename R = double>
		static constexpr void shiftup(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y, const size_t m = 1){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			
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
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			
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
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			
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
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			
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
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			static_assert(std::is_convertible_v<T, F>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
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
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			static_assert(std::is_convertible_v<T, F>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
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
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			static_assert(std::is_convertible_v<T, F>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			static_assert(std::is_convertible_v<T, L>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			static_assert(std::is_convertible_v<T, X>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
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
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			for(size_t j = 1; j <= M; ++j) Y(j,j) = u(j,1);
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
		//! @param[out]	y	出力ベクトル
		template<size_t P, size_t Q, typename R = double, size_t L = std::min(M,N)>
		static constexpr void getdiag(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& y){
			static_assert(Q == 1, "ArcsMat: Vector Error");	// 縦ベクトルチェック
			static_assert(P == L, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			for(size_t j = 1; j <= L; ++j) y(j,1) = U(j,j);
		}
				
		//! @brief 行列の対角要素を縦ベクトルとして取得する関数(戻り値渡し版)
		//! @tparam	L	対角要素の数
		//! @param[in]	U	入力行列
		//! @return	出力ベクトル
		template<size_t L = std::min(M,N)>
		static constexpr ArcsMat<L,1,T> getdiag(const ArcsMat<M,N,T>& U){
			ArcsMat<L,1,T> y;
			getdiag(U, y);
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
		
		//! @brief 行列要素の最大値の要素番号を返す関数(戻り値渡し版のみ)
		//! @param[in]	U	入力行列
		//! @return	結果(タプルで返す)
		static constexpr std::tuple<size_t,size_t> maxidx(const ArcsMat<M,N,T>& U){
			size_t k = 1, l = 1;
			for(size_t i = 1; i <= N; ++i){		// 横方向探索
				for(size_t j = 1; j <= M; ++j){	// 縦方向探索
					if( U(k,l) < U(j,i) ){
						k = j;
						l = i;
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
		
		//! @brief 行列要素の最小値の要素番号を返す関数(戻り値渡し版のみ)
		//! @param[in]	U	入力行列
		//! @return	結果(タプルで返す)
		static constexpr std::tuple<size_t,size_t> minidx(const ArcsMat<M,N,T>& U){
			size_t k = 1, l = 1;
			for(size_t i = 1; i <= N; ++i){		// 横方向探索
				for(size_t j = 1; j <= M; ++j){	// 縦方向探索
					if( U(k,l) > U(j,i) ){
						k = j;
						l = i;
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
		//! @tparam	P, Q, R	小行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void exp(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");		// 暗黙の型変換可能チェック
			static_assert(std::is_floating_point_v<T>, "ArcsMat: Type Error (Floating Point)");	// 型チェック
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) Y(j,i) = std::exp( U(j,i) );
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
		//! @tparam	P, Q, R	小行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void log(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");		// 暗黙の型変換可能チェック
			static_assert(std::is_floating_point_v<T>, "ArcsMat: Type Error (Floating Point)");	// 型チェック
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) Y(j,i) = std::log( U(j,i) );
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
		//! @tparam	P, Q, R	小行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void log10(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");		// 暗黙の型変換可能チェック
			static_assert(std::is_floating_point_v<T>, "ArcsMat: Type Error (Floating Point)");	// 型チェック
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) Y(j,i) = std::log10( U(j,i) );
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
		//! @tparam	P, Q, R	小行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void sin(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");		// 暗黙の型変換可能チェック
			static_assert(std::is_floating_point_v<T>, "ArcsMat: Type Error (Floating Point)");	// 型チェック
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) Y(j,i) = std::sin( U(j,i) );
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
		//! @tparam	P, Q, R	小行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void cos(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");		// 暗黙の型変換可能チェック
			static_assert(std::is_floating_point_v<T>, "ArcsMat: Type Error (Floating Point)");	// 型チェック
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) Y(j,i) = std::cos( U(j,i) );
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
		//! @tparam	P, Q, R	小行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void tan(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");		// 暗黙の型変換可能チェック
			static_assert(std::is_floating_point_v<T>, "ArcsMat: Type Error (Floating Point)");	// 型チェック
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) Y(j,i) = std::tan( U(j,i) );
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

		//! @brief 行列要素の平方根を計算する関数(引数渡し版)
		//! @tparam	P, Q, R	小行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		template<size_t P, size_t Q, typename R = double>
		static constexpr void sqrt(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");		// 暗黙の型変換可能チェック
			static_assert(std::is_floating_point_v<T>, "ArcsMat: Type Error (Floating Point)");	// 型チェック
			for(size_t i = 1; i <= N; ++i){
				for(size_t j = 1; j <= M; ++j) Y(j,i) = std::sqrt( U(j,i) );
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

		/*
		//! @brief 行列要素の絶対値を返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend ArcsMat abse(const ArcsMat& U){
			ArcsMat Y;
			for(size_t i = 0; i < U.N; ++i){
				for(size_t j = 0; j < U.M; ++j) Y.Data[i][j] = std::abs(U.Data[i][j]);
			}
			return Y;
		}
		//! @brief 行列要素のtanhを返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend ArcsMat tanhe(const ArcsMat& U){
			ArcsMat Y;
			for(size_t i = 0; i < U.N; ++i){
				for(size_t j = 0; j < U.M; ++j) Y.Data[i][j] = std::tanh(U.Data[i][j]);
			}
			return Y;
		}
*/

		//! @brief n列目を左端として右上の上三角部分のみを返す関数(下三角部分はゼロ)(引数渡し版)
		//! @tparam	P, Q, R	出力行列の高さ, 幅, 要素の型
		//! @param[in]	U	入力行列
		//! @param[out]	Y	出力行列
		//! @param[in]	n	切り出す左端位置 n列目(デフォルト値 = 1)
		template<size_t P, size_t Q, typename R = double>
		static constexpr void gettriup(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& Y, const size_t n = 1){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック

			for(size_t j = 1; j <= M; ++j){
				for(size_t i = j + n - 1; i <= N; ++i){
					Y(j,i) = U(j,i);
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
			static_assert(std::is_convertible_v<T, R>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック

			for(size_t i = 1; i <= N; ++i){
				for(size_t j = i + m - 1; j <= M; ++j){
					Y(j,i) = U(j,i);
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
			static_assert(std::is_convertible_v<T, TL>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			static_assert(std::is_convertible_v<T, TU>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			static_assert(std::is_convertible_v<T, TP>, "ArcsMat: Type Conversion Error");	// 暗黙の型変換可能チェック
			
			// 中間変数
			ArcsMat<M,N,T> X = A;		// 行列操作用にコピー
			size_t perm_count = 0;		// 入れ替えカウンタ
			T max_buff = 0;				// 最大値バッファ
			P = ArcsMat<M,N,T>::eye();	// 置換行列の準備

			// 列ごとに処理を実施
			size_t k = 0;
			for(size_t i = 1; i <= N; ++i){
				// i列目の中での最大値を探す
				k = i;									// 対角要素の縦位置で初期化
				max_buff = std::abs(X(i,i));			// 対角要素の値で初期化
				for(size_t j = i + 1; j <= M; ++j){
					if(max_buff < std::abs( X(j,i) )){	// スキャン中における最大値か判定
						k = j;							// 最大値の候補の縦位置
						max_buff = std::abs( X(j,i) );	// 最大値の候補であればその値を記憶
					}
				}

				// 対角要素が最大値でなければ、対角要素が最大となるように行丸ごと入れ替え
				if(k != i){
					swaprow(X, i, k);	// 対角要素の行と最大値の行を入れ替え
					swaprow(P, i, k);	// 置換行列も同じ様に入れ替え
					perm_count++;		// 入れ替えカウンタ
				}

				// LU分解のコア部分
				if( std::abs( X(i,i) ) < EPSILON ){
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
			gettrilo(X, L);	// 下三角のみを抽出
			gettriup(X, U);	// 上三角のみを抽出
			for(size_t j = 1; j <= M; ++j) L(j,j) = 1;	// 下三角行列の対角要素はすべて1
			
			// 入れ替え回数の判定と判定結果の保持
			if(perm_count % 2 == 0){
				P.Status = MatStatus::AMT_LU_EVEN;	// 偶数のとき
			}else{
				P.Status = MatStatus::AMT_LU_ODD;	// 奇数のとき
			}
		}

		//! @brief LU分解の結果と置換行列を返す関数(タプル返し版)
		//! @param[in]	A	入力行列
		//! @return	(L, U, P)	(下三角行列, 上三角行列, 置換行列)のタプル
		static constexpr std::tuple<ArcsMat<M,N,T>, ArcsMat<M,N,T>, ArcsMat<M,N,T>> LUP(const ArcsMat<M,N,T>& A){
			ArcsMat<M,N,T> L, U, P;
			LUP(A, L, U, P);
			return {L, U, P};
		}
		
		//! @brief LU分解の結果のみ返す関数(引数渡し版)
		//! @tparam	ML, NL, TL, MU, NU, TU, MP, NP, TP	L,U,P行列の高さ, 幅, 要素の型
		//! @param[in]	A	入力行列
		//! @param[out]	L	下三角行列
		//! @param[out]	U	上三角行列
		template<size_t ML, size_t NL, typename TL = double, size_t MU, size_t NU, typename TU = double>
		static constexpr void LU(const ArcsMat<M,N,T>& A, ArcsMat<ML,NL,TL>& L, ArcsMat<MU,NU,TU>& U){
			ArcsMat<M,N,T> P;
			LUP(A, L, U, P);
			L = tp(P)*L;
		}

		//! @brief LU分解の結果のみ返す関数(タプル返し版)
		//! @param[in]	A	入力行列
		//! @return	(L, U)	(下三角行列, 上三角行列)のタプル
		static constexpr std::tuple<ArcsMat<M,N,T>, ArcsMat<M,N,T>> LU(const ArcsMat<M,N,T>& A){
			ArcsMat<M,N,T> L, U;
			LU(A, L, U);
			return {L, U};
		}

/*		
		
		//! @brief 行列の非ゼロ要素数を返す関数
		//! @param[in]	U	入力行列
		//! @return 結果
		constexpr friend size_t nonzeroele(const ArcsMat& U){
			size_t ret = 0;
			for(size_t i = 0; i < U.N; ++i){
				for(size_t j = 0; j < U.M; ++j){
					if(epsilon < std::abs(U.Data[i][j])) ++ret;
				}
			}
			return ret;
		}
		
		//! @brief 行列のランクを返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend size_t rank(const ArcsMat& A){
			static_assert(N == M, "ArcsMat Size Error");	// 正方行列のみ対応
			ArcsMat<N,N,T> U, S, V;
			SVD(A, U, S, V);			// まず特異値分解して，
			return nonzeroele(diag(S));	// S行列の対角要素の非ゼロをカウントするとそれがランク
		}
		
		//! @brief 行列の無限大ノルムを返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend T infnorm(const ArcsMat& U){
			return max(sumcolumn(abse(tp(U))));
		}
		
		//! @brief ベクトルのユークリッドノルムを返す関数
		//! @param[in]	v	入力ベクトル
		//! @return	結果
		constexpr friend T euclidnorm(const ArcsMat<N,M,T>& v){
			ArcsMat<1,1,T> ret;
			if constexpr(std::is_same_v<T,std::complex<double>>){
				// 複素数型の場合
				if constexpr(N == 1){
					// 縦ベクトルの場合
					ret = Htp(v)*v;
				}else if constexpr(M == 1){
					// 横ベクトルの場合
					ret = v*Htp(v);
				}
			}else{
				// 実数型の場合
				if constexpr(N == 1){
					// 縦ベクトルの場合
					ret = tp(v)*v;
				}else if constexpr(M == 1){
					// 横ベクトルの場合
					ret = v*tp(v);
				}else{
					// 行列の場合
					ret = sumcolumn(sumrow(v & v));
				}
			}
			return std::sqrt(ret[1]);
		}
		
		//! @brief 修正コレスキー分解(LDL^T版)
		//! @param[in]	A	入力行列
		//! @param[out]	L	下三角行列
		//! @param[out]	D	対角行列
		constexpr friend void Cholesky(const ArcsMat& A, ArcsMat& L, ArcsMat& D){
			static_assert(A.N == A.M, "ArcsMat Size Error");	// Aが正方行列かチェック
			L.Data[0][0] = A.Data[0][0];
			D.Data[0][0] = 1.0/L.Data[0][0];
			
			for(size_t i = 1; i < A.N; ++i){
				for(size_t j = 0; j <= i; ++j){
					T lld = A.Data[j][i];
					for(size_t k = 0; k < j; ++k){
						lld -= L.Data[k][i]*L.Data[k][j]*D.Data[k][k];
					}
					L.Data[j][i] = lld;
				}
				D.Data[i][i] = 1.0/L.Data[i][i];
			}
		}
		
		//! @brief 修正コレスキー分解(LL^T版)
		//! @param[in]	A	入力行列
		//! @param[out]	L	下三角行列
		constexpr friend void Cholesky(const ArcsMat& A, ArcsMat& L){
			ArcsMat Lp, Dp;
			Cholesky(A, Lp, Dp);// まずコレスキー分解してから，
			L = Lp*sqrte(Dp);	// 対角行列の平方根を取って，下三角行列に掛けたものを出力
		}
		
		//! @brief QR分解
		//! 補足：実数型のときMATLABとはQとRの符号関係が逆の場合があるが正常なQR分解であることは確認済み
		//! 補足：複素数型のときMATLABとは全く違う値が得られるが，正常なQR分解であることは確認済み
		//! @param[in]	A	入力行列
		//! @param[out]	Q	ユニタリ行列 Q行列
		//! @param[out] R	上三角行列 R行列
		constexpr friend void QR(const ArcsMat<N,M,T>& A, ArcsMat<M,M,T>& Q, ArcsMat<N,M,T>& R){
			constexpr size_t K = std::min(N,M);
			ArcsMat<1,M,T> a;
			ArcsMat<1,M,T> e;
			ArcsMat<1,M,T> v;
			ArcsMat<M,M,T> I = ArcsMat<M,M,T>::eye();
			ArcsMat<M,M,T> H;
			ArcsMat<N,M,T> HA;
			HA = A;
			Q = I;
			
			if constexpr(std::is_same_v<T, std::complex<double>>){
				// Householder Reflections を使ったQR分解アルゴリズム(複素数版)
				ArcsMat<1,1,T> vHv;
				e[1] = std::complex(1.0, 0.0);
				for(size_t k = 1; k <= K; ++k){
					a = getcolumn(HA, k);
					a = shiftup(a, k - 1);
					v = a + sgn(a[1])*euclidnorm(a)*e;
					vHv = Htp(v)*v;
					if(k != 1) I.SetElement(M - (k - 2), M - (k - 2), 0);
					if(std::abs(vHv[1]) < epsilon) vHv[1] = epscomp;
					H = I - 2.0*v*Htp(v)/vHv[1];
					H = shiftdown(H, k - 1);
					H = shiftright(H, k - 1);
					for(size_t i = 1; i < k; ++i) H.SetElement(i,i, 1);
					HA = H*HA;
					Q = Q*H;
				}
				R = HA;
			}else{
				// Householder Reflections を使ったQR分解アルゴリズム(実数版)
				ArcsMat<1,1,T> vTv;
				e[1] = 1;
				for(size_t k = 1; k <= K; ++k){
					a = getcolumn(HA, k);
					a = shiftup(a, k - 1);
					v = a + sgn(a[1])*euclidnorm(a)*e;
					vTv = tp(v)*v;
					if(k != 1) I.SetElement(M - (k - 2), M - (k - 2), 0);
					if(       0 <= vTv[1] && vTv[1] < epsilon) vTv[1] =  epsilon;
					if(-epsilon <= vTv[1] && vTv[1] < 0      ) vTv[1] = -epsilon;
					H = I - 2.0*v*tp(v)/vTv[1];
					H = shiftdown(H, k - 1);
					H = shiftright(H, k - 1);
					for(size_t i = 1; i < k; ++i) H.SetElement(i,i, 1);
					HA = H*HA;
					Q = Q*H;
				}
				R = HA;
			}
		}
		
		//! @brief SVD特異値分解(引数渡し版)
		//! 補足：MATLABとはU,S,Vの符号関係が入れ替わっている場合があるが正常なSVDであることは確認済み
		//! @param[in]	A	入力行列
		//! @param[out]	U	U行列
		//! @param[out]	S	S行列
		//! @param[out]	V	V行列
		constexpr friend void SVD(const ArcsMat<N,M,T>& A, ArcsMat<M,M,T>& U, ArcsMat<N,M,T>& S, ArcsMat<N,N,T>& V){
			constexpr size_t LoopMax = 100*std::max(N,M);	// ループ打ち切り最大回数
			ArcsMat<N,M,T> Snm;
			ArcsMat<M,N,T> Smn;
			ArcsMat<M,M,T> Qm;
			ArcsMat<N,N,T> Qn;
			ArcsMat<M,N,T> e;
			double E = 0, F = 0;
			
			// 初期化
			U = ArcsMat<M,M,T>::eye();
			Snm = A;
			V = ArcsMat<N,N,T>::eye();
			T Error = 1;
			
			// ループ打ち切り最大回数に達するまでループ
			for(size_t l = 0; l < LoopMax; ++l){
				if(Error < epsilon) break;	// 誤差がイプシロンを下回ったらループ打ち切り
				QR(Snm, Qm, Snm);
				U = U*Qm;
				QR(tp(Snm), Qn, Smn);
				V = V*Qn;
				e = gettriup(Smn, 1);
				E = euclidnorm(e);
				F = euclidnorm(diag(Smn));
				if(-epsilon < F && F < epsilon) F = 1;
				Error = E/F;
				Snm = tp(Smn);
			}
			
			// 符号修正
			ArcsMat<1,std::min(N,M),T> Sd = diag(Smn);
			S = ArcsMat<N,M,T>::zeros();
			for(size_t k = 1; k <= std::min(N,M); ++k){
				T Sdn = Sd[k];
				S.SetElement(k,k, std::abs(Sdn));
				if(Sdn < 0){
					setcolumn(U, -getcolumn(U, k), k);
				}
			}
		}
		
		//! @brief SVD特異値分解(タプルで返す版)
		//! 補足：MATLABとはU,S,Vの符号関係が入れ替わっている場合があるが正常なSVDであることは確認済み
		//! @param[in]	A	入力行列
		//! @return [U行列, S行列, V行列]のタプル
		constexpr friend std::tuple<ArcsMat<M,M,T>, ArcsMat<N,M,T>, ArcsMat<N,N,T>> SVD(const ArcsMat<N,M,T>& A){
			ArcsMat<M,M,T> U;
			ArcsMat<N,M,T> S;
			ArcsMat<N,N,T> V;
			SVD(A, U, S, V);	// SVD特異値分解
			return {U, S, V};	// タプルで返す
		}
		
		//! @brief Schur分解
		//! @param[in]	A	入力行列
		//! @param[out]	Q	ユニタリ行列
		//! @param[out]	U	上三角行列または疑似上三角行列
		constexpr friend void Schur(const ArcsMat<N,M,T>& A, ArcsMat<N,M,T>& Q, ArcsMat<N,M,T>& U){
			static_assert(A.N == A.M, "ArcsMat Size Error");				// Aが正方行列かチェック
			
			// QR法によるシュール分解
			ArcsMat<N,M,T> Tk = A, Qk, Rk, Qa;
			Qa = ArcsMat<N,M,T>::eye();
			for(size_t k = 0; k < ITERATION_MAX; ++k){
				QR(Tk, Qk, Rk);
				Tk = Rk*Qk;
				Qa = Qa*Qk;
				if(std::abs(Tk.GetElement(1,A.M)) < A.epsilon) break;
			}
			Q = Qa;
			U = Tk;
		}
		
		//! @brief Ax = bの形の線形連立1次方程式をxについて解く関数(引数渡し版)
		//! @param[in]	A	係数行列
		//! @param[in]	b	係数ベクトル
		//! @param[out]	x	解ベクトル
		constexpr friend void solve(const ArcsMat& A, const ArcsMat<1,M,T>& b, ArcsMat<1,N,T>& x){
			static_assert(A.N == A.M, "ArcsMat Size Error");				// Aが正方行列かチェック
			static_assert(b.M == A.M, "ArcsMat and Vector Size Error");	// Aの高さとbの高さが同じかチェック
			static_assert(b.N == 1, "Input is NOT vector.");			// bは縦ベクトルかチェック
			
			if constexpr(M == 1){
				// スカラーの場合
				x[1] = b[1]/A.GetElement(1,1);	// スカラーのときはそのまま単純に除算
			}else{
				// 行列の場合
				// Ax = b において A をLU分解すると，(LU)x = b になって L(Ux) = b と表現できることを利用する。
				ArcsMat<A.N,A.N,T> L, U;
				ArcsMat<1,A.N,T> d, bb;
				ArcsMat<1,A.N,int> v;
				T buff = 0;
				LU(A, L, U, v);		// まず，LU分解(並べ替え有り)してから，Ux = d と勝手に置き換えて考える。
				bb = orderrow(b, v);// bベクトルも並べ替える
				// その次に Ld = b の方を d について解いて，
				// (下記では，Lの対角要素が1となるLU分解がなされているものとして計算する)
				d.Data[0][0] = bb.Data[0][0];
				for(size_t i = 1; i < A.N; ++i){
					for(size_t j = 0; j <= i - 1; ++j) buff += L.Data[j][i]*d.Data[0][j];
					d.Data[0][i] = (bb.Data[0][i] - buff);
					buff = 0;
				}
				// さらに Ux = d を x について解く。
				x.Data[0][A.N-1] = d.Data[0][A.N-1]/U.Data[A.N-1][A.N-1];
				for(int k = A.N - 2; 0 <= k; --k){
					for(size_t j = (size_t)k + 1; j < A.N; ++j){
						buff += U.Data[j][k]*x.Data[0][j];
					}
					x.Data[0][k] = (d.Data[0][k] - buff)/U.Data[k][k];
					buff = 0;
				}
			}
		}
		
		//! @brief Ax = bの形の線形連立1次方程式をxについて解く関数(戻り値として返す版)
		//! @param[in]	A	係数行列
		//! @param[in]	b	係数ベクトル
		//! @return	解ベクトル
		constexpr friend ArcsMat<1,M,T> solve(const ArcsMat& A, const ArcsMat<1,M,T>& b){
			ArcsMat<1,A.N,T> x;
			solve(A, b, x);
			return x;		// 最終的な答えのxベクトルを返す
		}
		
		//! @brief Uは上三角行列で，Ux = bの形の線形連立1次方程式をxについて解く関数(引数渡し版)
		//! @param[in]	U	係数行列(上三角行列)
		//! @param[in]	b	係数ベクトル
		//! @param[out]	x	解ベクトル
		constexpr friend void solve_upper_tri(const ArcsMat& U, const ArcsMat<1,M,T>& b, ArcsMat<1,N,T>& x){
			static_assert(U.N == U.M, "ArcsMat Size Error");				// Uが正方行列かチェック
			static_assert(b.M == U.M, "ArcsMat and Vector Size Error");	// Uの高さとbの高さが同じかチェック
			static_assert(b.N == 1, "Input is NOT vector.");			// bは縦ベクトルかチェック
			
			if constexpr(M == 1){
				// スカラーの場合
				x[1] = b[1]/U.GetElement(1,1);	// スカラーのときはそのまま単純に除算
			}else{
				// 行列の場合
				ArcsMat<1,U.N,int> v;
				T buff = 0;
				// 既にUは上三角行列なのでLU分解は不要
				// Ux = b を x について解く。
				x.Data[0][U.N-1] = b.Data[0][U.N-1]/U.Data[U.N-1][U.N-1];
				for(int k = U.N - 2; 0 <= k; --k){
					for(size_t j = (size_t)k + 1; j < U.N; ++j){
						buff += U.Data[j][k]*x.Data[0][j];
					}
					x.Data[0][k] = (b.Data[0][k] - buff)/U.Data[k][k];
					buff = 0;
				}
			}
		}
		
		//! @brief 行列式の値を返す関数
		//! @param[in]	A	入力行列
		//! @return	結果
		constexpr friend T det(const ArcsMat& A){
			static_assert(A.N == A.M, "ArcsMat Size Error");	// Aが正方行列かチェック
			ArcsMat<A.N,A.N,T> L, U;
			ArcsMat<1,A.N,int> v;
			int sign;		// 符号
			if(LU(A, L, U, v) == ArcsMat::ODD){	// LU分解と符号判定
				sign = -1;	// 奇数のとき
			}else{
				sign =  1;	// 偶数のとき
			}
			// |A| = |L||U| でしかも |L|と|U|は対角要素の総積に等しく，さらにLの対角要素は1なので|L|は省略可。
			// 最後にLU分解のときの並べ替え回数によって符号反転をする。
			return (T)sign*prod(U);
		}
		
		//! @brief 逆行列を返す関数 (正則チェック無し)
		//! @param[in]	A	入力行列
		//! @return	結果
		constexpr friend ArcsMat inv(const ArcsMat& A){
			static_assert(A.N == A.M, "ArcsMat Size Error");	// Aが正方行列かチェック
			ArcsMat I = ArcsMat<A.N,A.N,T>::ident();			// 単位行列の生成
			ArcsMat<1,A.N,T> x, b;
			ArcsMat<A.N,A.N,T> Ainv;
			for(size_t n = 1; n <= A.N; ++n){
				b = getcolumn(I, n);	// 単位行列のn列目を切り出してbベクトルとする
				solve(A, b, x);			// Ax = b の連立1次方程式をxについて解く
				setcolumn(Ainv, x, n);	// xはAの逆行列のn列目となるので、Ainvのn列目にxを書き込む
			}
			return Ainv;	// 最終的に得られる逆行列を返す
		}
		
		//! @brief 逆行列を返す関数 (正則チェック無し, 左上小行列のサイズ指定版)
		//! @param[in]	A	入力行列 (kより右と下は全部ゼロ埋めを想定)
		//! @param[in]	k	左上小行列のサイズ
		//! @return	結果
		constexpr friend ArcsMat inv(const ArcsMat& A, size_t k){
			arcs_assert(k <= A.N);					// kが範囲内かチェック
			ArcsMat I = ArcsMat<A.N,A.N,T>::ident();	// 単位行列の生成
			
			// 正則にするためにk列より右下の対角成分を「1」にする
			ArcsMat<A.N,A.N,T> A2 = A;
			for(size_t j = k + 1; j <= A.N; ++j){
				A2.SetElement(j, j, 1);
			}
			
			// k列までの逆行列を計算
			ArcsMat<1,A.N,T> x, b;
			ArcsMat<A.N,A.N,T> Ainv;
			for(size_t n = 1; n <= k; ++n){
				b = getcolumn(I, n);	// 単位行列のn列目を切り出してbベクトルとする
				solve(A2, b, x);		// Ax = b の連立1次方程式をxについて解く
				setcolumn(Ainv, x, n);	// xはAの逆行列のn列目となるので、Ainvのn列目にxを書き込む
			}
			return Ainv;	// 最終的に得られる逆行列を返す
		}
		
		//! @brief 逆行列を返す関数 (正則チェック有り)
		//! @param[in]	A	入力行列
		//! @return	結果
		constexpr friend ArcsMat inv_with_check(const ArcsMat& A){
			static_assert(A.N == A.M, "ArcsMat Size Error");	// Aが正方行列かチェック
			arcs_assert(A.epsilon < std::abs(det(A)));		// 正則かチェック
			return inv(A);	// 最終的に得られる逆行列を返す
		}
		
		//! @brief 上三角行列の逆行列を返す関数
		//! @param[in]	U	入力行列(上三角行列)
		//! @param[out]	Uinv	逆行列
		constexpr friend void inv_upper_tri(const ArcsMat& U, ArcsMat& Uinv){
			static_assert(U.N == U.M, "ArcsMat Size Error");			// Uが正方行列かチェック
			static_assert(Uinv.N == Uinv.M, "ArcsMat Size Error");	// Uinvが正方行列かチェック
			static_assert(U.N == Uinv.N, "ArcsMat Size Error");		// UとUinvが同じサイズかチェック
			ArcsMat I = ArcsMat<U.N,U.N,T>::ident();			// 単位行列の生成
			ArcsMat<1,U.N,T> x, b;
			for(size_t n = 1; n <= U.N; ++n){
				b = getcolumn(I, n);		// 単位行列のn列目を切り出してbベクトルとする
				solve_upper_tri(U, b, x);	// Ux = b の連立1次方程式をxについて解く
				setcolumn(Uinv, x, n);		// xはUの逆行列のn列目となるので、Uinvのn列目にxを書き込む
			}
			// Uinvが最終的に得られる逆行列
		}
		
		//! @brief 上三角行列の逆行列を返す関数(左上小行列のサイズ指定版)
		//! @param[in]	U	入力行列(上三角行列, kより右と下は全部ゼロ埋めを想定)
		//! @param[in]	k	左上小行列のサイズ
		//! @param[out]	Uinv	逆行列
		constexpr friend void inv_upper_tri(const ArcsMat& U, size_t k, ArcsMat& Uinv){
			static_assert(U.N == U.M, "ArcsMat Size Error");			// Uが正方行列かチェック
			static_assert(Uinv.N == Uinv.M, "ArcsMat Size Error");	// Uinvが正方行列かチェック
			static_assert(U.N == Uinv.N, "ArcsMat Size Error");		// UとUinvが同じサイズかチェック
			ArcsMat I = ArcsMat<U.N,U.N,T>::ident();			// 単位行列の生成
			
			// 正則にするためにk列より右下の対角成分を「1」にする
			ArcsMat<U.N,U.N,T> U2 = U;
			for(size_t j = k + 1; j <= U.N; ++j){
				U2.SetElement(j, j, 1);
			}
			
			// k列までの逆行列を計算
			ArcsMat<1,U.N,T> x, b;
			for(size_t n = 1; n <= k; ++n){
				b = getcolumn(I, n);		// 単位行列のn列目を切り出してbベクトルとする
				solve_upper_tri(U2, b, x);	// Ux = b の連立1次方程式をxについて解く
				setcolumn(Uinv, x, n);		// xはUの逆行列のn列目となるので、Uinvのn列目にxを書き込む
			}
			// Uinvが最終的に得られる逆行列
		}
		
		//! @brief 左擬似逆行列を返す関数 (Aが縦長行列の場合)
		//! @param[in]	A	入力行列
		//! @return	結果
		constexpr friend ArcsMat<M,N,T> lpinv(const ArcsMat& A){
			static_assert(A.N < A.M, "ArcsMat Size Error");	// 縦長行列かチェック
			return inv(tp(A)*A)*tp(A);
		}
		
		//! @brief 左擬似逆行列を返す関数 (Aが縦長行列の場合, 左上小行列のサイズ指定版)
		//! @param[in]	A	入力行列 (kより右と下は全部ゼロ埋めを想定)
		//! @param[in]	k	左上小行列のサイズ
		//! @return	結果
		constexpr friend ArcsMat<M,N,T> lpinv(const ArcsMat& A, size_t k){
			ArcsMat<M,N> At = tp(A);
			ArcsMat<N,N> A2 = At*A;
			return inv(A2, k)*At;
		}
		
		//! @brief 右擬似逆行列を返す関数 (Aが横長行列の場合)
		//! @param[in]	A	入力行列
		//! @return	結果
		constexpr friend ArcsMat<M,N,T> rpinv(const ArcsMat& A){
			static_assert(A.M < A.N, "ArcsMat Size Error");	// 横長行列かチェック
			return tp(A)*inv(A*tp(A));
		}
		
		//! @brief 右擬似逆行列を返す関数 (Aが横長行列の場合, 左上小行列のサイズ指定版)
		//! @param[in]	A	入力行列 (kより右と下は全部ゼロ埋めを想定)
		//! @param[in]	k	左上小行列のサイズ
		//! @return	結果
		constexpr friend ArcsMat<M,N,T> rpinv(const ArcsMat& A, size_t k){
			ArcsMat<M,N> At = tp(A);
			ArcsMat<M,M> A2 = A*At;
			return At*inv(A2, k);
		}
		
		//! @brief 行列指数関数 e^(U)
		//! @param[in]	U	入力行列
		//! @param[in]	Order	パデ近似の次数
		//! @return	結果
		constexpr friend ArcsMat expm(const ArcsMat& U, size_t Order){
			static_assert(U.N == U.M, "ArcsMat Size Error");	// 正方行列かチェック
			int e = 0;
			T c = 1;
			bool flag = false;
			// ノルムでスケーリング
			frexp(infnorm(U),&e);
			ArcsMat<U.N,U.M,T> A;
			if(0 < e){
				A = pow(0.5,e+1)*U;
			}else{
				e = 0;
				A = 0.5*U;
			}
			// 行列のパデ近似の計算
			ArcsMat<A.N,A.N,T> I = ident();// 単位行列の生成
			ArcsMat<A.N,A.N,T> L = I, R = I, X = I, cX;
			for(size_t i = 1; i <= Order; ++i){
				c = c*(T)(Order - i + 1)/(T)(i*(2*Order - i + 1));	// パデ近似係数の計算
				X = A*X;		// A^Mの計算
				cX = c*X;		// cM*A^Mの計算
				R += cX;		// R = I + c1*A + c2*A*A + c3*A*A*A + ... + cM*A^M
				if(flag == true){
					L += cX;	// L = I + c1*A + c2*A*A + c3*A*A*A + ... + cM*A^M の正の係数の場合
				}else{
					L -= cX;	// L = I + c1*A + c2*A*A + c3*A*A*A + ... + cM*A^M の負の係数の場合
				}
				flag = !flag;	// 正負係数の場合分け用フラグ
			}
			// スケールを元に戻す
			ArcsMat Y = inv(L)*R;
			for(size_t i = 0; i < (size_t)e + 1; ++i){
				Y = Y*Y;
			}
			return Y;	// 最終的に得られる行列指数を返す
		}
		
		//! @brief 指数行列の数値定積分[0,T]をする関数
		//! @param[in]	U	入力行列
		//! @param[in]	T	積分範囲の終わり
		//! @param[in]	DIV	分割数
		//! @param[in]	P	パデ近似の次数
		//! @return	結果
		constexpr friend ArcsMat integral_expm(const ArcsMat& U, const T T, const size_t DIV, const size_t P){
			static_assert(U.N == U.M, "ArcsMat Size Error");	// 正方行列かチェック
			const T h = T/((T)(2*DIV));	// 時間ステップ
			T t = 0;						// 時刻
			// シンプソン法による定積分の実行
			ArcsMat<U.N,U.M,T> S1, S2;
			for(size_t i = 1; i <= DIV; ++i){
				t = h*(T)(2*i - 1);
				S1 += expm(U*t, P);
			}
			for(size_t i = 1; i <= DIV - 1; ++i){
				t = h*(T)(2*i);
				S2 += expm(U*t, P);
			}
			return h/3.0*( ArcsMat<U.N,U.N,T>::eye() + 4.0*S1 + 2.0*S2 + expm(U*T,P) );	// 最終的な定積分結果を返す
		}
		
		
		
		//! @brief 複素数行列要素の実数部を返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend ArcsMat<N,M,double> reale(const ArcsMat& U){
			static_assert(std::is_same_v<T, std::complex<double>>, "ArcsMat Type Error");	// 複素数型のみ対応
			ArcsMat<N,M,double> Y;
			for(size_t i = 0; i < U.N; ++i){
				for(size_t j = 0; j < U.M; ++j) Y.Data[i][j] = U.Data[i][j].real();
			}
			return Y;
		}
		
		//! @brief 複素数行列要素の実数部に値をセットする関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr void real(const ArcsMat<N,M,double>& U){
			static_assert(N == N, "ArcsMat Size Error");
			static_assert(M == M, "ArcsMat Size Error");
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) Data[i][j] = std::complex(U.Data[i][j], 0.0);
			}
		}
		
		//! @brief 複素数行列要素の虚数部を返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend ArcsMat<N,M,double> image(const ArcsMat& U){
			static_assert(std::is_same_v<T, std::complex<double>>, "ArcsMat Type Error");	// 複素数型のみ対応
			ArcsMat<N,M,double> Y;
			for(size_t i = 0; i < U.N; ++i){
				for(size_t j = 0; j < U.M; ++j) Y.Data[i][j] = U.Data[i][j].imag();
			}
			return Y;
		}
		
		//! @brief 複素数行列要素の大きさを返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend ArcsMat<N,M,double> mage(const ArcsMat& U){
			static_assert(std::is_same_v<T, std::complex<double>>, "ArcsMat Type Error");	// 複素数型のみ対応
			ArcsMat<N,M,double> Y;
			for(size_t i = 0; i < U.N; ++i){
				for(size_t j = 0; j < U.M; ++j) Y.Data[i][j] = std::abs(U.Data[i][j]);
			}
			return Y;
		}
		
		//! @brief 複素数行列要素の偏角を返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend ArcsMat<N,M,double> arge(const ArcsMat& U){
			static_assert(std::is_same_v<T, std::complex<double>>, "ArcsMat Type Error");	// 複素数型のみ対応
			ArcsMat<N,M,double> Y;
			for(size_t i = 0; i < U.N; ++i){
				for(size_t j = 0; j < U.M; ++j) Y.Data[i][j] = std::arg(U.Data[i][j]);
			}
			return Y;
		}
		
		//! @brief 複素数行列要素の共役を返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend ArcsMat<N,M,std::complex<double>> conje(const ArcsMat& U){
			static_assert(std::is_same_v<T, std::complex<double>>, "ArcsMat Type Error");	// 複素数型のみ対応
			ArcsMat<N,M,std::complex<double>> Y;
			for(size_t i = 0; i < U.N; ++i){
				for(size_t j = 0; j < U.M; ++j) Y.Data[i][j] = std::conj(U.Data[i][j]);
			}
			return Y;
		}
		
		//! @brief エルミート転置行列を返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend ArcsMat<M,N,std::complex<double>> Htp(const ArcsMat<N,M,T>& U){
			static_assert(std::is_same_v<T, std::complex<double>>, "ArcsMat Type Error");	// 複素数型のみ対応
			return conje(tp(U));	// 転置して複素共役
		}
		
		//! @brief 固有値を返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend ArcsMat<1,N,std::complex<double>> eigen(const ArcsMat<N,M,T>& U){
			static_assert(N == M, "ArcsMat Size Error");		// 正方行列のみ対応
			constexpr size_t LoopMax = 100*N;					// ループ打ち切り最大回数
			auto I = ArcsMat<N,N,std::complex<double>>::eye();	// 単位行列
			ArcsMat<N,N,std::complex<double>> A, Q, R;
			std::complex<double> a, b, c, d, mu;
			
			if constexpr(std::is_same_v<T, std::complex<double>>){
				// 入力が複素数型の場合
				A = U;
			}else{
				// 入力が実数型の場合
				A.real(U);
			}
			
			// 複素数QR法による固有値計算
			for(size_t k = 1; k < LoopMax; ++k){
				// Aの最右下2×2小行列の固有値を求める
				a = A.GetElement(N - 1, N - 1);
				b = A.GetElement(N    , N - 1);
				c = A.GetElement(N - 1, N    );
				d = A.GetElement(N    , N    );
				mu = ( (a + d) + std::sqrt((a + d)*(a + d) - 4.0*(a*d - b*c)) )/2.0;
				
				// QR分解と収束計算
				QR(A - mu*I, Q, R);
				A = R*Q + mu*I;
				
				if(std::abs(std::abs(tr(Q)) - (double)N) < epsilon) break;	// 単位行列に漸近したらループ打ち切り
			}
			
			return diag(A);	// Aの対角要素が固有値
		}
		
		//! @brief 最大固有値の固有ベクトルを返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend ArcsMat<1,N,std::complex<double>> eigenvec(const ArcsMat<N,M,T>& U){
			static_assert(N == M, "ArcsMat Size Error");	// 正方行列のみ対応
			constexpr size_t LoopMax = 100*std::max(N,M);	// ループ打ち切り最大回数
			ArcsMat<N,N,std::complex<double>> A;
			
			if constexpr(std::is_same_v<T, std::complex<double>>){
				// 入力が複素数型の場合
				A = U;
			}else{
				// 入力が実数型の場合
				A.real(U);
			}
			
			// べき乗法による固有ベクトル計算
			auto x = ArcsMat<1,N,std::complex<double>>::ones();
			auto y = ArcsMat<1,N,std::complex<double>>::ones();
			for(size_t k = 1; k < LoopMax; ++k){
				y = A*x;
				x = y/euclidnorm(y);
			}
			
			return x;
		}
		
		//! @brief クロネッカー積
		//! @tparam	P	右側の幅
		//! @tparam Q	右側の高さ
		//! @param[in] Ul 演算子の左側
		//! @param[in] Ur 演算子の右側
		//! @return 結果
		template<size_t P, size_t Q>
		constexpr friend ArcsMat<N*P, M*Q, T> Kronecker(const ArcsMat<N,M,T>& Ul, const ArcsMat<P,Q,T>& Ur){
			ArcsMat<N*P, M*Q, T> Y;
			ArcsMat<P,Q,T> A;
			
			// 縦方向に小行列で埋めていく
			for(size_t j = 1; j <= M; ++j){
				// 横方向に小行列で埋めていく
				for(size_t i = 1; i <= N; ++i){
					A = Ul(j,i)*Ur;
					setsubmatrix(Y, (i - 1)*P + 1, (j - 1)*Q + 1, A);
				}
			}
			
			return Y;
		}
		
		//! @brief vec作用素(行列→縦ベクトル)
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend ArcsMat<1,N*M,T> vec(const ArcsMat<N,M,T>& U){
			ArcsMat<1,N*M,T> Y;
			
			// 横方向に走査
			size_t k = 0;
			for(size_t i = 1; i <= N; ++i){
				// 縦方向に走査
				for(size_t j = 1; j <= M; ++j){
					k++;
					Y(k,1) = U(j,i);
				}
			}
			
			return Y;
		}
		
		//! @brief vec作用素の逆(縦ベクトル→行列)
		//! @tparam	P	再構成後の幅
		//! @tparam Q	再構成後の高さ
		//! @param[in]		v	入力ベクトル
		//! @param[in,out]	Y	再構成後の行列
		template<size_t P, size_t Q>
		constexpr static void vecinv(const ArcsMat<N,M,T>& v, ArcsMat<P,Q,T>& Y){
			static_assert(N == 1);		// 入力は縦ベクトルかチェック
			static_assert(P*Q == M);	// 要素数が同じかチェック
			
			// 横方向に走査
			size_t k = 0;
			for(size_t i = 1; i <= P; ++i){
				// 縦方向に走査
				for(size_t j = 1; j <= Q; ++j){
					k++;
					Y(j,i) = v(k,1);
				}
			}
		}
		
		//! @brief vec作用素の逆(縦ベクトル→行列)
		//! @tparam	P	再構成後の幅
		//! @tparam Q	再構成後の高さ
		//! @param[in]		v	入力ベクトル
		//! @return	再構成後の行列
		template<size_t P, size_t Q>
		constexpr static ArcsMat<P,Q,T> vecinv(const ArcsMat<N,M,T>& v){
			ArcsMat<P,Q,T> Y;
			vecinv(v, Y);
			return Y;
		}
*/		
	private:
		// 基本定数
		static constexpr double EPSILON = 1e-12;		//!< 零とみなす閾値(実数版)
		static constexpr std::complex<double> EPSLCOMP = std::complex(1e-12, 1e-12);	//!< 零とみなす閾値(複素数版)
		static constexpr size_t ITERATION_MAX = 10000;	//!< 反復計算の最大値
		
		// 内部処理用
		size_t Nindex;	//!< 横方向カウンタ
		size_t Mindex;	//!< 縦方向カウンタ
		MatStatus Status = MatStatus::AMT_NA;	// 行列の状態

		// 行列データの実体
		// ArcsMatは行列の縦方向にメモリアドレスが連続しているので、縦ベクトル優先。
		std::array<std::array<T, M>, N> Data;//!< データ格納用変数 配列要素の順番は Data[N列(横方向)][M行(縦方向)]
/*		
		//! @brief 符号関数
		//! @param[in]	u	入力
		//! @return	符号結果
		static constexpr T sgn(T u){
			T ret = 0;
			if constexpr(std::is_same_v<T, std::complex<double>>){
				// 複素数型の場合
				if(u == std::complex(0.0, 0.0)){
					ret = std::complex(0.0, 0.0);
				}else{
					ret = u/std::abs(u);
				}
			}else{
				// 実数型の場合
				if((T)0 <= u){
					ret = (T)1;
				}else{
					ret = (T)(-1);
				}
			}
			return ret;
		}
*/		
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
	static void dispmatfmt_macro(const ArcsMat<M,N,T>& U, const std::string& format, const std::string& varname){
		printf("%s = \n", varname.c_str());
		U.Disp(format);
	}
	
	//! @brief 行列の要素を表示 (書式指定なし版，この関数はマクロを介して呼ばれることを想定している)
	//! @tparam	M, N, T	行列の高さ, 幅, 要素の型
	//! @param[in] U		表示する行列
	//! @param[in] varname	変数名
	template <size_t M, size_t N, typename T = double>
	static void dispmat_macro(const ArcsMat<M,N,T>& U, const std::string& varname){
		// データ型によって書式指定子を変える
		// 浮動小数点型の場合
		if constexpr(std::is_floating_point_v<T>){
			dispmatfmt_macro(U, "% g", varname);
			return;
		}
		
		// int型の場合
		if constexpr(std::is_same_v<T, int>){
			dispmatfmt_macro(U, "% d", varname);
			return;
		}
		
		// long型の場合
		if constexpr(std::is_same_v<T, long>){
			dispmatfmt_macro(U, "% ld", varname);
			return;
		}
		
		// size_t型の場合
		if constexpr(std::is_same_v<T, size_t>){
			dispmatfmt_macro(U, "% zu", varname);
			return;
		}
		
		// 複素数型の場合
		if constexpr(std::is_same_v<T, std::complex<double>> || std::is_same_v<T, std::complex<float>>){
			dispmatfmt_macro(U, "% g", varname);
			return;
		}
	}

	//! @brief MM行NN列の単位行列を返す関数
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
	//! @param[out]	y	出力ベクトル
	template<size_t M, size_t N, typename T = double, size_t P, size_t Q, typename R = double, size_t L = std::min(M,N)>
	constexpr void getdiag(const ArcsMat<M,N,T>& U, ArcsMat<P,Q,R>& y){
		ArcsMat<M,N,T>::getdiag(U, y);
	}
	
	//! @brief 行列の対角要素を縦ベクトルとして取得する関数(戻り値渡し版)
	//! @tparam	M, N, T, L	入力行列の高さ, 幅, 要素の型, 対角要素の数
	//! @param[in]	U	入力行列
	//! @return	出力ベクトル
	template<size_t M, size_t N, typename T = double, size_t L = std::min(M,N)>
	constexpr ArcsMat<L,1,T> getdiag(const ArcsMat<M,N,T>& U){
		return ArcsMat<M,N,T>::getdiag(U);
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
	//! @param[in]	A	入力行列
	//! @return	(L, U)	(下三角行列, 上三角行列)のタプル
	template<size_t M, size_t N, typename T = double>
	constexpr std::tuple<ArcsMat<M,N,T>, ArcsMat<M,N,T>> LU(const ArcsMat<M,N,T>& A){
		return ArcsMat<M,N,T>::LU(A);
	}

}

}


#endif


//! @file Matrix.hh
//! @brief 行列/ベクトル計算クラス(テンプレート版)
//! @date 2022/07/27
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

#ifndef MATRIX
#define MATRIX

#include <string>
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

#define PrintMatSize(a)  (PrintMatSize_Macro((a),#a))	//!< 行列サイズ表示マクロ
#define PrintMatrix(a,b) (PrintMatrix_Macro((a),b,#a))	//!< 行列要素表示マクロ(フォーマット指定あり版)
#define PrintMat(a) (PrintMat_Macro((a),#a))			//!< 行列要素表示マクロ(フォーマット指定なし版)

namespace ARCS {	// ARCS名前空間
//! @brief 行列/ベクトル計算クラス(テンプレート版)
//! @tparam	NN	行列の幅
//! @tparam MM	行列の高さ
//! @tparam TT	データ型(デフォルトはdouble型)
template <size_t NN, size_t MM, typename TT = double>
class Matrix {
	public:
		//! @brief LU分解の際の並べ替えが偶数回/奇数回発生したことを返すための定義
		enum LUperm {
			ODD,	//!< 奇数
			EVEN	//!< 偶数
		};
		
		//! @brief コンストラクタ
		constexpr Matrix()
			: Nindex(0), Mindex(0), Data({0})
		{
			static_assert(NN != 0);	// サイズゼロの行列は禁止
			static_assert(MM != 0);	// サイズゼロの行列は禁止
			FillAll(0);				// すべての要素を零で初期化
		}
		
		//! @brief コンストラクタ(任意初期値版)
		//! @param[in]	InitValue	行列要素の初期値
		constexpr explicit Matrix(const TT InitValue)
			: Nindex(0), Mindex(0), Data({0})
		{
			static_assert(NN != 0);	// サイズゼロの行列は禁止
			static_assert(MM != 0);	// サイズゼロの行列は禁止
			FillAll(InitValue);		// すべての要素を指定した値で初期化
		}
		
		//! @brief コンストラクタ(初期化リスト版)
		//! @param[in]	InitList	初期化リスト
		constexpr Matrix(std::initializer_list<TT> InitList)
			: Nindex(0), Mindex(0), Data({0})
		{
			static_assert(NN != 0);	// サイズゼロの行列は禁止
			static_assert(MM != 0);	// サイズゼロの行列は禁止
			const TT* ListVal = InitList.begin();		// 初期化リストの最初のポインタ位置
			size_t Ni = 0;			// 行カウンタ
			size_t Mi = 0;			// 列カウンタ
			for(size_t i = 0; i < InitList.size(); ++i){
				// 初期化リストを順番に読み込んでいく
				assert(Ni < N);		// 行カウンタが行の長さ以内かチェック
				assert(Mi < M);		// 列カウンタが列の長さ以内かチェック
				Data[Ni][Mi] = (TT)ListVal[i];	// 行列の要素を埋める
				Ni++;				// 行カウンタをインクリメント
				if(Ni == N){		// 行カウンタが最後まで行き着いたら，
					Ni = 0;			// 行カウンタを零に戻して，
					Mi++;			// その代わりに，列カウンタをインクリメント
				}
			}
		}
		
		/*
		//! @brief ムーブコンストラクタ
		Matrix(Matrix&& r)
			: Nindex(0), Mindex(0), Data(std::move(r.Data))
		{
			
		}
		*/
		
		//! @brief コピーコンストラクタ
		constexpr Matrix(const Matrix& right)
			: Nindex(0), Mindex(0), Data(right.Data)
		{
			
		}
		
		/* constexprコンストラクタのための「trivial destructor」
		//! @brief デストラクタ
		~Matrix(){
			
		}
		*/
		
		//! @brief 行列代入演算子
		//! @param[in] right 演算子の右側
		//! @return 結果
		constexpr auto& operator=(const Matrix& right){
			static_assert(this->N == right.N, "Matrix Size Error");	// 行列のサイズチェック
			static_assert(this->M == right.M, "Matrix Size Error");	// 行列のサイズチェック
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) this->Data[i][j] = right.Data[i][j];
			}
			return (*this);
		}
		
		//! @brief 単項プラス演算子
		//! @return 結果
		constexpr auto operator+(void) const{
			Matrix ret;
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) ret.Data[i][j] = Data[i][j];
			}
			return ret;
		}
		
		//! @brief 単項マイナス演算子
		//! @return 結果
		constexpr auto operator-(void) const{
			Matrix ret;
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) ret.Data[i][j] = -Data[i][j];
			}
			return ret;
		}
		
		//! @brief 行列加算演算子 (行列同士の加算の場合)
		//! @param[in] right 演算子の右側
		//! @return 結果
		constexpr auto operator+(const Matrix& right) const{
			static_assert(N == right.N, "Matrix Size Error");	// 行列のサイズチェック
			static_assert(M == right.M, "Matrix Size Error");	// 行列のサイズチェック
			Matrix ret;
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) ret.Data[i][j] = Data[i][j] + right.Data[i][j];
			}
			return ret;
		}
		
		//! @brief 行列加算演算子 (行列＋スカラーの場合)
		//! @param[in] right 演算子の右側
		//! @return 結果
		constexpr auto operator+(const TT& right) const{
			Matrix ret;
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) ret.Data[i][j] = Data[i][j] + right;
			}
			return ret;
		}
		
		//! @brief 行列減算演算子 (行列同士の減算の場合)
		//! @param[in] right 演算子の右側
		//! @return 結果
		constexpr auto operator-(const Matrix& right) const{
			static_assert(N == right.N, "Matrix Size Error");	// 行列のサイズチェック
			static_assert(M == right.M, "Matrix Size Error");	// 行列のサイズチェック
			Matrix ret;
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) ret.Data[i][j] = Data[i][j] - right.Data[i][j];
			}
			return ret;
		}
		
		//! @brief 行列減算演算子 (行列－スカラーの場合)
		//! @param[in] right 演算子の右側
		//! @return 結果
		constexpr auto operator-(const TT& right) const{
			Matrix ret;
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) ret.Data[i][j] = Data[i][j] - right;
			}
			return ret;
		}
		
		//! @brief 行列乗算演算子 (行列同士の乗算の場合)
		//! @param[in] right 演算子の右側
		//! @return 結果
		template <size_t Nright, size_t Mright, typename Tright>
		constexpr auto operator*(const Matrix<Nright,Mright,Tright>& right) const{
			static_assert(N == right.M, "Matrix Size Error");	// 行列のサイズチェック
			Matrix<right.N,M,TT> ret;
			for(size_t k = 0; k < right.N; ++k){
				for(size_t i = 0; i < N; ++i){
					for(size_t j = 0; j < M; ++j) ret.Data[k][j] += Data[i][j]*right.Data[k][i];
				}
			}
			return ret;
		}
		
		//! @brief 行列乗算演算子 (行列＊スカラーの場合)
		//! @param[in] right 演算子の右側
		//! @return 結果
		constexpr auto operator*(const TT& right) const{
			Matrix ret;
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) ret.Data[i][j] = Data[i][j]*right;
			}
			return ret;
		}
		
		//! @brief 行列スカラー除算演算子 (行列／スカラーの場合)
		//! @param[in] right 演算子の右側
		//! @return 結果
		constexpr auto operator/(const TT& right) const{
			Matrix ret;
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) ret.Data[i][j] = Data[i][j]/right;
			}
			return ret;
		}

		//! @brief 行列加算代入演算子 (行列＝行列＋行列の場合)
		//! @param[in] right 演算子の右側
		//! @return 結果
		constexpr auto& operator+=(const Matrix& right){
			static_assert(this->N == right.N, "Matrix Size Error");	// 行列のサイズチェック
			static_assert(this->M == right.M, "Matrix Size Error");	// 行列のサイズチェック
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) this->Data[i][j] = this->Data[i][j] + right.Data[i][j];
			}
			return (*this);
		}
		
		//! @brief 行列加算代入演算子 (行列＝行列＋スカラーの場合)
		//! @param[in] right 演算子の右側
		//! @return 結果
		constexpr auto& operator+=(const TT& right){
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) this->Data[i][j] = this->Data[i][j] + right;
			}
			return (*this);
		}
		
		//! @brief 行列減算代入演算子 (行列＝行列－行列の場合)
		//! @param[in] right 演算子の右側
		//! @return 結果
		constexpr auto& operator-=(const Matrix& right){
			static_assert(this->N == right.N, "Matrix Size Error");	// 行列のサイズチェック
			static_assert(this->M == right.M, "Matrix Size Error");	// 行列のサイズチェック
			for(size_t j = 0; j < M; ++j){
				for(size_t i = 0; i < N; ++i) this->Data[i][j] = this->Data[i][j] - right.Data[i][j];
			}
			return (*this);
		}
		
		//! @brief 行列減算代入演算子 (行列＝行列－スカラーの場合)
		//! @param[in] right 演算子の右側
		//! @return 結果
		constexpr auto& operator-=(const TT& right){
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) this->Data[i][j] = this->Data[i][j] - right;
			}
			return (*this);
		}
		
		//! @brief 行列べき乗演算子
		//! @param[in] right 演算子の右側
		//! @return 結果
		constexpr auto operator^(const size_t& right) const{
			Matrix A;
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) A.Data[i][j] = Data[i][j];
			}
			Matrix ret = A;
			for(size_t k = 1; k < right; ++k) ret = ret*A;
			return ret;
		}
		
		//! @brief 行列アダマール積演算子 (行列の要素ごとの乗算)
		//! @param[in] right 演算子の右側
		//! @return 結果
		constexpr auto operator&(const Matrix& right) const{
			static_assert(N == right.N, "Matrix Size Error");	// 行列のサイズチェック
			static_assert(M == right.M, "Matrix Size Error");	// 行列のサイズチェック
			Matrix ret;
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) ret.Data[i][j] = Data[i][j]*right.Data[i][j];
			}
			return ret;
		}
		
		//! @brief 行列アダマール除算演算子 (行列の要素ごとの除算)
		//! @param[in] right 演算子の右側
		//! @return 結果
		constexpr auto operator%(const Matrix& right) const{
			static_assert(N == right.N, "Matrix Size Error");	// 行列のサイズチェック
			static_assert(M == right.M, "Matrix Size Error");	// 行列のサイズチェック
			Matrix ret;
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) ret.Data[i][j] = Data[i][j]/right.Data[i][j];
			}
			return ret;
		}
		
		//! @brief 縦ベクトル添字演算子（縦ベクトルのm番目の要素の値を返す。GetElement(1,m)と同じ意味。ただしサイズチェックは行わない。）
		//! @param[in]	m	要素番号(「1」始まり)
		//! @return	要素の値
		constexpr TT operator[](size_t m) const{
			arcs_assert(m != 0);		// 「0」始まり防止チェック
			return Data[0][m-1];
		}
		
		//! @brief 縦ベクトル添字演算子（縦ベクトルのm番目の要素に値を設定する。SetElement(1,n,val)と同じ意味。ただしサイズチェックは行わない。）
		//! @param[in]	m	要素番号(「1」始まり)
		//! @return	設定後の縦ベクトル
		constexpr TT& operator[](size_t m){
			arcs_assert(m != 0);		// 「0」始まり防止チェック
			return Data[0][m-1];
		}
		
		//! @brief 行列括弧演算子（行列の(m,n)要素の値を返す。GetElem(m,n)と同じ意味。ただしサイズチェックは行わない。）
		//! @param[in]	m	m行目
		//! @param[in]	n	n列目
		//! @return	要素の値
		constexpr TT operator()(size_t m, size_t n) const{
			return Data[n-1][m-1];
		}
		
		//! @brief 行列括弧演算子（行列の(m,n)要素に値を設定する。SetElem(m,n,val)と同じ意味。ただしサイズチェックは行わない。）
		//! @param[in]	m	m行目
		//! @param[in]	n	n列目
		//! @return	設定後の行列
		constexpr TT& operator()(size_t m, size_t n){
			return Data[n-1][m-1];
		}
		
		//! @brief 行列加算演算子 (スカラー＋行列の場合)
		//! @param[in] left		左側のスカラー値
		//! @param[in] right	右側の行列
		constexpr friend Matrix operator+(const TT& left, const Matrix& right){
			Matrix ret;
			for(size_t i = 0; i < right.N; ++i){
				for(size_t j = 0; j < right.M; ++j) ret.Data[i][j] = left + right.Data[i][j];
			}
			return ret;
		}
		
		//! @brief 行列減算演算子 (スカラー－行列の場合)
		//! @param[in] left		左側のスカラー値
		//! @param[in] right	右側の行列
		constexpr friend Matrix operator-(const TT& left, const Matrix& right){
			Matrix<NN,MM,TT> ret;
			for(size_t i = 0; i < right.N; ++i){
				for(size_t j = 0; j < right.M; ++j) ret.Data[i][j] = left - right.Data[i][j];
			}
			return ret;
		}
		
		//! @brief 行列乗算演算子 (スカラー＊行列の場合)
		//! @param[in] left		左側のスカラー値
		//! @param[in] right	右側の行列
		constexpr friend Matrix operator*(const TT& left, const Matrix& right){
			Matrix<NN,MM,TT> ret;
			for(size_t i = 0; i < right.N; ++i){
				for(size_t j = 0; j < right.M; ++j) ret.Data[i][j] = right.Data[i][j]*left;
			}
			return ret;
		}
		
		//! @brief 行列要素に値を設定する関数
		//! @param[in] u1...u2 要素の値
		template<typename T1, typename... T2>		// 可変長引数テンプレート
		constexpr void Set(const T1& u1, const T2&... u2){	// 再帰で順番に可変長引数を読み込んでいく
			arcs_assert(Nindex < N);				// 行カウンタが行の長さ以内かチェック
			arcs_assert(Mindex < M);				// 列カウンタが列の長さ以内かチェック
			Data[Nindex][Mindex] = (TT)u1;	// 行列の要素を埋める
			Nindex++;			// 行カウンタをインクリメント
			if(Nindex == N){	// 行カウンタが最後まで行き着いたら，
				Nindex = 0;		// 行カウンタを零に戻して，
				Mindex++;		// その代わりに，列カウンタをインクリメント
			}
			Set(u2...);			// 自分自身を呼び出す(再帰)
		}
		constexpr void Set(){
			// 再帰の最後に呼ばれる関数
			Nindex = 0;	// すべての作業が終わったので，
			Mindex = 0;	// 行，列カウンタを零に戻しておく
		}
		
		//! @brief 行列要素から値を読み込む関数
		//! @param[in] u1...u2 要素の値
		template<typename T1, typename... T2>		// 可変長引数テンプレート
		constexpr void Get(T1& u1, T2&... u2){	// 再帰で順番に可変長引数を読み込んでいく
			arcs_assert(Nindex < N);	// 行カウンタが行の長さ以内かチェック
			arcs_assert(Mindex < M);	// 列カウンタが列の長さ以内かチェック
			u1 = Data[Nindex][Mindex];	// 行列の要素から読み込み
			Nindex++;			// 行カウンタをインクリメント
			if(Nindex == N){	// 行カウンタが最後まで行き着いたら，
				Nindex = 0;		// 行カウンタを零に戻して，
				Mindex++;		// その代わりに，列カウンタをインクリメント
			}
			Get(u2...);			// 自分自身を呼び出す(再帰)
		}
		constexpr void Get(){
			// 再帰の最後に呼ばれる関数
			Nindex = 0;	// すべての作業が終わったので，
			Mindex = 0;	// 行，列カウンタを零に戻しておく
		}
		
		//! @brief 1次元std::array配列を縦ベクトルとして読み込む関数
		//! @param[in]	Array	std::array配列(横1×縦MM)
		constexpr void LoadArray(const std::array<TT, MM>& Array){
			arcs_assert(N == 1);
			arcs_assert(M == MM);
			for(size_t j = 0; j < M; ++j) Data[0][j] = Array[j];
		}
		
		//! @brief 縦ベクトルを1次元std::array配列に書き込む関数
		//! @param[out]	Array	std::array配列(横1×縦MM)
		constexpr void StoreArray(std::array<TT, MM>& Array) const{
			arcs_assert(N == 1);
			arcs_assert(M == MM);
			for(size_t j = 0; j < M; ++j) Array[j] = Data[0][j];
		}
		
		//! @brief 2次元std::array配列を読み込む関数
		//! @param[in]	Array	std::array配列(横NN×縦MM)
		constexpr void LoadArray(const std::array<std::array<TT, NN>, MM>& Array){
			arcs_assert(N == NN);
			arcs_assert(M == MM);
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) Data[i][j] = Array[j][i];
			}
		}
		
		//! @brief 短い縦ベクトルを読み込み，且つ残りはゼロで埋める関数
		//! @param[in]	v	短い縦ベクトル
		template<size_t VM>
		constexpr void LoadShortVector(const Matrix<1,VM,TT>& v){
			static_assert(VM < MM);	// 読み込む縦ベクトルの方が短いかチェック
			for(size_t i = 0; i < VM; ++i){
				Data[0][i] = v.Data[0][i];
			}
		}

		//! @brief 指定した要素番号に値を設定する関数 (並び順逆版)
		//! @param[in]	m	m行目
		//! @param[in]	n	n列目
		//! @param[in]	val	valを代入
		constexpr void SetElem(size_t m, size_t n, TT val){
			arcs_assert(1<=n && n<=N && 1<=m && m<=M);	// サイズチェック
			Data[n-1][m-1] = val;
		}
		
		//! @brief 指定した要素番号の値を返す関数 (並び順逆版)
		//! @param[in]	m	m行目
		//! @param[in]	n	n列目
		//! @return	要素の値
		constexpr TT GetElem(size_t m, size_t n) const {
			arcs_assert(1<=n && n<=N && 1<=m && m<=M);	// サイズチェック
			return Data[n-1][m-1];
		}
		
		//! @brief 指定した要素番号に値を設定する関数
		//! @param[in]	m	m行目
		//! @param[in]	n	n列目
		//! @param[in]	val	valを代入
		constexpr void SetElement(size_t n, size_t m, TT val){
			arcs_assert(1<=n && n<=N && 1<=m && m<=M);	// サイズチェック
			Data[n-1][m-1] = val;
		}
		
		//! @brief 指定した要素番号の値を返す関数
		//! @param[in]	m	m行目
		//! @param[in]	n	n列目
		//! @return	要素の値
		constexpr TT GetElement(size_t n, size_t m) const {
			arcs_assert(1<=n && n<=N && 1<=m && m<=M);	// サイズチェック
			return Data[n-1][m-1];
		}
		
		//! @brief 指定した先頭位置から縦ベクトルを抜き出して返す関数
		//! @param[in]	m	m行目
		//! @param[in]	n	n列目
		//! @return	縦ベクトル
		template<size_t VM>
		constexpr Matrix<1,VM,TT> GetVerticalVec(size_t n, size_t m){
			arcs_assert(1<=n && n<=N && 1<=m && m<=M);	// サイズチェック
			Matrix<1,VM,TT> r;
			for(size_t i = 0; i < VM; ++i){
				r.Data[0][i] = Data[n-1][m-1+i];
			}
			return r;
		}
		
		//! @brief 指定した先頭位置から横ベクトルを抜き出して返す関数
		//! @param[in]	m	m行目
		//! @param[in]	n	n列目
		//! @return	横ベクトル
		template<size_t HN>
		constexpr Matrix<HN,1,TT> GetHorizontalVec(size_t n, size_t m){
			arcs_assert(1<=n && n<=N && 1<=m && m<=M);	// サイズチェック
			Matrix<HN,1,TT> r;
			for(size_t i = 0; i < HN; ++i){
				r.Data[i][0] = Data[n-1+i][m-1];
			}
			return r;
		}
		
		//! @brief 指定した先頭位置に縦ベクトルを埋め込む関数
		//! @param[in]	v	縦ベクトル
		//! @param[in]	m	m行目
		//! @param[in]	n	n列目
		template<size_t VM>
		constexpr void SetVerticalVec(Matrix<1,VM,TT> v, size_t n, size_t m){
			arcs_assert(1<=n && n<=N && 1<=m && m<=M);	// サイズチェック
			for(size_t i = 0; i < VM; ++i){
				Data[n-1][m-1+i] = v.Data[0][i];
			}
		}
		
		//! @brief 指定した先頭位置に横ベクトルを埋め込む関数
		//! @param[in]	h	横ベクトル
		//! @param[in]	m	m行目
		//! @param[in]	n	n列目
		template<size_t HN>
		constexpr void SetHorizontalVec(Matrix<HN,1,TT> h, size_t n, size_t m){
			arcs_assert(1<=n && n<=N && 1<=m && m<=M);	// サイズチェック
			for(size_t i = 0; i < HN; ++i){
				Data[n-1+i][m-1] = h.Data[i][0];
			}
		}
		
		//! @brief すべての要素を指定した値で埋める関数
		//! @param[in] u 埋める値
		constexpr void FillAll(TT u){
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) Data[i][j] = u;
			}
		}
		
		//! @brief すべての要素を指定したゼロで埋める関数
		constexpr void FillAllZero(void){
			FillAll(0);
		}

		//! @brief 行列の幅(列数)を返す関数
		//! @return 行列の幅
		constexpr size_t GetWidthLength(void) const {
			return N;
		}

		//! @brief 行列の高さ(行数)を返す関数
		//! @return 行列の幅
		constexpr size_t GetHeightLength(void) const {
			return M;
		}
		
		//! @brief MatrixクラスがSIMD命令が有効に設定されているかを返す関数
		//! @return SIMD有効 = true, SIMD無効 = false
		constexpr bool isEnabledSIMD(void){
			return false;
		}
		
		//! @brief 行列のサイズの表示 (この関数はマクロを介して呼ばれることを想定している)
		//! @param[in] u		表示する行列
		//! @param[in] varname	変数名
		friend void PrintMatSize_Macro(const Matrix& u, const std::string& varname){
			printf("%s = [ %zu x %zu ]\n", varname.c_str(), u.N, u.M);
		}
		
		//! @brief 行列の要素を表示 (書式指定あり版，この関数はマクロを介して呼ばれることを想定している)
		//! @param[in] u		表示する行列
		//! @param[in] format	表示形式 (%1.3e とか %5.3f とか printfと同じ)
		//! @param[in] varname	変数名
		friend void PrintMatrix_Macro(const Matrix& u, const std::string& format, const std::string& varname){
			printf("%s = \n", varname.c_str());
			for(size_t j = 0; j < u.M; ++j){
				printf("[ ");
				for(size_t i = 0; i < u.N; ++i){
					// データ型によって表示方法を変える
					if constexpr(std::is_same_v<TT, std::complex<double>> || std::is_same_v<TT, std::complex<float>>){
						// 複素数型の場合
						// 実数部の表示
						printf(format.c_str(), u.Data[i][j].real());
						// 虚数部の表示
						if(0.0 <= u.Data[i][j].imag()){
							printf(" + j");
						}else{
							printf(" - j");
						}
						printf( format.c_str(), std::abs(u.Data[i][j].imag()) );
					}else{
						// それ以外の場合
						printf(format.c_str(), u.Data[i][j]);
					}
					printf(" ");
				}
				printf("]\n");
			}
			printf("\n");
		}
		
		//! @brief 行列の要素を表示 (書式指定なし版，この関数はマクロを介して呼ばれることを想定している)
		//! @param[in] u		表示する行列
		//! @param[in] varname	変数名
		friend void PrintMat_Macro(const Matrix& u, const std::string& varname){
			// データ型によって書式指定子を変える
			// double型かfloat型の場合
			if constexpr(std::is_same_v<TT, double> || std::is_same_v<TT, float>){
				PrintMatrix_Macro(u, "% g", varname);
				return;
			}
			
			// int型の場合
			if constexpr(std::is_same_v<TT, int>){
				PrintMatrix_Macro(u, "% d", varname);
				return;
			}
			
			// long型の場合
			if constexpr(std::is_same_v<TT, long>){
				PrintMatrix_Macro(u, "% ld", varname);
				return;
			}
			
			// 複素数型の場合
			if constexpr(std::is_same_v<TT, std::complex<double>> || std::is_same_v<TT, std::complex<float>>){
				PrintMatrix_Macro(u, "% g", varname);
				return;
			}
		}
		
		//! @brief m行n列の零行列を返す関数
		//! @return 零行列
		static constexpr Matrix zeros(void){
			Matrix ret;
			return ret;
		}
		
		//! @brief m行n列の要素がすべて1の行列を返す関数
		//! @return 1行列
		static constexpr Matrix ones(void){
			Matrix ret;
			ret.FillAll(1);
			return ret;
		}
		
		//! @brief n行n列の単位行列を返す関数
		//! @return 単位行列
		static constexpr Matrix<NN,NN,TT> ident(void){
			static_assert(NN == MM, "Matrix Size Error");	// 行列のサイズチェック
			Matrix ret;
			if constexpr(std::is_same_v<TT, std::complex<double>> || std::is_same_v<TT, std::complex<float>>){
				// 複素数型の場合
				for(size_t i = 0; i < N; ++i) ret.SetElem(i+1,i+1, std::complex(1.0,0.0));	// 対角成分を 1 + j0 で埋める
			}else{
				// それ以外の場合
				for(size_t i = 0; i < N; ++i) ret.SetElem(i+1,i+1, 1);	// 対角成分を 1 で埋める
			}
			return ret;
		}
		
		//! @brief n行n列の単位行列を返す関数 (identのエイリアス)
		//! @return 単位行列
		static constexpr Matrix<NN,NN,TT> eye(void){
			return Matrix<NN,NN,TT>::ident();
		}
		
		//! @brief 単調増加の縦ベクトルを返す関数
		//! @return 縦ベクトル
		static constexpr Matrix ramp(void){
			static_assert(NN == 1, "Matrix Size Error");	// 行列のサイズチェック
			Matrix ret;
			for(size_t j = 1; j <= M; ++j) ret[j] = j;	// 単調増加を書き込む
			return ret;
		}
		
		//! @brief 転置行列を返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend Matrix<MM,NN,TT> tp(const Matrix<NN,MM,TT>& U){
			Matrix<U.M,U.N,TT> y;
			for(size_t j = 0; j < U.M; ++j){
				for(size_t i = 0; i < U.N; ++i) y.Data[j][i] = U.Data[i][j];
			}
			return y;
		}
		
		//! @brief 行列のトレースを返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend TT tr(const Matrix& U){
			static_assert(U.N == U.M, "Matrix Size Error");	// 正方行列であることのチェック
			TT y = 0;
			for(size_t i = 0; i < U.N; ++i) y += (TT)U.Data[i][i];
			return y;
		}
		
		//! @brief 行列の対角要素の総積を返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend TT prod(const Matrix& U){
			static_assert(U.N == U.M, "Matrix Size Error");	// 正方行列であることのチェック
			TT y = 1;
			for(size_t i = 0; i < U.N; ++i) y *= U.Data[i][i];
			return y;
		}
		
		//! @brief 行列の対角要素を縦ベクトルで返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend Matrix<1,std::min(NN,MM),TT> diag(const Matrix& U){
			constexpr size_t K = std::min(NN,MM);
			Matrix<1,K,TT> y;
			for(size_t i = 0; i < K; ++i) y.Data[0][i] = U.Data[i][i];
			return y;
		}
		
		//! @brief 行方向(横方向)へ加算して縦ベクトルを出力する関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend Matrix<1,MM,TT> sumrow(const Matrix& U){
			Matrix<1,U.M,TT> y;
			for(size_t i = 0; i < U.N; ++i){
				for(size_t j = 0; j < U.M; ++j) y.Data[0][j] += U.Data[i][j];
			}
			return y;
		}
		
		//! @brief 列方向(縦方向)へ加算して横ベクトルを出力する関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend Matrix<NN,1,TT> sumcolumn(const Matrix& U){
			Matrix<U.N,1,TT> y;
			for(size_t i = 0; i < U.N; ++i){
				for(size_t j = 0; j < U.M; ++j) y.Data[i][0] += U.Data[i][j];
			}
			return y;
		}
		
		//! @brief 行列の全要素を加算して出力する関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend double sumall(const Matrix& U){
			const Matrix<1,1> y = sumrow(sumcolumn(U));
			return y[1];
		}
		
		//! @brief ベクトル要素の最大値を返す関数
		//! @param[in]	u	入力ベクトル
		//! @return	結果
		constexpr friend TT max(const Matrix& u){
			static_assert( (u.N == 1)||(u.M == 1), "Input is NOT vector." );	// ベクトルかチェック
			size_t k = 0;
			TT y = 0;
			if(u.M == 1){
				// 行ベクトルの場合
				for(size_t i = 1; i < u.N; ++i)if(u.Data[k][0] < u.Data[i][0]) k = i;
				y = u.Data[k][0];
			}
			if(u.N == 1){
				// 列ベクトルの場合
				for(size_t i = 1; i < u.M; ++i)if(u.Data[0][k] < u.Data[0][i]) k = i;
				y = u.Data[0][k];
			}
			return y;
		}
		
		//! @brief ベクトル要素の絶対値の最大値を返す関数
		//! @param[in]	u	入力ベクトル
		//! @return	結果
		constexpr friend TT absmax(const Matrix& u){
			static_assert( (u.N == 1)||(u.M == 1), "Input is NOT vector." );	// ベクトルかチェック
			size_t k = 0;
			TT y = 0;
			if(u.M == 1){
				// 行ベクトルの場合
				for(size_t i = 1; i < u.N; ++i)if(fabs(u.Data[k][0]) < fabs(u.Data[i][0])) k = i;
				y = u.Data[k][0];
			}
			if(u.N == 1){
				// 列ベクトルの場合
				for(size_t i = 1; i < u.M; ++i)if(fabs(u.Data[0][k]) < fabs(u.Data[0][i])) k = i;
				y = u.Data[0][k];
			}
			return y;
		}
		
		//! @brief ベクトル要素の最大値の要素番号を返す関数
		//! @param[in]	u	入力ベクトル
		//! @return	結果
		constexpr friend size_t maxidx(const Matrix& u){
			static_assert( (u.N == 1)||(u.M == 1), "Input is NOT vector." );	// ベクトルかチェック
			size_t k = 0;
			if(u.M == 1){
				// 行ベクトルの場合
				for(size_t i = 1; i < u.N; ++i)if(u.Data[k][0] < u.Data[i][0]) k = i;
			}
			if(u.N == 1){
				// 列ベクトルの場合
				for(size_t i = 1; i < u.M; ++i)if(u.Data[0][k] < u.Data[0][i]) k = i;
			}
			return k + 1;
		}
		
		//! @brief ベクトル要素の絶対値の最大値の要素番号を返す関数
		//! @param[in]	u	入力ベクトル
		//! @return	結果
		constexpr friend size_t absmaxidx(const Matrix& u){
			static_assert( (u.N == 1)||(u.M == 1), "Input is NOT vector." );	// ベクトルかチェック
			size_t k = 0;
			if(u.M == 1){
				// 行ベクトルの場合
				for(size_t i = 1; i < u.N; ++i)if(fabs(u.Data[k][0]) < fabs(u.Data[i][0])) k = i;
			}
			if(u.N == 1){
				// 列ベクトルの場合
				for(size_t i = 1; i < u.M; ++i)if(fabs(u.Data[0][k]) < fabs(u.Data[0][i])) k = i;
			}
			return k + 1;
		}
		
		//! @brief 行列の非ゼロ要素数を返す関数
		//! @param[in]	U	入力行列
		//! @return 結果
		constexpr friend size_t nonzeroele(const Matrix& U){
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
		constexpr friend size_t rank(const Matrix& A){
			static_assert(NN == MM, "Matrix Size Error");	// 正方行列のみ対応
			Matrix<NN,NN,TT> U, S, V;
			SVD(A, U, S, V);			// まず特異値分解して，
			return nonzeroele(diag(S));	// S行列の対角要素の非ゼロをカウントするとそれがランク
		}
		
		//! @brief 指定した行から横ベクトルとして抽出する関数
		//! @param[in]	U	入力行列
		//! @param[in]	m	抽出したい行
		//! @return	結果
		constexpr friend Matrix<NN,1,TT> getrow(const Matrix& U, size_t m){
			arcs_assert(0 < m && m <= U.M);	// 指定した行が行列の高さ以内かチェック
			Matrix<U.N,1,TT> y;
			for(size_t i=0;i<U.N;i++) y.Data[i][0] = U.Data[i][m-1];
			return y;
		}
		
		//! @brief 指定した行を横ベクトルで上書きする関数
		//! @param[in,out]	U	入出力行列
		//! @param[in]	v	横ベクトル
		//! @param[in]	m	上書きしたい行
		constexpr friend void setrow(Matrix& U, const Matrix<NN,1,TT>& v, size_t m){
			static_assert(U.N == v.N, "Matrix Size Error");	// 行列の幅と横ベクトルの長さが等しいかチェック
			static_assert(v.M == 1, "Input is NOT vector.");// 行ベクトルかチェック
			arcs_assert(0 < m && m <= U.M);	// 指定した行が行列の高さ以内かチェック
			for(size_t i=0;i<U.N;i++) U.Data[i][m-1] = v.Data[i][0];
		}
		
		//! @brief 指定した行を横ベクトル(std::array)で上書きする関数
		//! @param[in,out]	U	入出力行列
		//! @param[in]	v	横ベクトル
		//! @param[in]	m	上書きしたい行
		constexpr friend void setrow(Matrix& U, const std::array<TT, NN>& v, size_t m){
			static_assert(U.N == NN, "Matrix Size Error");	// 行列の幅と横ベクトルの長さが等しいかチェック
			arcs_assert(0 < m && m <= U.M);	// 指定した行が行列の高さ以内かチェック
			for(size_t i=0;i<U.N;i++) U.Data[i][m-1] = v.at(i);
		}
		
		//! @brief 指定した行と行を入れ替える関数
		//! @param[in,out]	U	入出力行列
		//! @param[in]	m1	指定行1
		//! @param[in]	m2	指定行2
		constexpr friend void swaprow(Matrix& U, size_t m1, size_t m2){
			arcs_assert(m1 <= U.M);	// 行が高さ以内かチェック
			arcs_assert(m2 <= U.M);	// 行が高さ以内かチェック
			Matrix<U.N,1,TT> p, q;	// バッファ用ベクトル
			p = getrow(U,m1);	// m1行目を抽出
			q = getrow(U,m2);	// m2行目を抽出
			setrow(U, p, m2);	// m2行目に書き込み
			setrow(U, q, m1);	// m1行目に書き込み
		}
		
		//! @brief m行目のn1列目からm2列目を数値aで埋める関数  (n1 <= n2 であること)
		//! @param[in,out]	U	入出力行列
		//! @param[in]	a	埋める値
		//! @param[in]	m	指定行
		//! @param[in]	n1	開始列
		//! @param[in]	n2	終了列
		constexpr friend void fillrow(Matrix& U, TT a, size_t m, size_t n1, size_t n2){
			arcs_assert(0 < m && m <= U.M);		// 列が幅以内かチェック
			arcs_assert(0 < n1 && n1 <= U.N);	// 列が高さ以内かチェック
			arcs_assert(0 < n2 && n2 <= U.N);	// 列が高さ以内かチェック
			arcs_assert(n1 <= n2);				// 開始列と終了列が入れ替わらないかチェック
			for(size_t i=n1;i<=n2;i++) U.Data[i-1][m-1] = a;
		}
		
		//! @brief 指定した列から縦ベクトルとして抽出する関数
		//! @param[in]	U	入力行列
		//! @param[in]	n	指定列
		//! @return	結果
		constexpr friend Matrix<1,MM,TT> getcolumn(const Matrix& U, size_t n){
			arcs_assert(n <= U.N);	// 指定した列が行列の幅以内かチェック
			Matrix<1,U.M,TT> y;
			for(size_t i=0;i<U.M;i++) y.Data[0][i] = U.Data[n-1][i];
			return y;
		}
		
		//! @brief 指定した列を縦ベクトルで上書きする関数
		//! @param[in,out]	U	入出力行列
		//! @param[in]	v	縦ベクトル
		//! @param[in]	n	上書きしたい列
		constexpr friend void setcolumn(Matrix& U, const Matrix<1,MM,TT>& v, size_t n){
			static_assert(U.M == v.M, "Matrix Size Error");		// 行列の高さ列ベクトルの長さが等しいかチェック
			static_assert(v.N == 1, "Input is NOT vector.");	// 列ベクトルかチェック
			arcs_assert(0 < n && n <= U.N);	// 指定した列が行列の幅以内かチェック
			for(size_t i=0;i<U.M;i++) U.Data[n-1][i] = v.Data[0][i];
		}
		
		//! @brief 指定した列を縦ベクトル(std::array)で上書きする関数
		//! @param[in,out]	U	入出力行列
		//! @param[in]	v	縦ベクトル
		//! @param[in]	n	上書きしたい列
		constexpr friend void setcolumn(Matrix& U, const std::array<TT,MM>& v, size_t n){
			static_assert(U.M == MM, "Matrix Size Error");		// 行列の高さ列ベクトルの長さが等しいかチェック
			arcs_assert(0 < n && n <= U.N);	// 指定した列が行列の幅以内かチェック
			for(size_t i=0;i<U.M;i++) U.Data[n-1][i] = v.at(i);
		}
		
		//! @brief 指定した列と列を入れ替える関数
		//! @param[in,out]	U	入出力行列
		//! @param[in]	n1	指定列1
		//! @param[in]	n2	指定列2
		constexpr friend void swapcolumn(Matrix& U, size_t n1, size_t n2){
			arcs_assert(n1 <= U.N);	// 列が幅以内かチェック
			arcs_assert(n2 <= U.N);	// 列が幅以内かチェック
			Matrix<1,U.M,TT> p, q;	// バッファ用ベクトル
			p = getcolumn(U,n1);	// n1列目を抽出
			q = getcolumn(U,n2);	// n2列目を抽出
			setcolumn(U,p,n2);		// n2列目に書き込み
			setcolumn(U,q,n1);		// n1列目に書き込み
		}
		
		//! @brief n列目のm1行目からm2行目を数値aで埋める関数  (m1 <= m2 であること)
		//! @param[in,out]	U	入出力行列
		//! @param[in]	a	埋める値
		//! @param[in]	n	指定列
		//! @param[in]	m1	開始行
		//! @param[in]	m2	終了行
		constexpr friend void fillcolumn(Matrix& U, TT a, size_t n, size_t m1, size_t m2){
			arcs_assert(n <= U.N);	// 列が幅以内かチェック
			arcs_assert(m1 <= U.M);	// 行が高さ以内かチェック
			arcs_assert(m2 <= U.M);	// 行が高さ以内かチェック
			arcs_assert(m1 <= m2);	// 開始行と終了行が入れ替わらないかチェック
			for(size_t i=m1;i<=m2;i++) U.Data[n-1][i-1] = a;
		}
		
		//! @brief 指定した列を縦ベクトルで指定位置に上書きする関数
		//! @tparam	VM	縦ベクトルの長さ
		//! @param[in,out]	U	入出力行列
		//! @param[in]	v	縦ベクトル
		//! @param[in]	n	上書きしたい列
		//! @param[in]	m	縦方向の書き込む先頭位置
		template <size_t VM>
		constexpr friend void setvvector(Matrix<NN,MM,TT>& U, const Matrix<1,VM,TT>& v, size_t n, size_t m){
			arcs_assert(0 < n && n <= NN);	// 指定した列が行列の幅以内かチェック
			arcs_assert(VM + m - 1 <= MM);	// 縦ベクトルの下がハミ出ないかチェック
			for(size_t i = 0; i < VM; ++i) U.Data[n-1][m-1+i] = v.Data[0][i];
		}
		
		//! @brief 指定した列から縦ベクトルを指定位置から抽出する関数
		//! @tparam	VM	縦ベクトルの長さ
		//! @param[in]	U	入力行列
		//! @param[in]	n	読み込みたい列
		//! @param[in]	m	縦方向の読み込む先頭位置
		//! @param[out]	v	抽出した縦ベクトル
		template <size_t VM>
		constexpr friend void getvvector(const Matrix<NN,MM,TT>& U, size_t n, size_t m, Matrix<1,VM,TT>& v){
			arcs_assert(0 < n && n <= NN);	// 指定した列が行列の幅以内かチェック
			arcs_assert(VM + m - 1 <= MM);	// 縦ベクトルの下がハミ出ないかチェック
			for(size_t i = 0; i < VM; ++i) v.Data[0][i] = U.Data[n-1][m-1+i];
		}
		
		//! @brief 指定した行から横ベクトルを指定位置から抽出する関数
		//! @tparam	VM	縦ベクトルの長さ
		//! @param[in]	U	入力行列
		//! @param[in]	m	読み込みたい行
		//! @param[in]	n	横方向の読み込む先頭位置
		//! @param[out]	h	抽出した横ベクトル
		template <size_t VN>
		constexpr friend void gethvector(const Matrix<NN,MM,TT>& U, size_t m, size_t n, Matrix<VN,1,TT>& h){
			arcs_assert(0 < m && m <= MM);	// 指定した行が行列の高さ以内かチェック
			arcs_assert(VN + n - 1 <= NN);	// 横ベクトルの右がハミ出ないかチェック
			for(size_t i = 0; i < VN; ++i) h.Data[0][i] = U.Data[n-1+i][m-1];
		}
		
		//! @brief 行列から指定位置の小行列を抽出する関数
		//! @tparam	SN	小行列の幅
		//! @tparam	SM	小行列の高さ
		//! @param[in]	U	入力行列
		//! @param[in]	n	抽出する横位置（小行列の左）
		//! @param[in]	m	抽出する縦位置（小行列の上）
		//! @param[out]	Y	抽出した小行列
		template <size_t SN, size_t SM>
		constexpr friend void getsubmatrix(const Matrix<NN,MM,TT>& U, size_t n, size_t m, Matrix<SN,SM,TT>& Y){
			static_assert(SN <= NN);		// 小行列の方が幅が小さいかチェック
			static_assert(SM <= MM);		// 小行列の方が高さが小さいかチェック
			arcs_assert(SN + n - 1 <= NN);	// 右側がハミ出ないかチェック
			arcs_assert(SM + m - 1 <= MM);	// 下側がハミ出ないかチェック
			Matrix<1,SM,TT> v;
			for(size_t i = 1; i <= SN; ++i){
				getvvector(U, i + n - 1, m, v);
				setcolumn(Y, v, i);
			}
		}
		
		//! @brief 小行列を行列の指定位置に上書きする関数
		//! @tparam	SN	小行列の幅
		//! @tparam	SM	小行列の高さ
		//! @param[in,out]	U	入出力行列
		//! @param[in]	n	上書きしたい横位置（小行列の左）
		//! @param[in]	m	上書きしたい縦位置（小行列の上）
		//! @param[in]	A	書き込む小行列
		template <size_t SN, size_t SM>
		constexpr friend void setsubmatrix(Matrix<NN,MM,TT>& U, size_t n, size_t m, const Matrix<SN,SM,TT>& A){
			static_assert(SN <= NN);		// 小行列の方が幅が小さいかチェック
			static_assert(SM <= MM);		// 小行列の方が高さが小さいかチェック
			arcs_assert(SN + n - 1 <= NN);	// 右側がハミ出ないかチェック
			arcs_assert(SM + m - 1 <= MM);	// 下側がハミ出ないかチェック
			Matrix<1,SM,TT> v;
			for(size_t i = 1; i <= SN; ++i){
				v = getcolumn(A, i);
				setvvector(U, v, i + n - 1, m);
			}
		}
		
		//! @brief 並び替え記憶列ベクトルvの行番号に従って，入力行列Uの行を並び替える関数
		//! @param[in]	U	入力行列
		//! @param[in]	v	並び替え記憶列ベクトル
		//! @return	結果
		constexpr friend Matrix orderrow(const Matrix& U, const Matrix<1,MM,int>& v){
			static_assert(U.M == v.M, "Matrix and Vector Size Error");	// 行列の高さと列ベクトルの長さが同じかチェック
			Matrix Y;
			for(size_t i=0;i<U.M;i++){
				setrow(Y,  getrow(U, (size_t)v.Data[0][i]), i+1);
			}
			return Y;
		}
		
		//! @brief 並び替え記憶列ベクトルvが昇順になるように，入力行列Uの行を並び替えて元に戻す関数
		//! @param[in]	U	入力行列
		//! @param[in]	v	並び替え記憶列ベクトル
		//! @return	結果
		constexpr friend Matrix reorderrow(const Matrix& U, const Matrix<1,MM,int>& v){
			static_assert(U.M == v.M, "Matrix and Vector Size Error");	// 行列の高さと列ベクトルの長さが同じかチェック
			Matrix Y = U;
			Matrix<1,MM,int> p = v;
			for(size_t i=1;i<=U.M;i++){
				swaprow(Y,i,p.Data[0][i-1]);
				swaprow(p,i,p.Data[0][i-1]);
			}
			return Y;
		}
		
		//! @brief 行列の各要素を上に1行分シフトする関数(最下段の行はゼロになる)
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend Matrix shiftup(const Matrix& U){
			Matrix Y;
			for(size_t i = 2; i <= U.M; ++i){
				setrow(Y, getrow(U, i), i - 1);	// 横ベクトルを抽出して１つ上の行に書き込み
			}
			return Y;
		}
		
		//! @brief 行列の各要素を上にa行分シフトする関数(最下段の行はゼロになる)
		//! @param[in]	U	入力行列
		//! @param[in]	a	シフトする行数
		//! @return	結果
		constexpr friend Matrix shiftup(const Matrix& U, const size_t a){
			Matrix Y = U;
			for(size_t i = 0; i < a; ++i){
				Y = shiftup(Y);
			}
			return Y;
		}
		
		//! @brief 行列の各要素を下に1行分シフトする関数(最上段の行はゼロになる)
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend Matrix shiftdown(const Matrix& U){
			Matrix Y;
			for(size_t i = 2; i <= U.M; ++i){
				setrow(Y, getrow(U, i - 1), i);	// 横ベクトルを抽出して１つ下の行に書き込み
			}
			return Y;
		}
		
		//! @brief 行列の各要素を下にa行分シフトする関数(最上段の行はゼロになる)
		//! @param[in]	U	入力行列
		//! @param[in]	a	シフトする行数
		//! @return	結果
		constexpr friend Matrix shiftdown(const Matrix& U, const size_t a){
			Matrix Y = U;
			for(size_t i = 0; i < a; ++i){
				Y = shiftdown(Y);
			}
			return Y;
		}
		
		//! @brief 行列の各要素を右に1列分シフトする関数(最左段の列はゼロになる)
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend Matrix shiftright(const Matrix& U){
			Matrix Y;
			for(size_t i = 2; i <= U.M; ++i){
				setcolumn(Y, getcolumn(U, i - 1), i);	// 縦ベクトルを抽出して１つ右の行に書き込み
			}
			return Y;
		}
		
		//! @brief 行列の各要素を右にa列分シフトする関数(最左段の列はゼロになる)
		//! @param[in]	U	入力行列
		//! @param[in]	a	シフトする行数
		//! @return	結果
		constexpr friend Matrix shiftright(const Matrix& U, const size_t a){
			Matrix Y = U;
			for(size_t i = 0; i < a; ++i){
				Y = shiftright(Y);
			}
			return Y;
		}
		
		//! @brief 行列の各要素を左に1列分シフトする関数(最右段の列はゼロになる)
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend Matrix shiftleft(const Matrix& U){
			Matrix Y;
			for(size_t i = 2; i <= U.M; ++i){
				setcolumn(Y, getcolumn(U, i), i - 1);	// 縦ベクトルを抽出して１つ左の行に書き込み
			}
			return Y;
		}
		
		//! @brief 行列の各要素を左にa列分シフトする関数(最右段の列はゼロになる)
		//! @param[in]	U	入力行列
		//! @param[in]	a	シフトする行数
		//! @return	結果
		constexpr friend Matrix shiftleft(const Matrix& U, const size_t a){
			Matrix Y = U;
			for(size_t i = 0; i < a; ++i){
				Y = shiftleft(Y);
			}
			return Y;
		}
		
		//! @brief 行列のk番目より上の上三角部分を返す関数(下三角はゼロになる)
		//! @param[in]	U	入力行列
		//! @param[in]	k	k番目
		//! @return	結果
		constexpr friend Matrix gettriup(const Matrix& U, const size_t k){
			Matrix Y;
			for(size_t i = 0; i < MM; ++i){
				for(size_t j = i + k; j < NN; ++j){
					Y.Data[j][i] = U.Data[j][i];
				}
			}
			return Y;
		}
		
		//! @brief 行列の上三角部分を返す関数(下三角はゼロになる)
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend Matrix gettriup(const Matrix& U){
			return gettriup(U, 0);
		}
		
		//! @brief 行列の無限大ノルムを返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend TT infnorm(const Matrix& U){
			return max(sumcolumn(abse(tp(U))));
		}
		
		//! @brief ベクトルのユークリッドノルムを返す関数
		//! @param[in]	v	入力ベクトル
		//! @return	結果
		constexpr friend TT euclidnorm(const Matrix<NN,MM,TT>& v){
			Matrix<1,1,TT> ret;
			if constexpr(std::is_same_v<TT,std::complex<double>>){
				// 複素数型の場合
				if constexpr(NN == 1){
					// 縦ベクトルの場合
					ret = Htp(v)*v;
				}else if constexpr(MM == 1){
					// 横ベクトルの場合
					ret = v*Htp(v);
				}
			}else{
				// 実数型の場合
				if constexpr(NN == 1){
					// 縦ベクトルの場合
					ret = tp(v)*v;
				}else if constexpr(MM == 1){
					// 横ベクトルの場合
					ret = v*tp(v);
				}else{
					// 行列の場合
					ret = sumcolumn(sumrow(v & v));
				}
			}
			return std::sqrt(ret[1]);
		}
		
		//! @brief LU分解
		//! @param[in]	A	入力行列
		//! @param[out]	L	下三角行列
		//! @param[out]	U	上三角行列
		//! @param[out]	v	並べ替え記憶列ベクトル(int型)
		//! @return	並べ替え回数が奇数か偶数かを返す
		constexpr friend enum LUperm LU(const Matrix& A, Matrix& L, Matrix& U, Matrix<1,MM,int>& v){
			static_assert(A.N == A.M, "Matrix Size Error");	// 正方行列かチェック
			Matrix X = A;
			size_t perm_count = 0;
			double max_buff = 0;
			Matrix::LUperm ret = Matrix::EVEN;
			
			// 並べ替え記憶列ベクトルの準備
			for(size_t i = 0; i < X.N; ++i) v.Data[0][i] = i + 1;	// 並べ替え記憶列ベクトルに1から昇順の番号を書き込む
			
			size_t k = 0;
			for(size_t j = 0; j < X.N-1; ++j){
				// j列目の中での行要素の最大値を探す
				k = j;	// 対角要素番号で初期化
				max_buff = std::abs(X.Data[j][j]);	// 対角要素で初期化
				for(size_t i = j + 1; i < X.M; ++i){
					if(max_buff < std::abs(X.Data[j][i])){
						k = i;
						max_buff = std::abs(X.Data[j][i]);
					}
				}
				// 行入れ替えが必要なときは，
				if(k != j){
					swaprow(v,j+1,k+1);				// 並べ替え記憶列ベクトルのj行目とk行目を並べ替え
					swaprow(X,j+1,k+1);				// j行目とk行目を並べ替え
					perm_count++;
				}
				if( fabs(X.Data[j][j]) < X.epsilon )continue;
				{
					// 対角要素が零なら，j列目においてはLU分解が終わっているので，以下の処理はスキップ
					// ここからLU分解
					for(size_t i = j + 1; i < X.M; ++i){
						X.Data[j][i] /= X.Data[j][j];	// 対角要素で除算
						for(size_t l = j + 1; l < X.N; ++l){
							X.Data[l][i] -= X.Data[j][i]*X.Data[l][j];
						}
					}
				}
			}
			// 下三角行列と上三角行列に分離する
			for(size_t j = 0; j < X.N; ++j){
				for(size_t i = j; i < X.M; ++i){
					if(i == j){
						L.Data[j][i] = 1;				// 下三角行列の対角要素はすべて1
					}else{
						L.Data[j][i] = X.Data[j][i];	// 下三角のみコピー
					}
				}
				for(size_t i = 0; i <= j; ++i) U.Data[j][i] = X.Data[j][i];		// 上三角のみコピー
			}
			// 並べ替え回数の判定
			if(perm_count % 2 == 0){
				ret = Matrix::EVEN;	// 奇数のとき
			}else{
				ret = Matrix::ODD;	// 偶数のとき
			}
			return ret;
		}
		
		//! @brief 修正コレスキー分解(LDL^T版)
		//! @param[in]	A	入力行列
		//! @param[out]	L	下三角行列
		//! @param[out]	D	対角行列
		constexpr friend void Cholesky(const Matrix& A, Matrix& L, Matrix& D){
			static_assert(A.N == A.M, "Matrix Size Error");	// Aが正方行列かチェック
			L.Data[0][0] = A.Data[0][0];
			D.Data[0][0] = 1.0/L.Data[0][0];
			
			for(size_t i = 1; i < A.N; ++i){
				for(size_t j = 0; j <= i; ++j){
					TT lld = A.Data[j][i];
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
		constexpr friend void Cholesky(const Matrix& A, Matrix& L){
			Matrix Lp, Dp;
			Cholesky(A, Lp, Dp);// まずコレスキー分解してから，
			L = Lp*sqrte(Dp);	// 対角行列の平方根を取って，下三角行列に掛けたものを出力
		}
		
		//! @brief QR分解
		//! 補足：実数型のときMATLABとはQとRの符号関係が逆の場合があるが正常なQR分解であることは確認済み
		//! 補足：複素数型のときMATLABとは全く違う値が得られるが，正常なQR分解であることは確認済み
		//! @param[in]	A	入力行列
		//! @param[out]	Q	ユニタリ行列 Q行列
		//! @param[out] R	上三角行列 R行列
		constexpr friend void QR(const Matrix<NN,MM,TT>& A, Matrix<MM,MM,TT>& Q, Matrix<NN,MM,TT>& R){
			constexpr size_t K = std::min(NN,MM);
			Matrix<1,MM,TT> a;
			Matrix<1,MM,TT> e;
			Matrix<1,MM,TT> v;
			Matrix<MM,MM,TT> I = Matrix<MM,MM,TT>::eye();
			Matrix<MM,MM,TT> H;
			Matrix<NN,MM,TT> HA;
			HA = A;
			Q = I;
			
			if constexpr(std::is_same_v<TT, std::complex<double>>){
				// Householder Reflections を使ったQR分解アルゴリズム(複素数版)
				Matrix<1,1,TT> vHv;
				e[1] = std::complex(1.0, 0.0);
				for(size_t k = 1; k <= K; ++k){
					a = getcolumn(HA, k);
					a = shiftup(a, k - 1);
					v = a + sgn(a[1])*euclidnorm(a)*e;
					vHv = Htp(v)*v;
					if(k != 1) I.SetElement(MM - (k - 2), MM - (k - 2), 0);
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
				Matrix<1,1,TT> vTv;
				e[1] = 1;
				for(size_t k = 1; k <= K; ++k){
					a = getcolumn(HA, k);
					a = shiftup(a, k - 1);
					v = a + sgn(a[1])*euclidnorm(a)*e;
					vTv = tp(v)*v;
					if(k != 1) I.SetElement(MM - (k - 2), MM - (k - 2), 0);
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
		
		//! @brief SVD特異値分解(引数で返す版)
		//! 補足：MATLABとはU,S,Vの符号関係が入れ替わっている場合があるが正常なSVDであることは確認済み
		//! @param[in]	A	入力行列
		//! @param[out]	U	U行列
		//! @param[out]	S	S行列
		//! @param[out]	V	V行列
		constexpr friend void SVD(const Matrix<NN,MM,TT>& A, Matrix<MM,MM,TT>& U, Matrix<NN,MM,TT>& S, Matrix<NN,NN,TT>& V){
			constexpr size_t LoopMax = 100*std::max(NN,MM);	// ループ打ち切り最大回数
			Matrix<NN,MM,TT> Snm;
			Matrix<MM,NN,TT> Smn;
			Matrix<MM,MM,TT> Qm;
			Matrix<NN,NN,TT> Qn;
			Matrix<MM,NN,TT> e;
			double E = 0, F = 0;
			
			// 初期化
			U = Matrix<MM,MM,TT>::eye();
			Snm = A;
			V = Matrix<NN,NN,TT>::eye();
			TT Error = 1;
			
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
			Matrix<1,std::min(NN,MM),TT> Sd = diag(Smn);
			S = Matrix<NN,MM,TT>::zeros();
			for(size_t k = 1; k <= std::min(NN,MM); ++k){
				TT Sdn = Sd[k];
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
		constexpr friend std::tuple<Matrix<MM,MM,TT>, Matrix<NN,MM,TT>, Matrix<NN,NN,TT>> SVD(const Matrix<NN,MM,TT>& A){
			Matrix<MM,MM,TT> U;
			Matrix<NN,MM,TT> S;
			Matrix<NN,NN,TT> V;
			SVD(A, U, S, V);	// SVD特異値分解
			return {U, S, V};	// タプルで返す
		}
		
		//! @brief Schur分解
		//! @param[in]	A	入力行列
		//! @param[out]	Q	ユニタリ行列
		//! @param[out]	U	上三角行列または疑似上三角行列
		constexpr friend void Schur(const Matrix<NN,MM,TT>& A, Matrix<NN,MM,TT>& Q, Matrix<NN,MM,TT>& U){
			static_assert(A.N == A.M, "Matrix Size Error");				// Aが正方行列かチェック
			
			// QR法によるシュール分解
			Matrix<NN,MM,TT> Tk = A, Qk, Rk, Qa;
			Qa = Matrix<NN,MM,TT>::eye();
			for(size_t k = 0; k < ITERATION_MAX; ++k){
				QR(Tk, Qk, Rk);
				Tk = Rk*Qk;
				Qa = Qa*Qk;
				if(std::abs(Tk.GetElement(1,A.M)) < A.epsilon) break;
			}
			Q = Qa;
			U = Tk;
		}
		
		//! @brief Ax = bの形の線形連立1次方程式をxについて解く関数(引数で返す版)
		//! @param[in]	A	係数行列
		//! @param[in]	b	係数ベクトル
		//! @param[out]	x	解ベクトル
		constexpr friend void solve(const Matrix& A, const Matrix<1,MM,TT>& b, Matrix<1,NN,TT>& x){
			static_assert(A.N == A.M, "Matrix Size Error");				// Aが正方行列かチェック
			static_assert(b.M == A.M, "Matrix and Vector Size Error");	// Aの高さとbの高さが同じかチェック
			static_assert(b.N == 1, "Input is NOT vector.");			// bは縦ベクトルかチェック
			
			if constexpr(MM == 1){
				// スカラーの場合
				x[1] = b[1]/A.GetElement(1,1);	// スカラーのときはそのまま単純に除算
			}else{
				// 行列の場合
				// Ax = b において A をLU分解すると，(LU)x = b になって L(Ux) = b と表現できることを利用する。
				Matrix<A.N,A.N,TT> L, U;
				Matrix<1,A.N,TT> d, bb;
				Matrix<1,A.N,int> v;
				TT buff = 0;
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
		constexpr friend Matrix<1,MM,TT> solve(const Matrix& A, const Matrix<1,MM,TT>& b){
			Matrix<1,A.N,TT> x;
			solve(A, b, x);
			return x;		// 最終的な答えのxベクトルを返す
		}
		
		//! @brief Uは上三角行列で，Ux = bの形の線形連立1次方程式をxについて解く関数(引数で返す版)
		//! @param[in]	U	係数行列(上三角行列)
		//! @param[in]	b	係数ベクトル
		//! @param[out]	x	解ベクトル
		constexpr friend void solve_upper_tri(const Matrix& U, const Matrix<1,MM,TT>& b, Matrix<1,NN,TT>& x){
			static_assert(U.N == U.M, "Matrix Size Error");				// Uが正方行列かチェック
			static_assert(b.M == U.M, "Matrix and Vector Size Error");	// Uの高さとbの高さが同じかチェック
			static_assert(b.N == 1, "Input is NOT vector.");			// bは縦ベクトルかチェック
			
			if constexpr(MM == 1){
				// スカラーの場合
				x[1] = b[1]/U.GetElement(1,1);	// スカラーのときはそのまま単純に除算
			}else{
				// 行列の場合
				Matrix<1,U.N,int> v;
				TT buff = 0;
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
		constexpr friend TT det(const Matrix& A){
			static_assert(A.N == A.M, "Matrix Size Error");	// Aが正方行列かチェック
			Matrix<A.N,A.N,TT> L, U;
			Matrix<1,A.N,int> v;
			int sign;		// 符号
			if(LU(A, L, U, v) == Matrix::ODD){	// LU分解と符号判定
				sign = -1;	// 奇数のとき
			}else{
				sign =  1;	// 偶数のとき
			}
			// |A| = |L||U| でしかも |L|と|U|は対角要素の総積に等しく，さらにLの対角要素は1なので|L|は省略可。
			// 最後にLU分解のときの並べ替え回数によって符号反転をする。
			return (TT)sign*prod(U);
		}
		
		//! @brief 逆行列を返す関数 (正則チェック無し)
		//! @param[in]	A	入力行列
		//! @return	結果
		constexpr friend Matrix inv(const Matrix& A){
			static_assert(A.N == A.M, "Matrix Size Error");	// Aが正方行列かチェック
			Matrix I = Matrix<A.N,A.N,TT>::ident();			// 単位行列の生成
			Matrix<1,A.N,TT> x, b;
			Matrix<A.N,A.N,TT> Ainv;
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
		constexpr friend Matrix inv(const Matrix& A, size_t k){
			arcs_assert(k <= A.N);					// kが範囲内かチェック
			Matrix I = Matrix<A.N,A.N,TT>::ident();	// 単位行列の生成
			
			// 正則にするためにk列より右下の対角成分を「1」にする
			Matrix<A.N,A.N,TT> A2 = A;
			for(size_t j = k + 1; j <= A.N; ++j){
				A2.SetElement(j, j, 1);
			}
			
			// k列までの逆行列を計算
			Matrix<1,A.N,TT> x, b;
			Matrix<A.N,A.N,TT> Ainv;
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
		constexpr friend Matrix inv_with_check(const Matrix& A){
			static_assert(A.N == A.M, "Matrix Size Error");	// Aが正方行列かチェック
			arcs_assert(A.epsilon < std::abs(det(A)));		// 正則かチェック
			return inv(A);	// 最終的に得られる逆行列を返す
		}
		
		//! @brief 上三角行列の逆行列を返す関数
		//! @param[in]	U	入力行列(上三角行列)
		//! @param[out]	Uinv	逆行列
		constexpr friend void inv_upper_tri(const Matrix& U, Matrix& Uinv){
			static_assert(U.N == U.M, "Matrix Size Error");			// Uが正方行列かチェック
			static_assert(Uinv.N == Uinv.M, "Matrix Size Error");	// Uinvが正方行列かチェック
			static_assert(U.N == Uinv.N, "Matrix Size Error");		// UとUinvが同じサイズかチェック
			Matrix I = Matrix<U.N,U.N,TT>::ident();			// 単位行列の生成
			Matrix<1,U.N,TT> x, b;
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
		constexpr friend void inv_upper_tri(const Matrix& U, size_t k, Matrix& Uinv){
			static_assert(U.N == U.M, "Matrix Size Error");			// Uが正方行列かチェック
			static_assert(Uinv.N == Uinv.M, "Matrix Size Error");	// Uinvが正方行列かチェック
			static_assert(U.N == Uinv.N, "Matrix Size Error");		// UとUinvが同じサイズかチェック
			Matrix I = Matrix<U.N,U.N,TT>::ident();			// 単位行列の生成
			
			// 正則にするためにk列より右下の対角成分を「1」にする
			Matrix<U.N,U.N,TT> U2 = U;
			for(size_t j = k + 1; j <= U.N; ++j){
				U2.SetElement(j, j, 1);
			}
			
			// k列までの逆行列を計算
			Matrix<1,U.N,TT> x, b;
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
		constexpr friend Matrix<MM,NN,TT> lpinv(const Matrix& A){
			static_assert(A.N < A.M, "Matrix Size Error");	// 縦長行列かチェック
			return inv(tp(A)*A)*tp(A);
		}
		
		//! @brief 左擬似逆行列を返す関数 (Aが縦長行列の場合, 左上小行列のサイズ指定版)
		//! @param[in]	A	入力行列 (kより右と下は全部ゼロ埋めを想定)
		//! @param[in]	k	左上小行列のサイズ
		//! @return	結果
		constexpr friend Matrix<MM,NN,TT> lpinv(const Matrix& A, size_t k){
			Matrix<MM,NN> At = tp(A);
			Matrix<NN,NN> A2 = At*A;
			return inv(A2, k)*At;
		}
		
		//! @brief 右擬似逆行列を返す関数 (Aが横長行列の場合)
		//! @param[in]	A	入力行列
		//! @return	結果
		constexpr friend Matrix<MM,NN,TT> rpinv(const Matrix& A){
			static_assert(A.M < A.N, "Matrix Size Error");	// 横長行列かチェック
			return tp(A)*inv(A*tp(A));
		}
		
		//! @brief 右擬似逆行列を返す関数 (Aが横長行列の場合, 左上小行列のサイズ指定版)
		//! @param[in]	A	入力行列 (kより右と下は全部ゼロ埋めを想定)
		//! @param[in]	k	左上小行列のサイズ
		//! @return	結果
		constexpr friend Matrix<MM,NN,TT> rpinv(const Matrix& A, size_t k){
			Matrix<MM,NN> At = tp(A);
			Matrix<MM,MM> A2 = A*At;
			return At*inv(A2, k);
		}
		
		//! @brief 行列指数関数 e^(U)
		//! @param[in]	U	入力行列
		//! @param[in]	Order	パデ近似の次数
		//! @return	結果
		constexpr friend Matrix expm(const Matrix& U, size_t Order){
			static_assert(U.N == U.M, "Matrix Size Error");	// 正方行列かチェック
			int e = 0;
			TT c = 1;
			bool flag = false;
			// ノルムでスケーリング
			frexp(infnorm(U),&e);
			Matrix<U.N,U.M,TT> A;
			if(0 < e){
				A = pow(0.5,e+1)*U;
			}else{
				e = 0;
				A = 0.5*U;
			}
			// 行列のパデ近似の計算
			Matrix<A.N,A.N,TT> I = ident();// 単位行列の生成
			Matrix<A.N,A.N,TT> L = I, R = I, X = I, cX;
			for(size_t i = 1; i <= Order; ++i){
				c = c*(TT)(Order - i + 1)/(TT)(i*(2*Order - i + 1));	// パデ近似係数の計算
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
			Matrix Y = inv(L)*R;
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
		constexpr friend Matrix integral_expm(const Matrix& U, const TT T, const size_t DIV, const size_t P){
			static_assert(U.N == U.M, "Matrix Size Error");	// 正方行列かチェック
			const TT h = T/((TT)(2*DIV));	// 時間ステップ
			TT t = 0;						// 時刻
			// シンプソン法による定積分の実行
			Matrix<U.N,U.M,TT> S1, S2;
			for(size_t i = 1; i <= DIV; ++i){
				t = h*(TT)(2*i - 1);
				S1 += expm(U*t, P);
			}
			for(size_t i = 1; i <= DIV - 1; ++i){
				t = h*(TT)(2*i);
				S2 += expm(U*t, P);
			}
			return h/3.0*( Matrix<U.N,U.N,TT>::eye() + 4.0*S1 + 2.0*S2 + expm(U*T,P) );	// 最終的な定積分結果を返す
		}
		
		//! @brief 行列要素の指数関数を返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend Matrix expe(const Matrix& U){
			Matrix Y;
			for(size_t i = 0; i < U.N; ++i){
				for(size_t j = 0; j < U.M; ++j) Y.Data[i][j] = std::exp(U.Data[i][j]);
			}
			return Y;
		}
		
		//! @brief 行列要素の自然対数を返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend Matrix loge(const Matrix& U){
			Matrix Y;
			for(size_t i = 0; i < U.N; ++i){
				for(size_t j = 0; j < U.M; ++j) Y.Data[i][j] = std::log(U.Data[i][j]);
			}
			return Y;
		}
		
		//! @brief 行列要素の絶対値を返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend Matrix abse(const Matrix& U){
			Matrix Y;
			for(size_t i = 0; i < U.N; ++i){
				for(size_t j = 0; j < U.M; ++j) Y.Data[i][j] = std::abs(U.Data[i][j]);
			}
			return Y;
		}
		
		//! @brief 行列要素の平方根を返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend Matrix sqrte(const Matrix& U){
			Matrix Y;
			for(size_t i = 0; i < U.N; ++i){
				for(size_t j = 0; j < U.M; ++j) Y.Data[i][j] = std::sqrt(U.Data[i][j]);
			}
			return Y;
		}
		
		//! @brief 行列要素の平方根を参照で返す関数
		//! @param[in]	U	入力行列
		//! @param[out]	Y	結果
		constexpr friend void sqrte(const Matrix& U, Matrix& Y){
			for(size_t i = 0; i < U.N; ++i){
				for(size_t j = 0; j < U.M; ++j) Y.Data[i][j] = std::sqrt(U.Data[i][j]);
			}
		}
		
		//! @brief 行列要素のtanhを返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend Matrix tanhe(const Matrix& U){
			Matrix Y;
			for(size_t i = 0; i < U.N; ++i){
				for(size_t j = 0; j < U.M; ++j) Y.Data[i][j] = std::tanh(U.Data[i][j]);
			}
			return Y;
		}
		
		//! @brief 複素数行列要素の実数部を返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend Matrix<NN,MM,double> reale(const Matrix& U){
			static_assert(std::is_same_v<TT, std::complex<double>>, "Matrix Type Error");	// 複素数型のみ対応
			Matrix<NN,MM,double> Y;
			for(size_t i = 0; i < U.N; ++i){
				for(size_t j = 0; j < U.M; ++j) Y.Data[i][j] = U.Data[i][j].real();
			}
			return Y;
		}
		
		//! @brief 複素数行列要素の実数部に値をセットする関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr void real(const Matrix<NN,MM,double>& U){
			static_assert(N == NN, "Matrix Size Error");
			static_assert(M == MM, "Matrix Size Error");
			for(size_t i = 0; i < N; ++i){
				for(size_t j = 0; j < M; ++j) Data[i][j] = std::complex(U.Data[i][j], 0.0);
			}
		}
		
		//! @brief 複素数行列要素の虚数部を返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend Matrix<NN,MM,double> image(const Matrix& U){
			static_assert(std::is_same_v<TT, std::complex<double>>, "Matrix Type Error");	// 複素数型のみ対応
			Matrix<NN,MM,double> Y;
			for(size_t i = 0; i < U.N; ++i){
				for(size_t j = 0; j < U.M; ++j) Y.Data[i][j] = U.Data[i][j].imag();
			}
			return Y;
		}
		
		//! @brief 複素数行列要素の大きさを返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend Matrix<NN,MM,double> mage(const Matrix& U){
			static_assert(std::is_same_v<TT, std::complex<double>>, "Matrix Type Error");	// 複素数型のみ対応
			Matrix<NN,MM,double> Y;
			for(size_t i = 0; i < U.N; ++i){
				for(size_t j = 0; j < U.M; ++j) Y.Data[i][j] = std::abs(U.Data[i][j]);
			}
			return Y;
		}
		
		//! @brief 複素数行列要素の偏角を返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend Matrix<NN,MM,double> arge(const Matrix& U){
			static_assert(std::is_same_v<TT, std::complex<double>>, "Matrix Type Error");	// 複素数型のみ対応
			Matrix<NN,MM,double> Y;
			for(size_t i = 0; i < U.N; ++i){
				for(size_t j = 0; j < U.M; ++j) Y.Data[i][j] = std::arg(U.Data[i][j]);
			}
			return Y;
		}
		
		//! @brief 複素数行列要素の共役を返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend Matrix<NN,MM,std::complex<double>> conje(const Matrix& U){
			static_assert(std::is_same_v<TT, std::complex<double>>, "Matrix Type Error");	// 複素数型のみ対応
			Matrix<NN,MM,std::complex<double>> Y;
			for(size_t i = 0; i < U.N; ++i){
				for(size_t j = 0; j < U.M; ++j) Y.Data[i][j] = std::conj(U.Data[i][j]);
			}
			return Y;
		}
		
		//! @brief エルミート転置行列を返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend Matrix<MM,NN,std::complex<double>> Htp(const Matrix<NN,MM,TT>& U){
			static_assert(std::is_same_v<TT, std::complex<double>>, "Matrix Type Error");	// 複素数型のみ対応
			return conje(tp(U));	// 転置して複素共役
		}
		
		//! @brief 固有値を返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend Matrix<1,NN,std::complex<double>> eigen(const Matrix<NN,MM,TT>& U){
			static_assert(NN == MM, "Matrix Size Error");		// 正方行列のみ対応
			constexpr size_t LoopMax = 100*NN;					// ループ打ち切り最大回数
			auto I = Matrix<NN,NN,std::complex<double>>::eye();	// 単位行列
			Matrix<NN,NN,std::complex<double>> A, Q, R;
			std::complex<double> a, b, c, d, mu;
			
			if constexpr(std::is_same_v<TT, std::complex<double>>){
				// 入力が複素数型の場合
				A = U;
			}else{
				// 入力が実数型の場合
				A.real(U);
			}
			
			// 複素数QR法による固有値計算
			for(size_t k = 1; k < LoopMax; ++k){
				// Aの最右下2×2小行列の固有値を求める
				a = A.GetElement(NN - 1, NN - 1);
				b = A.GetElement(NN    , NN - 1);
				c = A.GetElement(NN - 1, NN    );
				d = A.GetElement(NN    , NN    );
				mu = ( (a + d) + std::sqrt((a + d)*(a + d) - 4.0*(a*d - b*c)) )/2.0;
				
				// QR分解と収束計算
				QR(A - mu*I, Q, R);
				A = R*Q + mu*I;
				
				if(std::abs(std::abs(tr(Q)) - (double)NN) < epsilon) break;	// 単位行列に漸近したらループ打ち切り
			}
			
			return diag(A);	// Aの対角要素が固有値
		}
		
		//! @brief 最大固有値の固有ベクトルを返す関数
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend Matrix<1,NN,std::complex<double>> eigenvec(const Matrix<NN,MM,TT>& U){
			static_assert(NN == MM, "Matrix Size Error");	// 正方行列のみ対応
			constexpr size_t LoopMax = 100*std::max(NN,MM);	// ループ打ち切り最大回数
			Matrix<NN,NN,std::complex<double>> A;
			
			if constexpr(std::is_same_v<TT, std::complex<double>>){
				// 入力が複素数型の場合
				A = U;
			}else{
				// 入力が実数型の場合
				A.real(U);
			}
			
			// べき乗法による固有ベクトル計算
			auto x = Matrix<1,NN,std::complex<double>>::ones();
			auto y = Matrix<1,NN,std::complex<double>>::ones();
			for(size_t k = 1; k < LoopMax; ++k){
				y = A*x;
				x = y/euclidnorm(y);
			}
			
			return x;
		}
		
		//! @brief クロネッカー積
		//! @tparam	PP	右側の幅
		//! @tparam QQ	右側の高さ
		//! @param[in] Ul 演算子の左側
		//! @param[in] Ur 演算子の右側
		//! @return 結果
		template<size_t PP, size_t QQ>
		constexpr friend Matrix<NN*PP, MM*QQ, TT> Kronecker(const Matrix<NN,MM,TT>& Ul, const Matrix<PP,QQ,TT>& Ur){
			Matrix<NN*PP, MM*QQ, TT> Y;
			Matrix<PP,QQ,TT> A;
			
			// 縦方向に小行列で埋めていく
			for(size_t j = 1; j <= MM; ++j){
				// 横方向に小行列で埋めていく
				for(size_t i = 1; i <= NN; ++i){
					A = Ul(j,i)*Ur;
					setsubmatrix(Y, (i - 1)*PP + 1, (j - 1)*QQ + 1, A);
				}
			}
			
			return Y;
		}
		
		//! @brief vec作用素(行列→縦ベクトル)
		//! @param[in]	U	入力行列
		//! @return	結果
		constexpr friend Matrix<1,NN*MM,TT> vec(const Matrix<NN,MM,TT>& U){
			Matrix<1,NN*MM,TT> Y;
			
			// 横方向に走査
			size_t k = 0;
			for(size_t i = 1; i <= NN; ++i){
				// 縦方向に走査
				for(size_t j = 1; j <= MM; ++j){
					k++;
					Y(k,1) = U(j,i);
				}
			}
			
			return Y;
		}
		
		//! @brief vec作用素の逆(縦ベクトル→行列)
		//! @tparam	PP	再構成後の幅
		//! @tparam QQ	再構成後の高さ
		//! @param[in]		v	入力ベクトル
		//! @param[in,out]	Y	再構成後の行列
		template<size_t PP, size_t QQ>
		constexpr static void vecinv(const Matrix<NN,MM,TT>& v, Matrix<PP,QQ,TT>& Y){
			static_assert(NN == 1);		// 入力は縦ベクトルかチェック
			static_assert(PP*QQ == MM);	// 要素数が同じかチェック
			
			// 横方向に走査
			size_t k = 0;
			for(size_t i = 1; i <= PP; ++i){
				// 縦方向に走査
				for(size_t j = 1; j <= QQ; ++j){
					k++;
					Y(j,i) = v(k,1);
				}
			}
		}
		
		//! @brief vec作用素の逆(縦ベクトル→行列)
		//! @tparam	PP	再構成後の幅
		//! @tparam QQ	再構成後の高さ
		//! @param[in]		v	入力ベクトル
		//! @return	再構成後の行列
		template<size_t PP, size_t QQ>
		constexpr static Matrix<PP,QQ,TT> vecinv(const Matrix<NN,MM,TT>& v){
			Matrix<PP,QQ,TT> Y;
			vecinv(v, Y);
			return Y;
		}
		
	private:
		static constexpr double epsilon = 1e-12;		//!< 零とみなす閾値(実数版)
		static constexpr std::complex<double> epscomp = std::complex(1e-12, 1e-12);	//!< 零とみなす閾値(複素数版)
		static constexpr size_t ITERATION_MAX = 10000;	//!< 反復計算の最大値
		size_t Nindex;	//!< 横用カウンタ
		size_t Mindex;	//!< 縦用カウンタ
		
		//! @brief 符号関数
		//! @param[in]	u	入力
		//! @return	符号結果
		static constexpr TT sgn(TT u){
			TT ret = 0;
			if constexpr(std::is_same_v<TT, std::complex<double>>){
				// 複素数型の場合
				if(u == std::complex(0.0, 0.0)){
					ret = std::complex(0.0, 0.0);
				}else{
					ret = u/std::abs(u);
				}
			}else{
				// 実数型の場合
				if((TT)0 <= u){
					ret = (TT)1;
				}else{
					ret = (TT)(-1);
				}
			}
			return ret;
		}
		
	public:
		static constexpr size_t N = NN;			//!< 行列の幅(列の数, 横)
		static constexpr size_t M = MM;			//!< 行列の高さ(行の数, 縦)
		std::array<std::array<TT, M>, N> Data;	//!< データ格納用変数 配列要素の順番は Data[N列(横)][M行(縦)] なので注意
};
}

#endif


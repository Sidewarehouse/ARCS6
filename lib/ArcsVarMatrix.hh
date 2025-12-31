//! @file ArcsVarMatrix.hh
//! @brief ARCS-Variable-Matrix 可変行列演算クラス
//!
//! 可変サイズの行列に関係する簡単な演算を実行するクラス
//! ArcsMatのヒープ領域版
//!
//! @date 2025/12/31
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2025 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.
//
// ・ヒープ領域を使うためArcsMatに比べて低速

#ifndef ARCSVARMATRIX
#define ARCSVARMATRIX

#include <string>
#include <cmath>
#include <cassert>
#include <vector>
#include <cstdint>

// ARCS組込み用マクロ
#ifdef ARCS_IN
	// ARCSに組み込まれる場合
	#include "ARCSassert.hh"
#else
	// ARCSに組み込まれない場合
	#define arcs_assert(a) (assert(a))
#endif

#include "ArcsMatrix.hh"

// ARCS名前空間
namespace ARCS {

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

		/*		
		
		//! @brief 行列加算演算子(行列＝行列＋行列の場合)
		//! @tparam	P, Q, R	演算子右側の行列の高さ, 幅, 要素の型
		//! @param[in] right 演算子の右側
		//! @return 結果
		template<size_t P, size_t Q, typename R = double>
		constexpr ArcsMat<M,N,T> operator+(const ArcsMat<P,Q,R>& right) const{
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsVarMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
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
			static_assert(ArcsVarMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
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
			static_assert(ArcsVarMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
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
			static_assert(ArcsVarMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
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
			static_assert(ArcsVarMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
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
			static_assert(ArcsVarMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
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
			static_assert(ArcsVarMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
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
			static_assert(ArcsVarMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
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
			static_assert(ArcsVarMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
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
			static_assert(ArcsVarMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			for(size_t j = 0; j < M; ++j) Data[0][j] = Array[j];
		}
		
		//! @brief 縦ベクトルを1次元std::array配列に書き込む関数
		//! @tparam	P, R	配列の長さ, 要素の型
		//! @param[out]	Array	std::array配列(縦MM×横1)
		template<size_t P, typename R = double>
		constexpr void StoreArray(std::array<R, P>& Array) const{
			static_assert(N == 1, "ArcsMat: Vector Error");	// 縦ベクトルチェック
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsVarMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			for(size_t j = 0; j < M; ++j) Array[j] = Data[0][j];
		}
		
		//! @brief 2次元std::array配列を行列として読み込む関数
		//! @tparam	P, Q, R	配列の高さ, 幅, 要素の型
		//! @param[in]	Array	std::array配列(縦MM×横NN)
		template<size_t P, size_t Q, typename R = double>
		constexpr void LoadArray(const std::array<std::array<R, P>, Q>& Array){
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsVarMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			Data = Array;
		}
		
		//! @brief 行列を2次元std::array配列に書き込む関数
		//! @tparam	P, Q, R	配列の高さ, 幅, 要素の型
		//! @param[in]	Array	std::array配列(縦MM×横NN)
		template<size_t P, size_t Q, typename R = double>
		constexpr void StoreArray(std::array<std::array<R, P>, Q>& Array) const{
			static_assert(M == P, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(N == Q, "ArcsMat: Size Error");	// 行列のサイズチェック
			static_assert(ArcsVarMatrix::IsApplicable<R>, "ArcsMat: Type Error");	// 対応可能型チェック
			Array = Data;
		}
		
*/		

};

}

#endif

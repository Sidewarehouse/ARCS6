//! @file RingBuffer.cc
//! @brief リングバッファクラス(テンプレート版)
//! @date 2022/03/25
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the BSD License.
// For details, see the License.txt file.
//
// リングバッファの動作例
//                                                                GetRelativeValueFromFirst(k=13)
//                                       SetFirstValue (進む方向→)             |
//                                               |                              |
// [0][j=i+1+k-N=1][2][3][4][5][6][7][8][9][10][i=11][j=i+1=12][13][j=i+1+k=14][15][N-1=16]  (N=17)
//          |                         |          |         |              |
//          |                         |          | GetFinalValue          |
//    GetRelativeValueFromEnd(k=6)    |          |             GetRelativeValueFromEnd(k=2)
//                                    |    GetFirstValue
//                         GetRelativeValueFromFirst(k=3)

#ifndef RINGBUFFER
#define RINGBUFFER

#include <cassert>
#include <pthread.h>
#include <array>
#include <string>
#include "ARCSassert.hh"
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
//! @brief リングバッファクラス
//! @tparam	T	型
//! @tparam	N	バッファサイズ
//! @tparam	M	Mutexロックを使うかどうか(デフォルトはtrue)
template <typename T, unsigned long N, bool M = true>
class RingBuffer {
	public:
		std::array<T, N> Buffer;	//!< リングバッファ
		
		//! @brief コンストラクタ
		// Length；リングバッファの大きさ
		RingBuffer()
			: Buffer(), i(0), BufMutex()
		{
			PassedLog();
			if constexpr(M == true) pthread_mutex_init(&BufMutex, NULL);// Mutexの初期化
			ClearBuffer();						// バッファクリア
		}
		
		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		RingBuffer(RingBuffer&& r)
			: i(r.i), Buffer(r.Buffer), BufMutex(r.BufMutex)
		{
			
		}
		
		//! @brief デストラクタ
		~RingBuffer(){
			PassedLog();
		}
		
		//! @brief 値をバッファの先頭に格納する関数
		//! @param[in]	u	入力値
		void SetFirstValue(const T& u){
			if constexpr(M == true) pthread_mutex_lock(&BufMutex);	// Mutexロック
			++i;				// リングバッファカウンタ
			if(N <= i) i = 0;	// バッファの最後まで行ったらカウンタリセット
			Buffer[i] = u;		// バッファに入力値を格納
			if constexpr(M == true) pthread_mutex_unlock(&BufMutex);// Mutexロック解除
		}
		
		//! @brief 値をバッファの現在の先頭から取り出す関数
		//! @return	出力値
		T GetFirstValue(void){
			T ret;
			if constexpr(M == true) pthread_mutex_lock(&BufMutex);	// Mutexロック
			ret = Buffer[i];		// バッファから読み出す
			if constexpr(M == true) pthread_mutex_unlock(&BufMutex);// Mutexロック解除
			return ret;
		}
		
		//! @brief バッファの最先端から k だけ後方側に戻ったところの値を取り出す関数
		//! @param[in]	k	相対要素番号
		T GetRelativeValueFromFirst(const unsigned long k){
			T ret;
			long j;
			if constexpr(M == true) pthread_mutex_lock(&BufMutex);	// Mutexロック
			j = i - k;			// 今のリングバッファカウンタから相対値を減算
			if(j < 0) j = j + N;// バッファ要素番号の最初を下回っていたら，先頭＋はみ出た分にセット
			arcs_assert((unsigned long)j < N);
			ret = Buffer[j];
			if constexpr(M == true) pthread_mutex_unlock(&BufMutex);// Mutexロック解除
			return ret;
		}
		
		//! @brief バッファの最後尾から k だけ先頭側に進んだところの値を取り出す関数
		//! @param[in]	k	相対要素番号
		T GetRelativeValueFromEnd(const unsigned long k){
			T ret;
			unsigned long j;
			if constexpr(M == true) pthread_mutex_lock(&BufMutex);	// Mutexロック
			j = i + 1 + k;			// 今のリングバッファカウンタに相対値を加算
			if(N <= j) j = j - N;	// バッファ要素番号の最後を超えていたら，先頭＋はみ出た分にセット
			arcs_assert(j < N);
			ret = Buffer[j];
			if constexpr(M == true) pthread_mutex_unlock(&BufMutex);// Mutexロック解除
			return ret;
		}
		
		//! @brief バッファから最後尾の値を取り出す関数
		//! @return	出力値
		T GetFinalValue(void){
			T ret;
			unsigned long j;
			if constexpr(M == true) pthread_mutex_lock(&BufMutex);	// Mutexロック
			if(i == N - 1){	// バッファカウンタがリングの終端に達していたら，
				j = 0;		// リングの先端へ
			}else{
				j = i + 1;	// それ以外は最後尾へ
			}
			ret = Buffer[j];
			if constexpr(M == true) pthread_mutex_unlock(&BufMutex);// Mutexロック解除
			return ret;
		}
		
		//! @brief カウンタを任意値に設定する関数
		//! @param[in]	j	任意のカウント値
		void SetCounter(const unsigned long j){
			if constexpr(M == true) pthread_mutex_lock(&BufMutex);	// Mutexロック
			i = j;					// カウンタを任意値にセット
			if(N <= i) i = N - 1;	// バッファの最後を超えてたら終端でリミット
			if constexpr(M == true) pthread_mutex_unlock(&BufMutex);// Mutexロック解除
		}
		
		//! @brief カウンタをリセットする関数
		void ResetCounter(void){
			if constexpr(M == true) pthread_mutex_lock(&BufMutex);	// Mutexロック
			i = 0;
			if constexpr(M == true) pthread_mutex_unlock(&BufMutex);// Mutexロック解除
		}
		
		//! @brief バッファを指定値で埋める
		//! @param[in]	u	指定値
		void FillBuffer(const T& u){
			if constexpr(M == true) pthread_mutex_lock(&BufMutex);	// Mutexロック
			for(unsigned long j = 0; j < N; ++j) Buffer[j] = u;
			if constexpr(M == true) pthread_mutex_unlock(&BufMutex);// Mutexロック解除
		}
		
		//! @brief バッファのゼロクリア
		void ClearBuffer(void){
			if constexpr(std::is_same_v<T, std::string>){
				// 文字列型のとき
				FillBuffer("");	// 空文字で埋める
			}else if constexpr(std::is_floating_point_v<T>){
				// 浮動小数点型のとき
				FillBuffer(0);	// ゼロで埋める
			}else if constexpr(std::is_integral_v<T>){
				// 整数型のとき
				FillBuffer(0);	// ゼロで埋める
			}else{
				// それ以外のときは行列型と判断
				for(size_t i = 0; i < N; ++i) Buffer.at(i).FillAllZero();	// バッファ内の全行列をクリア
			}
		}
		
	private:
		RingBuffer(const RingBuffer&) = delete;					//!< コピーコンストラクタ使用禁止
		const RingBuffer& operator=(const RingBuffer&) = delete;//!< 代入演算子使用禁止
		unsigned long i;			//!< リングバッファカウンタ
		pthread_mutex_t BufMutex;	//!< バッファアクセス用Mutex
};
}

#endif


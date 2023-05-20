//! @file DataMemory.hh
//! @brief データメモリクラス
//!
//! 実験データを保存してCSVファイルへの出力を行うクラス。
//!
//! @date 2020/06/11
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2023 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#ifndef DATAMEMORY
#define DATAMEMORY

#include <array>
#include <memory>
#include <cmath>
#include "ConstParams.hh"

namespace ARCS {	// ARCS名前空間
	//! @brief データメモリクラス
	class DataMemory {
		public:
			DataMemory();				//!< コンストラクタ
			~DataMemory();				//!< デストラクタ
			void Reset(void);			//!< リセットする関数
			void WriteCsvFile(void);	//!< CSVファイルを書き出す関数
			
			//! @brief データを格納する関数(可変長引数テンプレート)
			//! @param[in] u1...u2 インジケータの値
			template<typename T1, typename... T2>
			void SetData(const T1& u1, const T2&... u2){
				// 再帰で順番に可変長引数を読み込んでいく
				if(Nindex == 0){
					// 1個目の引数は周期として一時的に格納
					Tperiod = (double)u1;
				}else if(Nindex == 1){
					// 2個目の引数は時刻として一時的に格納
					Time = (double)u1;
					if(
						Time < ConstParams::DATA_START || ConstParams::DATA_END <= Time	||
						Tperiod < fmod(Time - ConstParams::DATA_START, ConstParams::DATA_RESO)
					){
						// 保存時間の範囲外，もしくは保存時刻でなかったら
						Nindex = 0;	// 再帰カウンタを零に戻しておいてから
						return;		// 可変長引数の読み込みを打ち切る
					}
				}else{
					// 3個目以降は変数値として処理
					if(0 <= Tindex && Tindex < ELEMENT_NUM){
						// 念のための要素番号リミッタ
						SaveBuffer->at(Tindex).at(Nindex - 1) = (double)u1;	// 2列目以降に変数値を保存
					}
				}
				++Nindex;		// 再帰カウンタをインクリメント
				SetData(u2...);	// 自分自身を呼び出す(再帰)
			}
			//! @brief 再帰の最後に呼ばれる関数
			void SetData(){
				if(0 <= Tindex && Tindex < ELEMENT_NUM){
					// 念のための要素番号リミッタ
					SaveBuffer->at(Tindex).at(0) = Time;// 1列目に時刻を保存
				}
				++Tindex;							// 時間用カウンタを進める
				Nindex = 0;		// 可変長引数の読み込み作業が終わったので，再帰カウンタを零に戻しておく
			}
			
		private:
			DataMemory(DataMemory&& r) = delete;					//!< ムーブコンストラクタ使用禁止
			DataMemory(const DataMemory&) = delete;					//!< コピーコンストラクタ使用禁止
			const DataMemory& operator=(const DataMemory&) = delete;//!< 代入演算子使用禁止
			
			//! @brief [-] データの要素数(丸め誤差対策のために念のため1個分メモリを多く確保しておく)
			static constexpr size_t ELEMENT_NUM = (ConstParams::DATA_END - ConstParams::DATA_START)/ConstParams::DATA_RESO + 1;
			
			//! @brief データバッファへのスマートポインタ
			//! 巨大なデータを扱えるようにヒープ領域に確保する
			std::unique_ptr<std::array<
				std::array<double, ConstParams::DATA_NUM>,
				ELEMENT_NUM
			>> SaveBuffer;
			size_t Tindex;		//!< 時間用カウンタ
			size_t Nindex;		//!< 変数用カウンタ
			double Tperiod;		//!< [s] 周期
			double Time;		//!< [s] 時刻
	};
}

#endif


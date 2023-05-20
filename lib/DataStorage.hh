//! @file DataStorage.hh
//! @brief データ格納/ファイル出力クラス(テンプレート版)
//!
//! データを一時的に格納し，CSV/TAB区切りDATファイルとして出力します。
//!
//! @date 2020/03/11
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2020 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef DATASTORAGE
#define DATASTORAGE

#include <cassert>
#include <array>
#include <fstream>

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
	//! @brief ファイル形式の定義
	enum class DataStorageType {
		FORMAT_CSV,	//!< CSVファイルとして保存
		FORMAT_DAT	//!< DAT(タブ区切り)ファイルとして保存
	};
	
	//! @brief データ格納/ファイル出力クラス(テンプレート版)
	//! @tparam	NumOfVar 変数の数，SmplTime [ns] 制御周期，SaveTime [s] 保存時間，Rate 間引き数
	template <size_t NumOfVar, unsigned long SmplTime, unsigned int SaveTime, unsigned int Rate>
	class DataStorage {
		public:
			
			//! @brief コンストラクタ
			DataStorage()
				: Data({0}), Buff({0}), Nindex(0), Tindex(0), Mcount(0)
			{
				PassedLog();	// イベントログにココを通過したことを記録
			}

			//! @brief デストラクタ
			~DataStorage(){
				PassedLog();	// イベントログにココを通過したことを記録
			}
			
			//! @brief データを格納する関数
			//! @param[in] u1...u2 インジケータの値
			template<typename T1, typename... T2>		// 可変長引数テンプレート
			void SetData(const T1& u1, const T2&... u2){
				// 再帰で順番に可変長引数を読み込んでいく
				if(Nindex < N) Buff.at(Nindex) = (double)u1;	// 配列の要素を埋める
				++Nindex;		// 再帰カウンタをインクリメント
				SetData(u2...);	// 自分自身を呼び出す(再帰)
			}
			void SetData(){
				// 再帰の最後に呼ばれる関数
				Nindex = 0;	// 可変長引数の読み込み作業が終わったので，再帰カウンタを零に戻しておく
				++Mcount;	// 間引き用カウンタ
				if(Mcount % Rate == 0){	// 間引き条件計算
					// データの格納
					if(Tindex < T){	// 設定時間の範囲内であって，
						for(size_t n = 0; n < N; ++n){
							Data[n][Tindex] = Buff.at(n);	// データを保持
						}
						++Tindex;	// 時間用カウンタを進める
					}
				}
			}
			
			//! @brief CSV/DATファイルの保存
			//! @param[in] FileName ファイル名
			//! @param[in] FileType ファイル形式
			void SaveDataFile(const std::string& FileName, const DataStorageType FileType) const {
				// ファイルストリーム
				std::ofstream fout(FileName.c_str(), std::ios::out | std::ios::trunc);
				
				// エラーチェック
				arcs_assert(fout.bad() == false);	// 致命的なエラーの場合
				arcs_assert(fout.fail() == false);	// ファイル開くのに失敗した場合
				
				// 書式設定
				fout.setf(std::ios::scientific);	// 指数表示
				fout.width(15);						// 総桁数
				fout.precision(14);					// 小数点桁数
				
				// 実験データの書き出し
				PassedLog();	// イベントログにココを通過したことを記録
				for(size_t j = 0; j < Tindex; ++j){		// 時間分だけ回す
					for(size_t i = 0; i < N; ++i){		// 変数の数だけ回す
						fout << Data[i][j];	// 数値の書き出し
						if(i < N-1){		// 最後の列以外のときは，
							if(FileType == DataStorageType::FORMAT_CSV) fout << ',';	// CSVの場合はコンマで区切る
							if(FileType == DataStorageType::FORMAT_DAT) fout << '\t';	// DATの場合はタブで区切る
						}
					}
					fout << std::endl;		// 改行
				}
				PassedLog();	// イベントログにココを通過したことを記録
			}
			
			//! @brief 要素数カウンタをクリア
			void ClearCounter(void){
				Tindex = 0;
				Mcount = 0;
			}
			
		private:
			DataStorage(const DataStorage&) = delete;					// コピーコンストラクタ使用禁止
			const DataStorage& operator=(const DataStorage&) = delete;	// 代入演算子使用禁止
			
			static constexpr size_t N = NumOfVar;					//!< 変数方向のデータ数
			static constexpr size_t T = SaveTime*1e9/SmplTime/Rate;	//!< 時間方向のデータ数
			std::array<std::array<double, T>, N> Data;	//!< データ格納用変数
			std::array<double, N> Buff;					//!< データバッファ
			size_t Nindex;	//!< 変数用カウンタ
			size_t Tindex;	//!< 時間用カウンタ
			size_t Mcount;	//!< 間引き用カウンタ
	};
}

#endif


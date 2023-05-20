//! @file CsvManipulator.hh
//! @brief CSVファイル操作クラス
//!
//! std::arrayやMatrixの値をCSVファイルとして読み書きしたりするクラス
//!
//! @date 2020/05/24
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef CSVMANIPULATOR
#define CSVMANIPULATOR

#include <cassert>
#include <array>
#include <memory>
#include <string>
#include <fstream>
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
//! @brief 数値表現の定義
enum class CsvExpression {
	NORMAL,		//!< 通常表記
	EXPONENTIAL,//!< 指数表示 (6.33e+4のような表記)
	HEXFLOAT	//!< 16進数浮動小数点 (完全に正確な浮動小数点表記)
};

//! @brief CSVファイル操作クラス
class CsvManipulator {
	public:
		//! @brief 1次元std::arrayをCSVファイルに書き出す関数
		//! @tparam	E	浮動小数点数の表現方法(デフォルトは指数表記)
		//! @tparam	M	配列の長さ
		//! @param[in]	Data	出力したい行列
		//! @param[in]	FileName	ファイル名
		template <CsvExpression E = CsvExpression::EXPONENTIAL, size_t M>
		static void SaveFile(const std::array<double, M>& Data, const std::string& FileName){	
			std::ofstream fout(FileName.c_str(), std::ios::out | std::ios::trunc);	// ファイル出力ストリーム
			CheckError(fout);		// エラーチェック
			SetExpression<E>(fout);	// 浮動小数点数の表現設定
			
			// CSVデータの書き出し
			for(size_t j = 0; j < M; ++j){	// 行数分だけ回す
				fout << Data.at(j);			// 数値の書き出し
				fout << std::endl;			// 改行
			}
		}
		
		//! @brief CSVファイルから1次元std::arrayに読み込む関数
		//! @tparam	T	データ型
		//! @tparam	M	配列の長さ
		//! @param[in]	Data	入力したい行列
		//! @param[in]	FileName	ファイル名
		template <typename T = double, size_t M>
		static void LoadFile(std::array<T,M>& Data, const std::string& FileName){
			std::ifstream fin(FileName.c_str());	// ファイル入力ストリーム
			CheckError(fin, FileName);				// エラーチェック
			
			// CSVデータの読み込み
			std::string ReadBuff = "";				// 読み込みバッファ
			size_t j = 0;							// 行カウンタ
			while(getline(fin, ReadBuff)){			// 1行読み込み
				++j;								// 行カウンタ
				arcs_assert(j <= M);				// 行カウンタ溢れチェック
				
				// データ型に従った処理
				if constexpr(std::is_same_v<T,double>){
					Data.at(j - 1) = std::stod(ReadBuff);	// 浮動小数点数値に変換してから配列データに書き込み
				}
				if constexpr(std::is_same_v<T,std::string>){
					Data.at(j - 1) = ReadBuff;		// そのまま配列データに書き込み
				}
				
			}
		}
		
		//! @brief 2次元std::arrayをCSVファイルに書き出す関数(書き出すサイズを指定する版)
		//! @tparam	E	浮動小数点数の表現方法(デフォルトは指数表記)
		//! @tparam N	配列の横幅
		//! @tparam	M	配列の縦の長さ
		//! @param[in]	Data	出力したい行列
		//! @param[in]	FileName	ファイル名
		//! @param[in]	NN	書き出す横幅
		//! @param[in]	MM	書き出す縦の長さ
		template <CsvExpression E = CsvExpression::EXPONENTIAL, size_t N, size_t M>
		static void SaveFile(const std::array<std::array<double, N>, M>& Data, const std::string& FileName, size_t NN, size_t MM){	
			std::ofstream fout(FileName.c_str(), std::ios::out | std::ios::trunc);	// ファイル出力ストリーム
			CheckError(fout);		// エラーチェック
			SetExpression<E>(fout);	// 浮動小数点数の表現設定
			
			// CSVデータの書き出し
			for(size_t j = 0; j < MM; ++j){		// 行数分だけ回す
				for(size_t i = 0; i < NN; ++i){	// 列数分だけ回す
					fout << Data.at(j).at(i);	// 数値の書き出し
					if(i < NN - 1){		// 最後の列以外のときは，
						fout << ',';	// コンマで区切る
					}
				}
				fout << std::endl;		// 改行
			}
		}
		
		//! @brief 2次元std::arrayをCSVファイルに書き出す関数
		//! @tparam	E	浮動小数点数の表現方法(デフォルトは指数表記)
		//! @tparam N	配列の横幅
		//! @tparam	M	配列の縦の長さ
		//! @param[in]	Data	出力したい行列
		//! @param[in]	FileName	ファイル名
		template <CsvExpression E = CsvExpression::EXPONENTIAL, size_t N, size_t M>
		static void SaveFile(const std::array<std::array<double, N>, M>& Data, const std::string& FileName){
			SaveFile<E, N, M>(Data, FileName, N, M);
		}
		
		//! @brief 2次元std::arrayをCSVファイルに書き出す関数(スマートポインタ＆書き出すサイズを指定する版)
		//! @tparam	E	浮動小数点数の表現方法(デフォルトは指数表記)
		//! @tparam N	配列の横幅
		//! @tparam	M	配列の縦の長さ
		//! @param[in]	Data	出力したい行列
		//! @param[in]	FileName	ファイル名
		//! @param[in]	NN	書き出す横幅
		//! @param[in]	MM	書き出す縦の長さ
		template <CsvExpression E = CsvExpression::EXPONENTIAL, size_t N, size_t M>
		static void SaveFile(std::unique_ptr< std::array<std::array<double, N>, M> >&& Data, const std::string& FileName, size_t NN, size_t MM){	
			std::ofstream fout(FileName.c_str(), std::ios::out | std::ios::trunc);	// ファイル出力ストリーム
			CheckError(fout);		// エラーチェック
			SetExpression<E>(fout);	// 浮動小数点数の表現設定
			
			// CSVデータの書き出し
			for(size_t j = 0; j < MM; ++j){		// 行数分だけ回す
				for(size_t i = 0; i < NN; ++i){	// 列数分だけ回す
					fout << Data->at(j).at(i);	// 数値の書き出し
					if(i < NN - 1){		// 最後の列以外のときは，
						fout << ',';	// コンマで区切る
					}
				}
				fout << std::endl;		// 改行
			}
		}
		
		//! @brief CSVファイルから2次元std::arrayに読み込む関数
		//! @tparam N	配列の横幅
		//! @tparam	M	配列の縦の長さ
		//! @param[in]	Data	入力したい行列
		//! @param[in]	FileName	ファイル名
		template <size_t N, size_t M>
		static void LoadFile(std::array<std::array<double, N>, M>& Data, const std::string& FileName){
			std::ifstream fin(FileName.c_str());	// ファイル入力ストリーム
			CheckError(fin, FileName);						// エラーチェック
			
			// CSVデータの読み込み
			std::string ReadBuff = "";				// 読み込みバッファ
			std::array<double, N> ValBuff;			// 数値配列バッファ
			size_t j = 0;							// 行カウンタ
			while(getline(fin, ReadBuff)){			// 1行読み込み
				++j;								// 行カウンタ
				arcs_assert(j <= M);				// 行カウンタ溢れチェック
				ParseCsvString(ReadBuff, ValBuff);	// CSVパース＆読み込み
				for(size_t i = 0; i < N; ++i){
					Data.at(j-1).at(i) = ValBuff.at(i);
				}
			}
		}
		
		//! @brief 行列をCSVファイルに書き出す関数
		//! @tparam	E	浮動小数点数の表現方法(デフォルトは指数表記)
		//! @tparam N	行列の横幅
		//! @tparam	M	行列の縦の長さ
		//! @param[in]	Data	出力したい行列
		//! @param[in]	FileName	ファイル名
		template <CsvExpression E = CsvExpression::EXPONENTIAL, size_t N, size_t M>
		static void SaveFile(const Matrix<N,M>& Data, const std::string& FileName){	
			std::ofstream fout(FileName.c_str(), std::ios::out | std::ios::trunc);	// ファイル出力ストリーム
			CheckError(fout);		// エラーチェック
			SetExpression<E>(fout);	// 浮動小数点数の表現設定
			
			// CSVデータの書き出し
			for(size_t j = 0; j < M; ++j){				// 行数分だけ回す
				for(size_t i = 0; i < N; ++i){			// 列数分だけ回す
					fout << Data.GetElement(i+1,j+1);	// 数値の書き出し
					if(i < N - 1){		// 最後の列以外のときは，
						fout << ',';	// コンマで区切る
					}
				}
				fout << std::endl;		// 改行
			}
		}
		
		//! @brief CSVファイルから行列に読み込む関数
		//! @tparam N	行列の横幅
		//! @tparam	M	行列の縦の長さ
		//! @param[in]	Data	入力したい行列
		//! @param[in]	FileName	ファイル名
		template <size_t N, size_t M>
		static void LoadFile(Matrix<N,M>& Data, const std::string& FileName){
			std::ifstream fin(FileName.c_str());	// ファイル入力ストリーム
			CheckError(fin, FileName);						// エラーチェック
			
			// CSVデータの読み込み
			std::string ReadBuff = "";				// 読み込みバッファ
			std::array<double, N> ValBuff;			// 数値配列バッファ
			size_t j = 0;						// 行カウンタ
			while(getline(fin, ReadBuff)){			// 1行読み込み
				++j;								// 行カウンタ
				arcs_assert(j <= M);				// 行カウンタ溢れチェック
				ParseCsvString(ReadBuff, ValBuff);	// CSVパース＆読み込み
				setrow(Data, ValBuff, j);			// 行を書き込む
			}
		}
		
	private:
		CsvManipulator() = delete;					//!< コンストラクタ使用禁止
		CsvManipulator(CsvManipulator&& r) = delete;//!< ムーブコンストラクタ使用禁止
		~CsvManipulator() = delete;					//!< デストラクタ使用禁止
		CsvManipulator(const CsvManipulator&) = delete;					//!< コピーコンストラクタ使用禁止
		const CsvManipulator& operator=(const CsvManipulator&) = delete;//!< 代入演算子使用禁止
		
		//! @brief ファイル書き込みエラーのチェック
		//! @param[in]	fileout	ファイル出力ストリーム
		static void CheckError(const std::ofstream& fileout){
			arcs_assert(fileout.bad() == false);	// 致命的なエラーの場合
			arcs_assert(fileout.fail() == false);	// ファイルを開くのに失敗した場合
		}
		
		//! @brief ファイル読み込みエラーのチェック
		//! @param[in]	filein	ファイル入力ストリーム
		//! @param[in]	filename	ファイル名
		static void CheckError(const std::ifstream& filein, const std::string& filename){
			EventLog("Checking CSV File Read Error of " + filename + "...");
			arcs_assert(filein.is_open() == true);	// ファイルを開くのに失敗した場合
			EventLog("Checking CSV File Read Error of " + filename + "...Done");
		}
		
		//! @brief 浮動小数点数の表現方法の設定
		//! @tparam	E	浮動小数点数の表現方法
		//! @param[in]	fileout	ファイル出力ストリーム
		template <CsvExpression E>
		static void SetExpression(std::ofstream& fileout){
			if constexpr(E == CsvExpression::NORMAL){
				// 通常表記の場合
				fileout.setf(std::ios::fixed);
			}
			if constexpr(E == CsvExpression::EXPONENTIAL){
				// 指数表記の場合
				fileout.setf(std::ios::scientific);	// 指数表示
				fileout.width(15);					// 総桁数
				fileout.precision(14);				// 小数点桁数
			}
			if constexpr(E == CsvExpression::HEXFLOAT){
				// 16進数浮動小数点表記の場合
				fileout.setf(std::ios::scientific | std::ios::fixed);
			}
		}
		
		//! @brief 1行のCSV文字列データを浮動小数点std::arrayに読み込む関数
		//! @tparam	N	配列の長さ
		//! @param[in]	CsvLine	CSV文字列データ
		//! @param[out]	ValArray	数値配列
		template <size_t N>
		static void ParseCsvString(const std::string& CsvLine, std::array<double, N>& ValArray){
			size_t FirstIndex = 0;		// 抽出開始位置
			size_t LastIndex = 0;		// 抽出終了位置
			std::string ParseBuff = "";	// パース用バッファ
			double ValBuff = 0;			// 浮動小数点用バッファ
			size_t j = 0;			// 列番号インデックス
			while(LastIndex != std::string::npos){	// 「,」が出現しなくなるまでループ
				arcs_assert(j < N);				// 列番号溢れチェック
				LastIndex = CsvLine.find_first_of(",", FirstIndex);				// 次の「,」の位置を探す
				ParseBuff = CsvLine.substr(FirstIndex, LastIndex - FirstIndex);	// 「,」の手前まで文字列抽出
				ValBuff = std::stod(ParseBuff);	// 浮動小数点数値に変換
				ValArray.at(j) = ValBuff;		// 数値配列に書き込み
				FirstIndex = LastIndex + 1;		// 次のループでの「,」の探索開始位置
				++j;							// 列番号インデックス
			}
		}
		
		
};
}

#endif


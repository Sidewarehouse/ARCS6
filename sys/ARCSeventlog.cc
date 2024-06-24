//! @file ARCSeventlog.cc
//! @brief ARCS イベントログクラス
//!
//! ARCS用のイベントログクラス
//! 注意：リアルタイム空間ではこのライブラリが提供する関数は可能な限り使用しないこと。ジッタに影響するので。デバッグ時なら使用OK。
//!
//! @date 2024/06/24
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#include <fstream>
#include "ARCSeventlog.hh"
#include "ARCScommon.hh"
#include "ARCSscreen.hh"
#include "ARCSparams.hh"
#include "ConstParams.hh"

using namespace ARCS;

// 静的メンバ変数の実体
ARCSscreen* ARCSeventlog::ARCSscreenPtr = nullptr;	//!< ARCS画面ポインタ

//! @brief コンストラクタ
ARCSeventlog::ARCSeventlog(void){
	// イベントログファイルの準備とヘッダの書き出し
	std::ofstream EventLogFile(ARCSparams::EVENTLOG_NAME, std::ios::out | std::ios::trunc);
	EventLogFile << "ARCS EVENT LOG FILE" << std::endl;
	EventLogFile << "DATE: " << ARCScommon::GetNowTime();
	EventLogFile << "CTRLNAME: " << ConstParams::CTRLNAME << std::endl;
	EventLogFile << "ARCS_REVISION: " << ARCSparams::ARCS_REVISION << std::endl;
	EventLogFile << std::endl;
	EventLogFile << "CPU:TIME: FILE:LINE: MESSAGE" << std::endl;
	PassedLog();
}

//! @brief デストラクタ
ARCSeventlog::~ARCSeventlog(){
	PassedLog();
}

//! @brief ARCS画面ポインタの設定
//! @param[in]	ScrPtr	ARCS画面へのポインタ
void ARCSeventlog::SetScreenPtr(ARCSscreen* ScrPtr){
	ARCSscreenPtr = ScrPtr;
}

//! @brief イベントログを残す関数 (任意メッセージ記録版)
//! @param[in] str 記録したいメッセージ
//! @param[in] file 記録したいファイル名
//! @param[in] line 記録したい行番号
//! @param[in] cpu	CPUコア
//! @param[in] time	CPU時間
void ARCSeventlog::EventLog_from_macro(
	const std::string& str, const std::string& file, const int line, const int cpu, const clock_t time
){
	WriteEventLog(str, file, line, cpu, time);		// イベントログをファイルに書き出す
	if(ARCSscreenPtr != nullptr){
		// 画面が準備できていたら，イベントログデータを画面バッファに書き込む
		ARCSscreenPtr->WriteEventLogBuffer(file + " " + std::to_string(line) + ": " + str);
	}
}

//! @brief 変数用イベントログ
//! @param[in] u 表示する変数
//! @param[in] varname 変数名
//! @param[in] file 記録したいファイル名
//! @param[in] line 記録したい行番号
//! @param[in] cpu	CPUコア
//! @param[in] time	CPU時間
void ARCSeventlog::EventLogVar_from_macro(
	const double u, const std::string& varname, const std::string& file, const int line, const int cpu, const clock_t time
){
	std::string str = varname + " = " + ARCScommon::DoubleToString(u,"%g");
	EventLog_from_macro(str, file, line, cpu, time);
}

//! @brief 通過確認用ログを残す関数（ファイルと行番号のみ記録版）
//! @param[in] str 記録したいメッセージ
//! @param[in] file 記録したいファイル名
//! @param[in] line 記録したい行番号
//! @param[in] cpu	CPUコア
//! @param[in] time	CPU時間
void ARCSeventlog::PassedLog_from_macro(
	const std::string& file, const int line, const int cpu, const clock_t time
){
	const std::string str = "PASSED";			// 任意メッセージは"PASSED"
	EventLog_from_macro(str, file, line, cpu, time);	// イベントログをファイルに書き出す＆画面表示
}

//! @brief イベントログをファイルに書き出す
//! @param[in] str ログに書き残したい文字列
//! @param[in] file ファイル名
//! @param[in] line 行番号
void ARCSeventlog::WriteEventLog(const std::string& str, const std::string& file, const int line){
	std::ofstream EventLogFile(ARCSparams::EVENTLOG_NAME, std::ios::out | std::ios::app);
	EventLogFile << file << " " << line << ": " << str <<std::endl;
}

//! @brief イベントログをファイルに書き出す(CPU時間＆コアあり版)
//! @param[in] str ログに書き残したい文字列
//! @param[in] file ファイル名
//! @param[in] line 行番号
//! @param[in] cpu	CPUコア
//! @param[in] time	CPU時間
void ARCSeventlog::WriteEventLog(const std::string& str, const std::string& file, const int line, const int cpu, const clock_t time){
	std::ofstream EventLogFile(ARCSparams::EVENTLOG_NAME, std::ios::out | std::ios::app);
	EventLogFile << cpu << ":" << time << ": " << file << ":" << line << ": " << str <<std::endl;
}

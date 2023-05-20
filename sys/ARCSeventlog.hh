//! @file ARCSeventlog.hh
//! @brief ARCS イベントログクラス
//!
//! ARCS用のイベントログクラス
//! 注意：リアルタイム空間ではこのライブラリが提供する関数は可能な限り使用しないこと。ジッタに影響するので。デバッグ時なら使用OK。
//!
//! @date 2020/03/06
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef ARCSEVENTLOG
#define ARCSEVENTLOG

#include <sched.h>
#include <string>
#include <ctime>

// 前方宣言
namespace ARCS{
	class ARCSscreen;
}

// 関数呼び出し用マクロ
#define PassedLog() (ARCSeventlog::PassedLog_from_macro(__FILE__,__LINE__,sched_getcpu(),std::clock()))				//!< イベントログ用マクロ（ファイルと行番号のみ記録版）
#define EventLog(a) (ARCSeventlog::EventLog_from_macro(a,__FILE__,__LINE__,sched_getcpu(),std::clock()))			//!< イベントログ用マクロ (任意メッセージ記録版)
#define EventLogVar(a) (ARCSeventlog::EventLogVar_from_macro(a,#a,__FILE__,__LINE__,sched_getcpu(),std::clock()))	//!< 変数用イベントログマクロ  a : 表示する変数

namespace ARCS {
	//! @brief ARCS イベントログクラス
	class ARCSeventlog {
		public:
			ARCSeventlog();		//!< コンストラクタ
			~ARCSeventlog();	//!< デストラクタ
			void SetScreenPtr(ARCSscreen* ScrPtr);	//!< ARCS画面ポインタの設定
			
			// 下記の関数はマクロから呼ばれることを想定
			static void EventLog_from_macro(
				const std::string& str,
				const std::string& file, const int line, const int cpu, const clock_t time);	//!< イベントログを残す関数 (任意メッセージ記録版)
			static void EventLogVar_from_macro(
				const double u, const std::string& varname,
				const std::string& file, const int line, const int cpu, const clock_t time);	//!< 変数用イベントログ  u : 表示する変数, varname : 変数名
			static void PassedLog_from_macro(
				const std::string& file, const int line, const int cpu, const clock_t time);	//!< イベントログを残す関数（ファイルと行番号のみ記録版）
			
			static void WriteEventLog(const std::string& str, const std::string& file, const int line);		//!< イベントログをファイルに書き出す
			static void WriteEventLog(
				const std::string& str,
				const std::string& file, const int line, const int cpu, const clock_t time);	//!< イベントログをファイルに書き出す(CPU時間＆コアあり版)
			
		private:
			ARCSeventlog(const ARCSeventlog&) = delete;					//!< コピーコンストラクタ使用禁止
			const ARCSeventlog& operator=(const ARCSeventlog&) = delete;//!< 代入演算子使用禁止
			
			static ARCSscreen* ARCSscreenPtr;	//!< ARCS画面ポインタ
	};
}

#endif


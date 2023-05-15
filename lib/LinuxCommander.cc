//! @file LinuxCommander.cc
//! @brief Linuxコマンダー
//!
//! 任意のLinuxシェルコマンドを実行するクラス
//!
//! @date 2020/03/10
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include <cassert>
#include <stdio.h>
#include "LinuxCommander.hh"

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

using namespace ARCS;

//! @brief Linuxコマンドを実行する関数
//! @param[in]	CommandLine	Linuxコマンド
void LinuxCommander::Execute(std::string CommandLine){
	EventLog(CommandLine);
	FILE* fp = popen(CommandLine.c_str(), "w");
	arcs_assert(fp != nullptr);
	pclose(fp);
}

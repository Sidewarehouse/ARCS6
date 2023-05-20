//! @file LinuxCommander.hh
//! @brief Linuxコマンダー
//!
//! 任意のLinuxシェルコマンドを実行するクラス
//!
//! @date 2019/07/23
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef LINUXCOMMANDER
#define LINUXCOMMANDER

#include <string>

namespace ARCS {	// ARCS名前空間
	//! @brief Linuxコマンダー
	class LinuxCommander {
		public:
			static void Execute(std::string CommandLine);	//!< Linuxコマンドを実行する関数
			
		private:
			LinuxCommander() = delete;						//!< コンストラクタ使用禁止
			LinuxCommander(LinuxCommander&& right) = delete;//!< ムーブコンストラクタ使用禁止
			~LinuxCommander() = delete;						//!< デストラクタ使用禁止
			LinuxCommander(const LinuxCommander&) = delete;					//!< コピーコンストラクタ使用禁止
			const LinuxCommander& operator=(const LinuxCommander&) = delete;//!< 代入演算子使用禁止
	};
}

#endif


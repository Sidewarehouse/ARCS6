//! @file CPUSettings.hh
//! @brief CPU設定クラス
//!
//! CPUのコア，ポリシーなどを設定するクラス
//!
//! @date 2019/07/25
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef CPUSETTINGS
#define CPUSETTINGS

#include <pthread.h>

namespace ARCS {	// ARCS名前空間
	//! @brief CPU設定クラス
	class CPUSettings {
		public:
			static void SetCPUandPolicy(const pthread_t FuncName, const unsigned int CPUnum, const int Policy);	//!< スレッドを動作させるCPUコアとポリシーを設定する関数
			
		private:
			CPUSettings() = delete;						//!< コンストラクタ使用禁止
			CPUSettings(CPUSettings&& right) = delete;	//!< ムーブコンストラクタ使用禁止
			~CPUSettings() = delete;					//!< デストラクタ使用禁止
			CPUSettings(const CPUSettings&) = delete;					//!< コピーコンストラクタ使用禁止
			const CPUSettings& operator=(const CPUSettings&) = delete;	//!< 代入演算子使用禁止
	};
}

#endif


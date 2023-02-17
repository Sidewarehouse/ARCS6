//! @file ClassTemplate.hh
//! @brief クラステンプレート
//!
//! クラス(テンプレート版)を追加する場合は，このクラステンプレートを基に作ってネ。
//!
//! @date 20XX/XX/XX
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-20XX Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef CLASSTEMPLATE
#define CLASSTEMPLATE

#include <cassert>

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
//! @brief クラステンプレート
//! @tparam 
//template <>
class ClassTemplate {
	public:
		//! @brief コンストラクタ
		ClassTemplate()
			// :
		{
			
		}

		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		ClassTemplate(ClassTemplate&& r)
			// :
		{
			
		}

		//! @brief デストラクタ
		~ClassTemplate(){
			
		}
		
	private:
		ClassTemplate(const ClassTemplate&) = delete;					//!< コピーコンストラクタ使用禁止
		const ClassTemplate& operator=(const ClassTemplate&) = delete;	//!< 代入演算子使用禁止
};
}

#endif


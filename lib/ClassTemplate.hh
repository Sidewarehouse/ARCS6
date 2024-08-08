//! @file ClassTemplate.hh
//! @brief クラステンプレート
//!
//! クラス(テンプレート版)を追加する場合は，このクラステンプレートを基に作ってネ。
//!
//! @date 20XX/XX/XX
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-20XX Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

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
		ClassTemplate() noexcept
			// :
		{
			
		}

		//! @brief ムーブコンストラクタ
		//! @param[in]	r	右辺値
		ClassTemplate(ClassTemplate&& r) noexcept
			// :
		{
			
		}

		//! @brief ムーブ代入演算子
		//! @param[in]	r	右辺値
		ClassTemplate& operator=(ClassTemplate&& r) noexcept {
			return *this;
		}
		
		//! @brief デストラクタ
		~ClassTemplate() noexcept {
			
		}
		
	private:
		ClassTemplate(const ClassTemplate&) = delete;					//!< コピーコンストラクタ使用禁止
		const ClassTemplate& operator=(const ClassTemplate&) = delete;	//!< コピー代入演算子使用禁止
};
}

#endif


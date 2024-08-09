//! @file ClassBase.cc
//! @brief クラスベースコード
//!
//! クラス(非テンプレート版)を追加する場合は，このクラスベースコードを基に作ってネ。
//!
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-20XX Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#include <cassert>
#include "ClassBase.hh"

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

//! @brief コンストラクタ
ClassBase::ClassBase() noexcept
	// :
{
	
}

//! @brief ムーブコンストラクタ
//! @param[in]	r	演算子右側
ClassBase::ClassBase(ClassBase&& r) noexcept
	// :
{
	
}

//! @brief ムーブ代入演算子
//! @param[in]	r	演算子右側
ClassBase& ClassBase::operator=(ClassBase&& r) noexcept {
	return *this;
}
		
//! @brief デストラクタ
ClassBase::~ClassBase() noexcept {
	
}


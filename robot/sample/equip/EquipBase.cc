//! @file EquipBase.cc
//! @brief クラスベースコード(Equipディレクトリ版)
//!
//! クラス(非テンプレート版)を追加する場合は，このクラスベースコードを基に作ってネ。
//!
//! @date 20XX/XX/XX
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-20XX Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include <cassert>
#include "EquipBase.hh"

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
EquipBase::EquipBase()
	// :
{
	
}

//! @brief ムーブコンストラクタ
//! @param[in]	r	右辺値
EquipBase::EquipBase(EquipBase&& r)
	// :
{
	
}

//! @brief デストラクタ
EquipBase::~EquipBase(){
	
}


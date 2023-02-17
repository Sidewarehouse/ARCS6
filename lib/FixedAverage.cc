//! @file FixedAverage.cc
//! @brief FixedAverageクラス
//!
//! 時刻ゼロから現在時刻までの信号の平均を計算する。
//!
//! @date 2019/07/29
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include "FixedAverage.hh"

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
FixedAverage::FixedAverage()
	: count(0), ave(0)
{
	PassedLog();
}

//! @brief ムーブコンストラクタ
//! @param[in]	r	右辺値
FixedAverage::FixedAverage(FixedAverage&& r)
	: count(r.count), ave(r.ave)
{
	
}

//! @brief デストラクタ
FixedAverage::~FixedAverage(){
	PassedLog();
}

//! @brief ひたすら平均値を計算し続ける関数
//! @param[in]	in	入力信号
//! @return	平均値
double FixedAverage::GetSignal(double in){
	count++;	// 呼ばれた回数をカウント
	ave += in;	// ひたすら加算
	return ave/(double)count;	// 平均値を返す
}

//! @brief 平均値とカウント数をリセットする関数
void FixedAverage::Reset(void){
	count = 0;
	ave = 0;
}


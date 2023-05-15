//! @file TriangleWave.cc
//! @brief 三角波発振器
//! @date 2022/03/09
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2022 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include <cassert>
#include <cmath>
#include "TriangleWave.hh"

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

//! @brief 三角波発振器
//! @param[in]	freq	[Hz] 周波数
//! @param[in]	time	[s]  時刻
//! @return	振幅±1の三角波出力
double ARCS::TriangleWave(const double freq, const double time){
	const double Tp = 1.0/freq;
	const double a = 2.0/Tp;
	double t, y;
	
	t = fmod(time,Tp);	// 時刻を0～Tpの時間範囲に収める
	
	// 三角波の生成
	if(0 <= t && t < Tp/2.0){
		y = a*t;		// 正の傾き
	}else{
		y = -a*t + 2.0;	// 負の傾き
	}
	
	return 2.0*y - 1.0;	// ±1になるようにする
}

//! @brief 三角波発振器（指定開始時刻まではゼロ出力）
//! @param[in]	freq		[Hz] 周波数
//! @param[in]	time		[s]  時刻
//! @param[in]	starttime	[s]  開始時刻
//! @return	振幅±1の三角波出力
double ARCS::TriangleWave(const double freq, const double time, const double starttime){
	double ret = 0;
	if(time < starttime){
		ret = 0;
	}else{
		ret = TriangleWave(freq, time - starttime);
	}
	return ret;
}

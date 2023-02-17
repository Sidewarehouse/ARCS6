//! @file ChirpGenerator.cc
//! @brief チャープ信号生成器クラス
//!
//! チャープ信号生成器
//!
//! @date 2019/07/29
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include <cmath>
#include "ChirpGenerator.hh"

using namespace ARCS;

//! @brief コンストラクタ
//! @param[in]	freq_min	開始周波数
//! @param[in]	freq_max	終了周波数
//! @param[in]	time_max	終了時刻
ChirpGenerator::ChirpGenerator(double freq_min, double freq_max, double time_max)
	: fmin(freq_min), fmax(freq_max), tmax(time_max), A(0)
{
	A = (fmax-fmin)/tmax;
}

//! @brief デストラクタ
ChirpGenerator::~ChirpGenerator(){
	
}

//! @brief チャープ信号の取得
//! @param[in]	time	現在時刻
double ChirpGenerator::GetSignal(double time) const {
	return sin( 2.0*M_PI*(A*time*time/2.0 + fmin*time) );
}

//! @brief 現在の周波数 [Hz] の取得
//! @param[in]	time	現在時刻
//! @return	周波数 [Hz]
double ChirpGenerator::GetFrequency(double time) const {
	return A*time + fmin;
}


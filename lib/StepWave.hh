//! @file StepWave.cc
//! @brief ステップ波形生成器
//!
//! ステップ波形を生成する。
//!
//! @date 2020/04/10
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef STEPWAVE
#define STEPWAVE

namespace ARCS {	// ARCS名前空間
	double StepWave(const double Tstp, const double t);	//!< 1回だけステップ波形を出力する関数
}

#endif


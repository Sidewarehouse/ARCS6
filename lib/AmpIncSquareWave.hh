//! @file AmpIncSquareWave.hh
//! @brief 振幅増加方形波生成器
//!
//! 時間に従って振幅が増加していく方形波を生成する関数
//!
//! @date 2021/12/09
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2021 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef AMPINCSQUAREWAVE
#define AMPINCSQUAREWAVE

namespace ARCS {	// ARCS名前空間
	double AmpIncSquareWave(const double Time, const double Freq, const double Ystp, const double Nstp);
}

#endif


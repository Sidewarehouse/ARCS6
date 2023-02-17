//! @file TriangleWave.cc
//! @brief 三角波発振器
//! @date 2022/03/09
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2022 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef TRIANGLEWAVE
#define TRIANGLEWAVE

namespace ARCS {	// ARCS名前空間
	double TriangleWave(const double freq, const double time);							//!< 三角波発振器
	double TriangleWave(const double freq, const double time, const double starttime);	//!< 三角波発振器（指定開始時刻まではゼロ出力）
}

#endif

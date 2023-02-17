//! @file TwoStairsWave.cc
//! @brief 2段階段波発生器
//!
//! @date 2022/11/14
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2022 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef TWOSTAIRSWAVE
#define TWOSTAIRSWAVE

namespace ARCS {	// ARCS名前空間
double TwoStairsWave(const double freq, const double time);							//!< 2段階段波を出力する関数
double TwoStairsWave(const double freq, const double time, const double starttime);	//!< 2段階段波を出力する関数（開始時刻指定版）

}

#endif

//! @file Sigmoid.hh
//! @brief シグモイド関数
//! y = 1/(1 + exp(-a*u))
//! @date 2017/06/09
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2017 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef SIGMOID
#define SIGMOID

namespace ARCS {	// ARCS名前空間
	double Sigmoid(double u, double a);	//!< シグモイド関数
}

#endif


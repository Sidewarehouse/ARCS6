//! @file MotorFrameTransform.cc
//! @brief モータ座標変換関数群
//!
//! モータの座標変換に関連する関数群
//!
//! @date 2019/09/18
//! @author Yusuke Asai & Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yusuke Asai & Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef DQTRANSFORM
#define DQTRANSFORM

namespace ARCS {
	std::tuple<double, double> UVWtoAB(const double u, const double v, const double w);						//!< UVW相からαβ軸に変換する関数(タプル版)
	std::tuple<double, double> UVWtoDQ(const double u, const double v, const double w, const double theta);	//!< UVW相からdq軸に変換する関数(タプル版)
	std::tuple<double, double, double> DQtoUVW(const double d, const double q, const double theta);			//!< dq軸からUVW相に変換する関数(タプル版)
	void dq_uvw(double d, double q, double theta, double *u, double *v, double *w);		// dq→uvw変換
	void dq_uvw_inv(double u, double v, double w, double theta, double *d, double *q);	// 逆dq→uvw変換
	void dq_uv(double d, double q, double theta, double *u, double *v);	// dq-uv 変換行列 位置の単位に注意　degではなく『rad』
}

#endif


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

#include <cmath>
#include <tuple>
#include "MotorFrameTransform.hh"

//! @brief UVW相からαβ軸に変換する関数(タプル版)
//! @param[in]	u	U相
//! @param[in]	v	V相
//! @param[in]	w	W相
//! @return	αβ軸 (α, β)のタプル
std::tuple<double, double> ARCS::UVWtoAB(const double u, const double v, const double w){
	double a = sqrt(2.0/3.0)*( 1.0*u -       1.0/2.0*v -       1.0/2.0*w );
	double b = sqrt(2.0/3.0)*(         sqrt(3.0)/2.0*v - sqrt(3.0/2.0)*w );
	return std::forward_as_tuple(a, b);
}


//! @brief UVW相からdq軸に変換する関数(タプル版)
//! @param[in]	u	U相
//! @param[in]	v	V相
//! @param[in]	w	W相
//! @return	dq軸 (d, q)のタプル
std::tuple<double, double> ARCS::UVWtoDQ(const double u, const double v, const double w, const double theta){
	double d = sqrt(2.0/3.0)*(   u*cos(theta) + v*cos(theta - 2.0/3.0*M_PI) + w*cos(theta + 2.0/3.0*M_PI) );
	double q = sqrt(2.0/3.0)*( - u*sin(theta) - v*sin(theta - 2.0/3.0*M_PI) - w*sin(theta + 2.0/3.0*M_PI) );
	return std::forward_as_tuple(d, q);		// タプルで返す
}

//! @brief dq軸からUVW相に変換する関数(タプル版)
//! @param[in]	d	d軸
//! @param[in]	q	q軸
//! @return	UVW相 (u, v, w)のタプル
std::tuple<double, double, double> ARCS::DQtoUVW(const double d, const double q, const double theta){
	double u = sqrt(2.0/3.0)*( d*cos(theta               ) - q*sin(theta               ) );	// U相の計算
	double v = sqrt(2.0/3.0)*( d*cos(theta - 2.0/3.0*M_PI) - q*sin(theta - 2.0/3.0*M_PI) );	// V相の計算
	double w = sqrt(2.0/3.0)*( d*cos(theta + 2.0/3.0*M_PI) - q*sin(theta + 2.0/3.0*M_PI) );	// W相の計算
	return std::forward_as_tuple(u, v, w);	// タプルで返す
}

void ARCS::dq_uvw(double d, double q, double theta, double *u, double *v, double *w){
	// dq-uvw 変換行列 位置の単位に注意　degではなく『rad』
	*u=sqrt(2.0/3.0)*(d*cos(theta)-q*sin(theta));
	*v=sqrt(2.0/3.0)*(d*cos(theta-(2.0/3.0*M_PI))-q*sin(theta-(2.0/3.0*M_PI)));
	*w=-*u-*v;	
}

void ARCS::dq_uvw_inv(double u, double v, double w, double theta, double *d, double *q){
	// 逆行列
	*d=sqrt(2.0/3.0)*(u*cos(theta)+v*cos(theta-2.0/3.0*M_PI)+w*cos(theta+2.0/3.0*M_PI));
	*q=-sqrt(2.0/3.0)*(u*sin(theta)+v*sin(theta-2.0/3.0*M_PI)+w*sin(theta+2.0/3.0*M_PI));
}

void ARCS::dq_uv(double d, double q, double theta, double *u, double *v){
	// dq-uv 変換行列 位置の単位に注意　degではなく『rad』
	*u=sqrt(2.0/3.0)*(d*cos(theta)-q*sin(theta));
	*v=sqrt(2.0/3.0)*(d*cos(theta-(2.0/3.0*M_PI))-q*sin(theta-(2.0/3.0*M_PI)));
}


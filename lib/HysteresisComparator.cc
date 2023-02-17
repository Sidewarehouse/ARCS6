//! @file HysteresisComparator.cc
//! @brief ヒステリシス比較器
//! @date 2019/02/22
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include "HysteresisComparator.hh"

using namespace ARCS;

//! @brief コンストラクタ
//! @param[in]	Left	ヒステリシス特性の左側の比較値
//! @param[in]	Right	ヒステリシス特性の右側の比較値
HysteresisComparator::HysteresisComparator(double Left, double Right)
	: R(Right), L(Left), yZ1(0)
{
	
}

//! @brief デストラクタ
HysteresisComparator::~HysteresisComparator(){
	
}

//! @brief 出力信号の取得
//! @param[in]	u	入力信号
//! @return ヒステリシス比較出力
double HysteresisComparator::GetSignal(double u){
	double y = 0;
	
	if(yZ1 <= 0 && R < u){		// ヒステリシス特性の右側の比較
		y = 1;
	}else if(0 < yZ1 && u < L){	// ヒステリシス特性の左側の比較
		y = -1;
	}else{
		y = yZ1;	// 変化なし
    }
	
	// 状態変数の更新
	yZ1 = y;
	
	return y;
}


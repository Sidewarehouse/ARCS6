//! @file TwoStepWave.cc
//! @brief 2段ステップ信号生成器
//!
//! 2段ステップの信号を出力する。
//!
//! @date 2019/02/26
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include "TwoStepWave.hh"

using namespace ARCS;

TwoStepWave::TwoStepWave(
	double Stp1Hgt, double Stp2Hgt,
	double Stp1StaT, double Stp2StaT, double Stp2EndT, double Stp1EndT)
	// コンストラクタ  ステップ1の高さ, ステップ2の高さ，ステップ1の開始時刻, ステップ2の開始時刻, ステップ2の終了時刻, ステップ1の終了時刻
	: Height1(Stp1Hgt), Height2(Stp2Hgt), T1sta(Stp1StaT), T2sta(Stp2StaT), T2end(Stp2EndT), T1end(Stp1EndT)
{
}

double TwoStepWave::GetSignal(double t) const {
	// 2段ステップの信号を出力する t : 時刻
	double ret = 0;
	if(t < T1sta){
		ret = 0;
	}else if(T1sta <= t && t < T2sta){
		ret = Height1;
	}else if(T2sta <= t && t < T2end){
		ret = Height2;
	}else if(T2end <= t && t < T1end){
		ret = Height1;
	}else if(T1end <= t){
		ret = 0;
	}
	return ret;
}

TwoStepWave::~TwoStepWave(){
	// デストラクタ
}

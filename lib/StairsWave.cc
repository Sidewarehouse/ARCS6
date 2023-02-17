//! @file StairsWave.cc
//! @brief 階段波発生器
//! @date 2017/04/21
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2017 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include "StairsWave.hh"
#include <math.h>

//! @brief 階段波発生器
//! @param[in]	Time	[s] 現在時刻
//! @param[in]	Tini	[s] 階段の初期時刻
//! @param[in]	Ystp	[-] 階段1つ分の高さ
//! @param[in]	Tstp	[s] 階段1つ分の時間長さ
//! @param[in]	Nstp	[-] 階段の数
//! @return		階段波出力
double ARCS::StairsWave(double Time, double Tini, double Ystp, double Tstp, double Nstp){
	int n = (int)((Time-Tini)/Tstp);	// 現在の階段の段数を計算
	double y;
	
	if(n <= Nstp){
		y = Ystp*(double)n;
	}else{
		y = 0;
	}
	
	return y;
}


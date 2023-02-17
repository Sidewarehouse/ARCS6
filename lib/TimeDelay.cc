//! @file TimeDelay.cc
//! @brief 遅延器クラス
//!
//! 遅延器 G(s)=e^(-s*T*n) -> G(z)=z^(-n)
//! (但し T は制御周期で，制御周期単位でしか遅延時間を設定できないことに注意)
//!
//! @date 2019/02/26
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2019 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include "TimeDelay.hh"

using namespace ARCS;

TimeDelay::TimeDelay(const long MaxDelay)
	// コンストラクタ MaxDelay；最大遅延時間
	: dmem_max(MaxDelay),
	  num(0), Wcount(0), Rcount(0), dmem(0)
{
	dmem=new double [MaxDelay];
	ClearDelayMemory();	// 遅延メモリのゼロクリア
}

TimeDelay::~TimeDelay(){
	// デストラクタ
	delete [] dmem;
}

double TimeDelay::GetSignal(const double u){
	// 出力信号の取得 u；入力信号
	// 下記アルゴリズムはリングバッファと同一
	double y;
	
	// カウンタを進める
	Rcount++;
	Wcount=Rcount+num;
	
	// カウンタ値の場合分け
	if(dmem_max<=Wcount)Wcount=Wcount-dmem_max;		// 最大値を超える場合
	if(dmem_max<=Rcount)Rcount=Rcount-dmem_max;		// 最大値を超える場合
	
	// 要素番号リミッタ(念のため)
	if(Wcount<0)Wcount=0;
	if(dmem_max<=Wcount)Wcount=dmem_max;
	if(Rcount<0)Rcount=0;
	if(dmem_max<=Rcount)Rcount=dmem_max;

	dmem[Wcount]=u;	// 遅延メモリへの書き込み
	y=dmem[Rcount];	// 遅延メモリから読み出し

	return y; 
}

void TimeDelay::SetDelayTime(const long DelayTime){
	// 遅延時間の設定 DelayTime；遅延時間 (最大遅延時間を越えないこと)
	if(DelayTime<dmem_max){
		num=DelayTime;	// 設定が最大遅延時間を越えてないとき
	}else{
		num=dmem_max-1;	// 設定が最大遅延時間を越えているとき
	}
}

void TimeDelay::ClearDelayMemory(void){
	// 遅延メモリのゼロクリア
	for(long i=0;i<dmem_max;i++)dmem[i]=0;
}





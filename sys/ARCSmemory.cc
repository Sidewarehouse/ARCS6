//! @file ARCSmemory.cc
//! @brief データメモリクラス
//!
//! 実験データ保存用のデータメモリクラス。CSVファイルへの出力も行う。
//!
//! @date 2024/05/06
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2024 Yokokura, Yuki
// MIT License. For details, see the LICENSE file.

#include <cassert>
#include "ARCSmemory.hh"
#include "CsvManipulator.hh"

// ARCS組込み用マクロ
#ifdef ARCS_IN
	// ARCSに組み込まれる場合
	#include "ARCSassert.hh"
	#include "ARCSeventlog.hh"
#else
	// ARCSに組み込まれない場合
	#define arcs_assert(a) (assert(a))
	#define PassedLog()
	#define EventLog(a)
	#define EventLogVar(a)
#endif

using namespace ARCS;

//! @brief コンストラクタ
ARCSmemory::ARCSmemory()
	: SaveBuffer(nullptr),
	  Tindex(0),
	  Nindex(0),
	  Tperiod(0),
	  Time(0)
{
	PassedLog();
	
	// 開始時刻と終了時刻が入れ替わってないかのチェック
	static_assert(ConstParams::DATA_START < ConstParams::DATA_END);
	
	// データバッファのメモリ確保とゼロ埋め
	SaveBuffer = std::make_unique< std::array<std::array<double, ConstParams::DATA_NUM >, ARCSmemory::ELEMENT_NUM> >();
	for(size_t j = 0; j < ARCSmemory::ELEMENT_NUM; ++j){
		SaveBuffer->at(j).fill(0);
	}
	
	PassedLog();
}

//! @brief デストラクタ
ARCSmemory::~ARCSmemory(){
	PassedLog();
}

//! @brief リセットする関数
void ARCSmemory::Reset(void){
	Tindex = 0;	// 時間用カウンタをゼロに戻す
}

//! @brief CSVファイルを書き出す関数
void ARCSmemory::WriteCsvFile(void){
	CsvManipulator::SaveFile(std::move(SaveBuffer), ConstParams::DATA_NAME, ConstParams::DATA_NUM, Tindex);
}

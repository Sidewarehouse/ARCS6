//! @file SA2-RasPi2.cc
//! @brief ServoAmplifier2-RaspberryPi2インターフェースクラス
//!
//! ServoAmplifier2を制御するためのラズベリーパイ2の色々な入出力機能を提供します。
//! (即値マジックナンバーが埋め込まれているが，そのうちキレイにする予定)
//!
//! @date 2020/04/10
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2020 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include <cassert>
#include <cmath>
#include "SA2-RasPi2.hh"
#include "Limiter.hh"

// ARCS組込み用マクロ
#ifdef ARCS_IN
	// ARCSに組み込まれる場合
	#include "ARCSassert.hh"
	#include "ARCSeventlog.hh"
	#include "ARCSprint.hh"
#else
	// ARCSに組み込まれない場合
	#define arcs_assert(a) (assert(a))
	#define PassedLog()
	#define EventLog(a)
	#define EventLogVar(a)
#endif

using namespace ARCS;

//! @brief コンストラクタ
SA2RasPi2::SA2RasPi2()
	: gpio(),
	  DATA_PWMU(DutyToPwmData(0.5)), DATA_PWMV(DutyToPwmData(0.5)), DATA_PWMW(DutyToPwmData(0.5)),
	  DATA_DAC1(VoltToDacData(0)), DATA_DAC2(VoltToDacData(0)),
	  DATA_IU(0), DATA_IV(0), DATA_IW(0), DATA_ABS(0), DATA_INC(0),
	  RawCountPrev(0), UpperCount(0), AbsMultCount(0),
	  IncRawCountPrev(0), IncUpperCount(0)
{
	SetBusModeIn();			// バス衝突を避けるためにデフォルトでデータバスを入力モードにしておく
	ClearReadSelect();		// RD_SEL信号を念のためのバス衝突防止のために停止
}

//! @brief ムーブコンストラクタ
//! @param[in]	r	右辺値
SA2RasPi2::SA2RasPi2(SA2RasPi2&& r)
	: gpio(std::move(r.gpio)),
	  DATA_PWMU(DutyToPwmData(0.5)), DATA_PWMV(DutyToPwmData(0.5)), DATA_PWMW(DutyToPwmData(0.5)),
	  DATA_DAC1(VoltToDacData(0)), DATA_DAC2(VoltToDacData(0)),
	  DATA_IU(0), DATA_IV(0), DATA_IW(0), DATA_ABS(0), DATA_INC(0),
	  RawCountPrev(0), UpperCount(0), AbsMultCount(0),
	  IncRawCountPrev(0), IncUpperCount(0)
{
	
}

//! @brief デストラクタ
SA2RasPi2::~SA2RasPi2(){
	SetBusModeIn();			// バス衝突を避けるために最後も念のためにデータバスを入力モードにしておく
	ClearReadSelect();		// RD_SEL信号を念のためのバス衝突防止のために停止
}

//! @brief 開始処理
void SA2RasPi2::StartProcess(void){
	WaitPwmSyncHigh();				// PWM同期信号の立上りが来るまで待機
	SendConversionStart();			// 変換開始信号を出力
	//WaitEndOfConversionIuvw();	// 電流センサのAD変換完了を待機 (これを入れると何かの拍子にAD7899がEOCモードになりデッドロックする)
	DATA_IU = ReadDataBus(ADDR_IU);	// U相電流の読み込み
	DATA_IV = ReadDataBus(ADDR_IV);	// V相電流の読み込み
	DATA_IW = ReadDataBus(ADDR_IW);	// W相電流の読み込み
	WaitEndOfConversionInc();			// インクリメンタルエンコーダの変換処理を待機
	DATA_INC = ReadDataBus(ADDR_INCL);	// インクリメンタルエンコーダ下位の読み込み
	WaitEndOfConversionAbs();					// 絶対値エンコーダのSP変換を待機
	DATA_ABS  =  ReadDataBus(ADDR_ABSL);		// 絶対値エンコーダ下位の読み込み
	DATA_ABS |= (ReadDataBus(ADDR_ABSH) << 14);	// 絶対値エンコーダ上位の読み込みと下位との結合
	AbsMultCount = GetMultiTurnCount(DATA_ABS);	// 絶対値多回転カウント値の取得
}

//! @brief 終了処理
void SA2RasPi2::EndProcess(void){
	WaitPwmSyncLow();					// PWM同期信号の立下りが来るまで待機
	SetBusModeOut();					// データバスを出力モードに変更
	WriteDataBus(ADDR_PWMU, DATA_PWMU);	// PWM U相へのデータ書き込み
	WriteDataBus(ADDR_PWMV, DATA_PWMV);	// PWM V相へのデータ書き込み
	WriteDataBus(ADDR_PWMW, DATA_PWMW);	// PWM W相へのデータ書き込み
	WriteDataBus(ADDR_DAC1, DATA_DAC1);	// DAC1へのデータ書き込み
	WriteDataBus(ADDR_DAC2, DATA_DAC2);	// DAC2へのデータ書き込み
	SetBusModeIn();						// データバスを入力モードに変更
}

//! @brief PWMの電圧指令を設定する関数
//! @param[in]	Vu	U相 電圧指令 [V]
//! @param[in]	Vv	V相 電圧指令 [V]
//! @param[in]	Vw	W相 電圧指令 [V]
//! @param[in]	Vdc	直流リンク電圧 [V]
void SA2RasPi2::SetVoltage(double Vu, double Vv, double Vw, double Vdc){
	SetDutyRatioPWM(Vu/Vdc + 0.5, Vv/Vdc + 0.5, Vw/Vdc + 0.5);
}

//! @brief PWMのデューティ比を設定する関数
//! @param[in]	DutyU	PWM U相デューティ比 [-]
//! @param[in]	DutyV	PWM V相デューティ比 [-]
//! @param[in]	DutyW	PWM W相デューティ比 [-]
void SA2RasPi2::SetDutyRatioPWM(double DutyU, double DutyV, double DutyW){
	SetDutyRatioPWMU(DutyU);
	SetDutyRatioPWMV(DutyV);
	SetDutyRatioPWMW(DutyW);
}

//! @brief 単相Hブリッジインバータとして使う場合のデューティ比を設定する関数
//! @param[in]	Duty	PWMデューティ比 [-]
void SA2RasPi2::SetSinglePhaseDutyRatio(double Duty){
	SetDutyRatioPWMU(0);			// U相はそのまま
	SetDutyRatioPWMV(Duty);			// V相はゼロ
	SetDutyRatioPWMW(1.0 - Duty);	// W相は反対向き
}

//! @brief PWM U相のデューティ比を設定する関数
//! @param[in]	Duty	PWMデューティ比 [-]
void SA2RasPi2::SetDutyRatioPWMU(double Duty){
	DATA_PWMU = DutyToPwmData(Duty);	// PWM U相 のデータバッファに書き込み
}

//! @brief PWM V相のデューティ比を設定する関数
//! @param[in]	Duty	PWMデューティ比 [-]
void SA2RasPi2::SetDutyRatioPWMV(double Duty){
	DATA_PWMV = DutyToPwmData(Duty);	// PWM V相 のデータバッファに書き込み
}

//! @brief PWM W相のデューティ比を設定する関数
//! @param[in]	Duty	PWMデューティ比 [-]
void SA2RasPi2::SetDutyRatioPWMW(double Duty){
	DATA_PWMW = DutyToPwmData(Duty);	// PWM W相 のデータバッファに書き込み
}

//! @brief D/Aコンバータの電圧を設定する関数
//! @param[in]	Vdac1	DAC CH1出力電圧 [V]
//! @param[in]	Vdac2	DAC CH2出力電圧 [V]
void SA2RasPi2::SetVoltageDAC(double Vdac1, double Vdac2){
	SetVoltageDAC1(Vdac1);
	SetVoltageDAC2(Vdac2);
}

//! @brief D/AコンバータCH1の電圧を設定する関数
//! @param[in]	Vdac	DAC出力電圧 [V]
void SA2RasPi2::SetVoltageDAC1(double Vdac){
	DATA_DAC1 = VoltToDacData(Vdac);	// DAC CH1 のデータバッファに書き込み
}

//! @brief D/AコンバータCH2の電圧を設定する関数
//! @param[in]	Vdac	DAC出力電圧 [V]
void SA2RasPi2::SetVoltageDAC2(double Vdac){
	DATA_DAC2 = VoltToDacData(Vdac);	// DAC CH2 のデータバッファに書き込み
}

//! @brief 相電流を取得する関数
//! @param[out]	Iu	U相電流 [A]
//! @param[out]	Iv	V相電流 [A]
//! @param[out]	Iw	W相電流 [A]
void SA2RasPi2::GetPhaseCurrent(double& Iu, double& Iv, double& Iw){
	Iu = AdcDataToAmpere(DATA_IU);	// [A]に換算して出力
	Iv = AdcDataToAmpere(DATA_IV);	// [A]に換算して出力
	Iw = AdcDataToAmpere(DATA_IW);	// [A]に換算して出力
}

//! @brief 単相電流を取得する関数
//! @param[out]	I	単相電流 [A]
void SA2RasPi2::GetSinglePhaseCurrent(double& I){
	I = (AdcDataToAmpere(DATA_IV) - AdcDataToAmpere(DATA_IW))*0.5;	// [A]に換算してVW相の平均取って出力
}

//! @brief モータ機械角と電気角を取得する関数
//! @param[out]	theta_rm	モータ機械角 [rad]
//! @param[out]	theta_re	モータ電気角 [rad]
void SA2RasPi2::GetPosition(double& theta_rm, double& theta_re){
	theta_rm = AbsDataToRadian(AbsMultCount);	// 機械角に換算して出力
	theta_re = ConvElectAngle(AbsMultCount);	// 電気角に換算して出力
}

//! @brief インクリメンタルエンコーダのバイナリデータを返す関数
//! @return	エンコーダデータ
uint32_t SA2RasPi2::GetIncEncoder(void){
	return INC_MASK & DATA_INC;
}

//! @brief インクリメンタルエンコーダの多回転バイナリデータを返す関数
//! @return 多回転エンコーダデータ
int32_t SA2RasPi2::GetIncMultiTurnEncoder(void){
	return GetIncMultiTurnCount(GetIncEncoder());
}

//! @brief インクリメンタルエンコーダの位置応答を返す関数
//! @return 位置応答 [rad]
double SA2RasPi2::GetIncPosition(void){
	return ((double)GetIncMultiTurnEncoder())*2.0*M_PI/4096.0;	// 12bitエンコーダをradに変換
}

//! @brief PWM_SYNC が High になるまでブロッキング待機する関数
void SA2RasPi2::WaitPwmSyncHigh(void){
	while(1){
		if(gpio.BitGet(27) == true) break;
	}
}

//! @brief PWM_SYNC が Low になるまでブロッキング待機する関数
void SA2RasPi2::WaitPwmSyncLow(void){
	while(1){
		if(gpio.BitGet(27) == false) break;
	}
}

//! @brief 電流センサのAD変換が完了するまでブロッキング待機する関数
void SA2RasPi2::WaitEndOfConversionIuvw(void){
	while(1){
		if(gpio.BitGet(23) == false) break;
	}
}

//! @brief 絶対値エンコーダのSP変換が完了するまでブロッキング待機する関数
void SA2RasPi2::WaitEndOfConversionAbs(void){
	while(1){
		if(gpio.BitGet(24) == false) break;
	}
}

//! @brief インクリメンタルエンコーダの変換処理が完了するまでブロッキング待機する関数
void SA2RasPi2::WaitEndOfConversionInc(void){
	while(1){
		if(gpio.BitGet(22) == false) break;
	}
}

//! @brief CONVST信号を出力する関数
void SA2RasPi2::SendConversionStart(void){
	gpio.BitSet(21);
	WaitUsingNop(100);
	gpio.BitClear(21);
}

//! @brief WR_SEL信号を出力する関数
void SA2RasPi2::SendWriteSelect(void){
	WaitUsingNop(5000);
	gpio.BitSet(20);
	WaitUsingNop(100);
	gpio.BitClear(20);
	WaitUsingNop(1000);
}

//! @brief RD_SEL信号を出力する関数
void SA2RasPi2::SetReadSelect(void){
	WaitUsingNop(100);
	gpio.BitSet(19);
	WaitUsingNop(100);
}

//! @brief RD_SEL信号を停止する関数
void SA2RasPi2::ClearReadSelect(void){
	WaitUsingNop(100);
	gpio.BitClear(19);
	WaitUsingNop(100);
}

//! @brief データバスを入力モードに設定する関数
void SA2RasPi2::SetBusModeIn(void){
	//                      GPIO--999888777666555444333222111000 (000 = 入力，001 = 出力)
	gpio.SetConfig_GPIO9to0(  0b00000000000000000000000000000000);	// GPIO9(D9)～2(D2)を出力に設定
	gpio.SetConfig_GPIO19to10(0b00001001001001000000000000000000);	// GPIO19(RD_SEL),18(A2),17(A1),16(A0)を出力, 13(D13)～10(D10)を出力に設定
	gpio.SetConfig_GPIO29to20(0b00000000000000000000000000001001);	// GPIO26(D1),25(D0)を出力, 21(CONVST),20(WR_SEL)を出力に設定
}

//! @brief データバスを出力モードに設定する関数
void SA2RasPi2::SetBusModeOut(void){
	//                      GPIO--999888777666555444333222111000 (000 = 入力，001 = 出力)
	gpio.SetConfig_GPIO9to0(  0b00001001001001001001001001000000);	// GPIO9(D9)～2(D2)を出力に設定
	gpio.SetConfig_GPIO19to10(0b00001001001001000000001001001001);	// GPIO19(RD_SEL),18(A2),17(A1),16(A0)を出力, 13(D13)～10(D10)を出力に設定
	gpio.SetConfig_GPIO29to20(0b00000000000001001000000000001001);	// GPIO26(D1),25(D0)を出力, 21(CONVST),20(WR_SEL)を出力に設定
}

//! @brief NOP命令を使って待機する関数
//! @param[in]	N	NOP命令を実行する数
void SA2RasPi2::WaitUsingNop(const unsigned int N) const{
	for(unsigned int i = 0; i < N; i++) asm("nop");
}

//! @brief 指定したアドレスでデータバスに書き込む関数
//! @param[in]	Address	アドレス
//! @param[in]	Data	データ
void SA2RasPi2::WriteDataBus(uint32_t Address, uint32_t Data){
	uint32_t BusBuff    = 0;					// バスバッファ
	uint32_t NowState = gpio.GetBitDataLo();	// 現在のGPIOのビット列を取得
	// GPIO番号                       3         2         1         0
	//                               10987654321098765432109876543210
	BusBuff  |= ((Address << 16) & 0b00000000000001110000000000000000);	// 書き込んではいけないビットをANDでマスクして，アドレス部分のビットをORで書き込む
	BusBuff  |=            (Data & 0b00000000000000000011111111111100);	// 書き込んではいけないビットをANDでマスクして，データ部分のビットをORで書き込む
	BusBuff  |=    ((Data << 25) & 0b00000110000000000000000000000000);	// 書き込んではいけないビットをANDでマスクして，データ部分のビットをORで書き込む
	NowState &=                    0b11111001111110001100000000000011;	// 現在のGPIOのビット列の触っていい部分をANDでクリアしてから，
	NowState |= BusBuff;			// 現在のGPIOのビット列にORで書き込んで，
	gpio.SetBitDataLo(NowState);	// GPIOからビット列出力
	SendWriteSelect();				// WR_SEL信号の出力
}

//! @brief 指定したアドレスでデータバスから読み込む関数
//! @param[in]	Address	アドレス
//! @return	データバスのビット列
uint32_t SA2RasPi2::ReadDataBus(uint32_t Address){
	// アドレスバスから出力
	uint32_t BusBuff    = 0;					// バスバッファ
	uint32_t NowState = gpio.GetBitDataLo();	// 現在のGPIOのビット列を取得
	// GPIO番号                       3         2         1         0
	//                               10987654321098765432109876543210
	BusBuff  |= ((Address << 16) & 0b00000000000001110000000000000000);	// 書き込んではいけないビットをANDでマスクして，アドレス部分のビットをORで書き込む
	NowState &=                    0b11111001111110001100000000000011;	// 現在のGPIOのビット列の触っていい部分をANDでクリアしてから，
	NowState |= BusBuff;			// 現在のGPIOのビット列にORで書き込んで，
	gpio.SetBitDataLo(NowState);	// GPIOからビット列出力
	
	//データバスから入力
	SetReadSelect();				// RD_SEL信号を出力
	WaitUsingNop(1000);				// dsPICのバスが揃うのが遅いので待機
	NowState = gpio.GetBitDataLo();	// 現在のGPIOのビット列を取得
	ClearReadSelect();				// RD_SEL信号を停止
	uint32_t D13_2 = NowState & 0b00000000000000000011111111111100;			// D13～D2のみ切り出す
	uint32_t D1_0  = (NowState & 0b00000110000000000000000000000000) >> 25;	// D1～D0のみ切り出して右にシフト
	return D13_2 | D1_0;			// D13～D2とD1～D0を結合して返す
}

//! @brief 絶対値多回転カウント値を取得する関数
//! @param[in]	RawCount	絶対値エンコーダの生カウント値 [-]
//! @return	絶対値多回転カウント値 [-]
int32_t SA2RasPi2::GetMultiTurnCount(uint32_t RawCount){
	//     3         2         1         0
	//    10987654321098765432109876543210
	//  0b00000000000000000000000000000000
	//    --UpperCount---^^^^^RawCount^^^^                                --UpperCount---^^^^^RawCount^^^^
	if( 0b00000000000000011000000000000000 < RawCountPrev && RawCount < 0b00000000000000000111111111111111 ){
		// オーバーフローを検出したら
		UpperCount++;	// 上位カウンタをカウントアップ
	}
	//                   --UpperCount---^^^^^RawCount^^^^      --UpperCount---^^^^^RawCount^^^^
	if( RawCountPrev < 0b00000000000000000111111111111111 && 0b00000000000000011000000000000000 < RawCount ){
		// アンダーフローを検出したら
		UpperCount--;	// 上位カウンタをカウントダウン
	}
	RawCountPrev = RawCount;	// 次回のために保持
	return (int32_t)((UpperCount << 17) | RawCount);	// 下位ビットと連結して返す
}

//! @brief インクリメンタル多回転カウント値を取得する関数
//! @param[in]	RawCount	インクリメンタルエンコーダの生カウント値 [-]
//! @return	インクリメンタル多回転カウント値 [-]
int32_t SA2RasPi2::GetIncMultiTurnCount(uint32_t RawCount){
	//     3         2         1         0
	//    10987654321098765432109876543210
	//  0b00000000000000000000000000000000
	//    -----UpperCount-----^^RawCount^^                                   -----UpperCount-----^^RawCount^^
	if( 0b00000000000000000000110000000000 < IncRawCountPrev && RawCount < 0b00000000000000000000001111111111 ){
		// オーバーフローを検出したら
		IncUpperCount++;	// 上位カウンタをカウントアップ
	}
	//                      -----UpperCount-----^^RawCount^^      -----UpperCount-----^^RawCount^^
	if( IncRawCountPrev < 0b00000000000000000000001111111111 && 0b00000000000000000000110000000000 < RawCount ){
		// アンダーフローを検出したら
		IncUpperCount--;	// 上位カウンタをカウントダウン
	}
	IncRawCountPrev = RawCount;	// 次回のために保持
	return (int32_t)((IncUpperCount << 12) | RawCount);	// 下位ビットと連結して返す
}

//! @brief 電圧[V]からDACバイナリデータへの換算関数
//! @param[in]	Vdac	DAC出力電圧 [V]
//! @return	DACバイナリデータ
uint32_t SA2RasPi2::VoltToDacData(double Vdac){
	return (uint32_t)(2047.0*(-1.0/10.0*Limiter(Vdac,10) + 1.0));
}

//! @brief デューティ比からPWMバイナリデータへの換算関数
//! @param[in]	Duty	デューティ比 [-]
//! @return	PWMバイナリデータ
uint32_t SA2RasPi2::DutyToPwmData(double Duty){
	return (uint32_t)( (double)(0x0BB7 >> 1)/0.5*Limiter(Duty,1,0) );
}

//! @brief ADCバイナリデータから電流[A]への換算関数
//! @param[in]	AdcData	ADCバイナリデータ(2の補数)
//! @return	検出電流 [A]
double SA2RasPi2::AdcDataToAmpere(uint32_t AdcData){
	// 符号拡張
	if(AdcData & (1 << 13)){
		// 符号ビットが立っていたら負数として符号拡張
		AdcData |= 0b11111111111111111100000000000000;
	}
	return ((int32_t)AdcData)/8192.0*7.5;	// ±8192(14bitの2の補数)を±1に換算して，センサの最大電流7.5Aを掛ける
}

//! @brief 絶対値エンコーダデータから位置[rad]への換算関数
//! @param[in]	AbsData	エンコーダデータ
//! @return	位置 [rad]
double SA2RasPi2::AbsDataToRadian(int32_t AbsMultData){
	// 注意事項： 多摩川の場合，エンコーダ値がゼロのとき電気角が-πである。
	return (double)AbsMultData*2.0*M_PI/131072.0 - M_PI/4.0;	// 17bitエンコーダを rad に換算
}

//! @brief モータ電気角 [rad] へ換算する関数
//! @param[in]	AbsData	エンコーダカウント値
//! @return	電気角 [rad] (0～2πの値域循環)
double SA2RasPi2::ConvElectAngle(uint32_t AbsData){
	// 注意事項： 多摩川の場合，エンコーダ値がゼロのとき電気角が-πである。
	//            なので，0x0000～0x7FFF → -π～π であるので，オフセットバイナリとみなせる。
	int32_t EleAngleData = (int32_t)(AbsData % (131072/4)) - (131072/4/2);	// 17bitエンコーダで極対数4の場合
	return (double)EleAngleData*M_PI/((double)(131072/4/2));
}


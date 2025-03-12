//! @file WEF-6A.cc
//! @brief ワコーテック社製 DynPick WEF-6A 6軸力覚センサクラス
//!
//! シリアル通信(RS422)によりWEF-6Aのセンサ情報を読み取ります。
//!
//! @date 2022/03/11
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2022 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#include <unistd.h>
#include <string>
#include <array>
#include <tuple>
#include "ARCSeventlog.hh"
#include "WEF-6A.hh"
#include "ARCSprint.hh"
#include "ArcsMatrix.hh"

using namespace ARCS;

//! @brief コンストラクタ
//! @param[in]	RS422CHx RS422シリアル通信ボード(PCI46610x)へのポインタ
WEF6A::WEF6A(std::unique_ptr<PCI46610x> RS422CHx)
	: RS422(std::move(RS422CHx)), VersionInfo(), RecordNumprev(0),
	  SFx(0), SFy(0), SFz(0), SMx(0), SMy(0), SMz(0),
	  Fxprev(0), Fyprev(0), Fzprev(0), Mxprev(0), Myprev(0), Mzprev(0),
	  IsCalibrated(false)
{
	PassedLog();
	GetVersionInfo(VersionInfo);	// バージョン情報の取得
	EventLog(VersionInfo);			// バージョン情報をイベントログへ書き込む
	GetSensitivity(SFx, SFy, SFz, SMx, SMy, SMz);	// センサ主軸感度の取得
	EventLogVar(SFx);	// センサ感度をイベントログに書き込む
	EventLogVar(SFy);	// センサ感度をイベントログに書き込む
	EventLogVar(SFz);	// センサ感度をイベントログに書き込む
	EventLogVar(SMx);	// センサ感度をイベントログに書き込む
	EventLogVar(SMy);	// センサ感度をイベントログに書き込む
	EventLogVar(SMz);	// センサ感度をイベントログに書き込む
	// オフセットキャリブレーションの実行はコンストラクタから削除しました。(2022/01/24)
}

//! @brief ムーブコンストラクタ
WEF6A::WEF6A(WEF6A&& right)
	: RS422(std::move(right.RS422)), VersionInfo(right.VersionInfo), RecordNumprev(right.RecordNumprev),
	  SFx(right.SFx), SFy(right.SFy), SFz(right.SFz), SMx(right.SMx), SMy(right.SMy), SMz(right.SMz),
	  Fxprev(right.Fxprev), Fyprev(right.Fyprev), Fzprev(right.Fzprev),
	  Mxprev(right.Mxprev), Myprev(right.Myprev), Mzprev(right.Mzprev),
	  IsCalibrated(right.IsCalibrated)
{
	right.RS422 = nullptr;	// ムーブ元はヌルポにしておく
}

//! @brief デストラクタ
WEF6A::~WEF6A(){
	PassedLog();
}

//! @brief 力覚センサのバージョン情報を取得する関数
//! @param[out]	VerInfo バージョン情報の文字列
void WEF6A::GetVersionInfo(std::string& VerInfo){
	RS422->SetTransData((uint8_t)'V');	// 力覚センサのファームウェアバージョンのリクエスト送信
	usleep(1000);						// 受信が完了したであろうと思われるまで待機(文字数が決まってなくて受信完了判断がつかないので，仕方なくスリープ)
	RS422->GetAllRecvData(VerInfo);		// 受信データ取得
}

//! @brief センサ主軸感度を取得する関数
//! @param[out]	SensFx-SensMz 確実のセンサ感度
void WEF6A::GetSensitivity(double& SensFx, double& SensFy, double& SensFz, double& SensMx, double& SensMy, double& SensMz){
	std::string SensText;				// 主軸感度の10進テキストデータ
	RS422->SetTransData((uint8_t)'p');	// 力覚センサの主軸感度のリクエスト送信
	usleep(1000);						// 受信が完了したであろうと思われるまで待機(文字数が決まってなくて受信完了判断がつかないので，仕方なくスリープ)
	RS422->GetAllRecvData(SensText);	// 受信データ取得
	sscanf(SensText.c_str(), "%lf,%lf,%lf,%lf,%lf,%lf  ", &SensFx, &SensFy, &SensFz, &SensMx, &SensMy, &SensMz);	// コンマで分解して浮動小数点に変換
}

//! @brief 6軸力覚センサ値のリクエストを送る関数
void WEF6A::SendForceRequest(void){
	RS422->ClearFIFOMemory();			// 念のため受信FIFOメモリをクリア
	RS422->SetTransData((uint8_t)'R');	// 力覚センサ値のリクエスト送信
}

//! @brief 6軸力覚センサ値の受信を待機する関数(ブロッキング動作)
void WEF6A::WaitForceData(void){
	// データが全部来るまで待機
	while(RS422->GetReceiveCount() < NUM_FORCEDATA){
		asm("nop");	// 最適化防止のためのNOP
	}
}

//! @brief センサオフセットのゼロキャリブレーションを実行する関数
void WEF6A::ZeroCalibration(void){
	// WEF-6Aの仕様に従って，３回分のオフセットキャリブレーションコマンドを送信
	RS422->SetTransData((uint8_t)'O');	// 力覚センサ値のオフセットキャリブレーションコマンド送信
	usleep(1000);						// 念のための待機
	RS422->SetTransData((uint8_t)'O');	// 力覚センサ値のオフセットキャリブレーションコマンド送信
	usleep(1000);						// 念のための待機
	RS422->SetTransData((uint8_t)'O');	// 力覚センサ値のオフセットキャリブレーションコマンド送信
	usleep(1000);						// 念のための待機
}

//! @brief センサオフセットのゼロキャリブレーションを実行する関数(1回限定版)
void WEF6A::ZeroCalibrationOneTime(void){
	if(IsCalibrated == false){
		// キャリブレーションが未実施だったら，
		RS422->SetTransData((uint8_t)'O');	// 力覚センサ値のオフセットキャリブレーションコマンド送信
		RS422->SetTransData((uint8_t)'O');	// 力覚センサ値のオフセットキャリブレーションコマンド送信
		RS422->SetTransData((uint8_t)'O');	// 力覚センサ値のオフセットキャリブレーションコマンド送信
		IsCalibrated = true;	// キャリブレーションフラグを立てておく
	}
}

//! @brief センサ内蔵の移動平均フィルタを設定する関数
//! @param[in]	FS	フィルタ設定
void WEF6A::SetInternalFilter(const FilterSettings FS){
	switch(FS){
		case FILT_DISABLE:
			RS422->SetTransData((uint8_t)'1');	// 移動平均1点(＝フィルタ無効)
			break;
		case FILT_2AVE:
			RS422->SetTransData((uint8_t)'2');	// 移動平均2点
			break;
		case FILT_4AVE:
			RS422->SetTransData((uint8_t)'4');	// 移動平均4点
			break;
		case FILT_8AVE:
			RS422->SetTransData((uint8_t)'8');	// 移動平均8点
			break;
		default:
			RS422->SetTransData((uint8_t)'1');	// 例外の場合は移動平均1点
			break;
	}
	RS422->SetTransData((uint8_t)'F');	// フィルタ設定コマンドを送出
}

//! @brief 6軸力覚センサ値を取得する関数
//! @param[out]	Fx-Fy 各軸の並進力 [N]
//! @param[out]	Mx-My 各軸のトルク [Nm]
//! @return	true = 正常終了, false = 異常検知
bool WEF6A::Get6axisForce(double& Fx, double& Fy, double& Fz, double& Mx, double& My, double& Mz){
	bool IsSuccessful = false;						// 正常終了フラグ
	uint8_t RecordNum, DATA_CR, DATA_LF;			// レコード番号, キャリッジリターン, ラインフィード
	double Fxbuf, Fybuf, Fzbuf, Mxbuf, Mybuf, Mzbuf;// [N],[Nm] 力覚センサ値のバッファ
	
	// 力覚データの構造 [N 1byte][Fx 4bytes][Fy 4bytes][Fz 4bytes][Mx 4bytes][My 4bytes][Mz 4bytes][CR][LF]
	// 注意！バイナリではなく文字データとして来る！なんでそんな仕様にしたんだ！
	// データ読み込みの順番が重要なので，下記の行を入れ替えることは不可なので注意！
	RecordNum = RS422->GetRecvData() - 0x30;// [-] レコード番号を取得してアスキーコードから数値へ変換
	Fxbuf = GetConv1axisForce(SFx);			// [N] X軸並進力の取得
	Fybuf = GetConv1axisForce(SFy);			// [N] Y軸並進力の取得
	Fzbuf = GetConv1axisForce(SFz);			// [N] Z軸並進力の取得
	Mxbuf = GetConv1axisForce(SMx);			// [Nm] X軸トルクの取得
	Mybuf = GetConv1axisForce(SMy);			// [Nm] Y軸トルクの取得
	Mzbuf = GetConv1axisForce(SMz);			// [Nm] Z軸トルクの取得
	DATA_CR = RS422->GetRecvData();			// キャリッジリターンの取得
	DATA_LF = RS422->GetRecvData();			// ラインフィードの取得
	
	// データの健全性のチェック
	if(RecordNum <=9 && DATA_CR == 0x0d && DATA_LF == 0x0a){
		// レコード番号が0～9で，最後にちゃんと CR LF が来ていれば，正常として処理するので，
		std::tie(Fx, Fy, Fz, Mx, My, Mz)
		 = std::forward_as_tuple(Fxbuf, Fybuf, Fzbuf, Mxbuf, Mybuf, Mzbuf);	// 力覚センサ値を更新
		std::tie(Fxprev, Fyprev, Fzprev, Mxprev, Myprev, Mzprev)
		 = std::forward_as_tuple(Fxbuf, Fybuf, Fzbuf, Mxbuf, Mybuf, Mzbuf);	// 次回用に力センサ値を保持しておく
		RecordNumprev = RecordNum;	// 次回用にレコード番号も記録しておく
		IsSuccessful = true;		// フラグを正常終了にセット
	}else{
		// レコード番号が変だったり，最後が CR LF じゃなかったら何かがおかしいってことで，異常として処理するので，
		std::tie(Fx, Fy, Fz, Mx, My, Mz)
		 = std::forward_as_tuple(Fxprev, Fyprev, Fzprev, Mxprev, Myprev, Mzprev);	// 力覚センサ値は更新せずに前回値を出力する
		IsSuccessful = false;		// フラグを異常検知にセット
	}
	
	return IsSuccessful;	// 正常終了か異常検知を返す
}

//! @brief 6軸力覚センサ値を取得する関数(ArcsMat版)
//! @param[out]	Force 各軸の並進力 [N], 各軸のトルク [Nm]
//! @return	true = 正常終了, false = 異常検知
bool WEF6A::Get6axisForce(ArcsMat<6,1>& Force){
	bool IsSuccessful = false;	// 正常終了フラグ
	ArcsMat<6,1> Fi={0};	// 力覚センサ値
	IsSuccessful = Get6axisForce(Fi[1], Fi[2], Fi[3], Fi[4], Fi[5], Fi[6]);	// 値を取得
	Force = Fi;					// 力覚センサ値を出力
	return IsSuccessful;		// 正常終了か異常検知を返す
}

//! @brief 1軸分のデータ取得と換算を行う関数
//! @param[in]	Sensitivity センサ主軸感度 [LSB/N],[LSB/Nm]
//! @return	1軸分の並進力かトルク [Nm]
double WEF6A::GetConv1axisForce(const double Sensitivity){
	uint8_t HHbuf, HLbuf, LHbuf, LLbuf;	// 4bytes分のバイナリデータバッファ
	uint16_t ForceData;					// 力覚データのバイナリデータ
	
	// 4バイト分を受信FIFOメモリから取得
	HHbuf = ConvAsciiToHex(RS422->GetRecvData());	// 上位の上位4bitデータ
	HLbuf = ConvAsciiToHex(RS422->GetRecvData());	// 上位の下位4bitデータ
	LHbuf = ConvAsciiToHex(RS422->GetRecvData());	// 下位の上位4bitデータ
	LLbuf = ConvAsciiToHex(RS422->GetRecvData());	// 下位の下位4bitデータ
	
	// 4bitデータ4コ分を結合して16bitのオフセットバイナリデータに変換
	ForceData = (HHbuf << 12) | (HLbuf << 8) | (LHbuf << 4) | LLbuf;
	
	// オフセットバイナリから浮動小数点に変換して，さらにSI単位系に換算して返す
	return (double)((int16_t)ForceData - 0x2000)/Sensitivity;
}

//! @brief アスキーコードから16進バイナリデータに変換する関数
//! @param[in]	Ascii	アスキーコード（0～9, 大文字A～Fのみ対応）
//! @return	16進バイナリデータ
uint8_t WEF6A::ConvAsciiToHex(const uint8_t Ascii){
	// sscanfでの変換は遅いので独自にテーブルで実装
	uint8_t Hex;
	switch(Ascii){
		case 0x30:	Hex = 0x0;	break;	// 0
		case 0x31:	Hex = 0x1;	break;	// 1
		case 0x32:	Hex = 0x2;	break;	// 2
		case 0x33:	Hex = 0x3;	break;	// 3
		case 0x34:	Hex = 0x4;	break;	// 4
		case 0x35:	Hex = 0x5;	break;	// 5
		case 0x36:	Hex = 0x6;	break;	// 6
		case 0x37:	Hex = 0x7;	break;	// 7
		case 0x38:	Hex = 0x8;	break;	// 8
		case 0x39:	Hex = 0x9;	break;	// 9
		case 0x41:	Hex = 0xA;	break;	// A
		case 0x42:	Hex = 0xB;	break;	// B
		case 0x43:	Hex = 0xC;	break;	// C
		case 0x44:	Hex = 0xD;	break;	// D
		case 0x45:	Hex = 0xE;	break;	// E
		case 0x46:	Hex = 0xF;	break;	// F
		default:	Hex = 0x0;	break;	// それ以外はとりあえずゼロ
	}
	return Hex;
}


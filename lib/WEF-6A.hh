//! @file WEF-6A.hh
//! @brief ワコーテック社製 DynPick WEF-6A 6軸力覚センサクラス
//!
//! RS422シリアル通信によりWEF-6Aのセンサ情報を読み取ります。
//!
//! @date 2022/03/11
//! @author Yokokura, Yuki
//
// Copyright (C) 2011-2022 Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the FreeBSD License.
// For details, see the License.txt file.

#ifndef WEF_6A
#define WEF_6A

#include <array>
#include <string>
#include <memory>
#include "PCI-46610x.hh"

namespace ARCS {	// ARCS名前空間
//! @brief ワコーテック社製 DynPick WEF-6A 6軸力覚センサクラス
class WEF6A {
	public:
		//! @brief センサ内蔵の移動平均フィルタ設定の定義
		enum FilterSettings {
			FILT_DISABLE,	//!< フィルタ無効
			FILT_2AVE,		//!< データ2点の移動平均
			FILT_4AVE,		//!< データ4点の移動平均
			FILT_8AVE		//!< データ8点の移動平均
		};
		
		explicit WEF6A(std::unique_ptr<PCI46610x> RS422CHx);//!< コンストラクタ
		WEF6A(WEF6A&& right);								//!< ムーブコンストラクタ
		~WEF6A();											//!< デストラクタ
		void GetVersionInfo(std::string& VerInfo);			//!< 力覚センサのバージョン情報を取得する関数
		void GetSensitivity(
			double& SensFx, double& SensFy, double& SensFz,
			double& SensMx, double& SensMy, double& SensMz
		);	//!< センサ主軸感度を取得する関数
		void SendForceRequest(void);		//!< 6軸力覚センサ値のリクエストを送る関数
		void WaitForceData(void);			//!< 6軸力覚センサ値の受信を待機する関数(ブロッキング動作)
		void ZeroCalibration(void);			//!< センサオフセットのゼロキャリブレーションを実行する関数
		void ZeroCalibrationOneTime(void);	//!< センサオフセットのゼロキャリブレーションを実行する関数(1回限定版)
		void SetInternalFilter(const FilterSettings FS);	//!< センサ内蔵の移動平均フィルタを設定する関数
		bool Get6axisForce(
			double& Fx, double& Fy, double& Fz,
			double& Mx, double& My, double& Mz
		);													//!< 6軸力覚センサ値を取得する関数
		bool Get6axisForce(std::array<double, 6>& Force);	//!< 6軸力覚センサ値を取得する関数(配列版)
		
	private:
		WEF6A(const WEF6A&) = delete;						//!< コピーコンストラクタ使用禁止
		const WEF6A& operator=(const WEF6A&) = delete;		//!< 代入演算子使用禁止
		const unsigned int NUM_FORCEDATA = 27;				//!< [bytes] 力覚センサ値のバイト数
		std::unique_ptr<PCI46610x> RS422;					//!< RS422シリアル通信ボードへのポインタ
		std::string VersionInfo;							//!< 力覚センサのファームウェアバージョン
		uint8_t RecordNumprev;	//!< レコード番号(前回値)
		double SFx;				//!< [LSB/N] X軸並進力センサ感度
		double SFy;				//!< [LSB/N] Y軸並進力センサ感度
		double SFz;				//!< [LSB/N] Z軸並進力センサ感度
		double SMx;				//!< [LSB/Nm] X軸モーメントセンサ感度
		double SMy;				//!< [LSB/Nm] Y軸モーメントセンサ感度
		double SMz;				//!< [LSB/Nm] Z軸モーメントセンサ感度
		double Fxprev;			//!< [N] X軸並進力(前回値)
		double Fyprev;			//!< [N] Y軸並進力(前回値)
		double Fzprev;			//!< [N] Z軸並進力(前回値)
		double Mxprev;			//!< [Nm] X軸トルク(前回値)
		double Myprev;			//!< [Nm] Y軸トルク(前回値)
		double Mzprev;			//!< [Nm] Z軸トルク(前回値)
		bool IsCalibrated;		//!< ゼロキャリブレーションがされたかどうかのフラグ
		double GetConv1axisForce(const double Sensitivity);	//!< 1軸分のデータ取得と換算を行う関数
		static uint8_t ConvAsciiToHex(const uint8_t Ascii);	//!< アスキーコードから16進バイナリデータに変換する関数
};
}

#endif


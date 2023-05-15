//! @file PCIe-AC01.hh
//! @brief PCIe-AC01 インタフェースクラス
//! @date 2016/12/14
//! @author Yuki YOKOKURA
//
// Copyright (C) 2011-2016 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the BSD License.
// For details, see the License.txt file.

#ifndef PCIEAC01
#define PCIEAC01

#include <unistd.h>
#include <stdint.h>
#include <string>

namespace ARCS {	// ARCS名前空間
	//! @brief PCIe-AC01 インタフェースクラス
	//!
	//! サーボアンプ(SRVAMP1)とPCI Express経由で通信するためのクラス
	//!
	class PCIeAC01 {
		public:
			PCIeAC01(const off_t BaseAddress, unsigned int CPUNum);	// コンストラクタ
			~PCIeAC01();											// デストラクタ
			static const unsigned int MaxAxis = 16;					//!< PCIe-AC01が対応できる最大の軸数
			void GetEncorderData(uint32_t EncData[MaxAxis]);		// エンコーダデータを読み込む関数
			void GetQaxisCurrent(double Iq[MaxAxis]);				// q軸電流を読み込む関数
			void SetQaxisCurrentRef(double Iqref[MaxAxis]);			// q軸電流指令を書き込む関数
			
		private:
			PCIeAC01(const PCIeAC01&);								// コピーコンストラクタ使用禁止
			const PCIeAC01& operator=(const PCIeAC01&);				// 代入演算子使用禁止
			static void ExecLinuxCommand(std::string CommandLine);	// Linuxコマンドを実行する関数
			//! @brief readデータの定義
			enum ReadDataList {
				Read_Status			= 0x0001,	//!< ステータスレジスタ
				Read_DaxisCurrent	= 0x0002,	//!< d軸電流
				Read_QaxisCurrent	= 0x0004,	//!< q軸電流
				Read_EncoderData	= 0x0008,	//!< エンコーダデータ
				Read_EchoBack		= 0x0010,	//!< エコーバック
			};
			//! @brief writeデータの定義
			enum WriteDataList {
				Write_Command			= 0x0001,	//!< コマンドレジスタ
				Write_DaxisCurrentRef	= 0x0002,	//!< d軸電流指令
				Write_QaxisCurrentRef	= 0x0004,	//!< q軸電流指令
				Write_EchoForward		= 0x0010,	//!< エコーフォワード
			};
			//! @brief uint32_t←→float 変換用共用体
			union uint32_t_float {
				uint32_t binary_data[MaxAxis];
				float float_data[MaxAxis];
			};
			const std::string ShellPath;	//!< シェルまでのファイルパス
			const std::string StartPath;	//!< 開始スクリプトまでのファイルパス
			const std::string StopPath;		//!< 停止スクリプトまでのファイルパス
			int fd;							//!< ファイルデスクリプタ
			cpu_set_t cpuset;				//!< CPU設定用変数
	};
}

#endif


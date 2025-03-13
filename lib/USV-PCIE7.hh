//! @file USV_PCIE7.cc
//! @brief USV_PCIE7インターフェースクラス
//!
//! ユニパルス社製(ロボテック社製)ユニサーボの通信規格UNINET
//! 通信用USV-PCIE7ボードとのインターフェースを提供します
//!
//! @date 2025/03/13
//! @author Juan Padron and Yokokura, Yuki
//
// Copyright (C) 2011-2025 Juan Padron and Yokokura, Yuki
// This program is free software;
// you can redistribute it and/or modify it under the terms of the BSD License.
// For details, see the License.txt file.

#ifndef USVPCIE7
#define USVPCIE7

#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <cstdint>
#include <cmath>
#include <array>
#include <cassert>

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

namespace ARCS {
//! @brief USV_PCIE7インターフェースクラス
//! ユニパルスUSV_PCIE7コントローラのための入出力機能を提供します。
//! @tparam N	使用する軸数
template <size_t N>
class USV_PCIE7{
	public:
		//! @brief 制御モードの定義
		enum ControlMode {
			ACCELERATION_CTRL,	//!< 加速度制御モード
			CURRENT_CTRL		//!< 電流制御モード
		};

		//! @brief 空コンストラクタ
		USV_PCIE7(void)
			: ADDR_BASE(0), fd(0), Memptr(nullptr)
		{
			// 空コンストラクタ
		}
		
		//! @brief コンストラクタ
		//! @param[in] Addr ベースアドレス(lspciコマンドで表示される「Region1」のアドレス)
		USV_PCIE7(const size_t Addr)
			: ADDR_BASE(Addr), fd(0), Memptr(nullptr)
		{
			PassedLog();
			fd = open("/dev/mem", O_RDWR|O_SYNC);	// メモリ空間へのfile descriptor取得
			arcs_assert(fd != -1);					// エラーがないか確認
			Memptr = (uint32_t*)mmap(nullptr, MEMMAP_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, fd, ADDR_BASE); //メモリ空間へのメモリーマップポインターを取得
			close(fd);
			arcs_assert(Memptr != MAP_FAILED);		//ちゃんとポインタ取得できたか確認
		}

		//! @brief ムーブコンストラクタ
		//! @param[in] right 右辺
		USV_PCIE7(USV_PCIE7&& right)
			: ADDR_BASE(right.ADDR_BASE), fd(0), Memptr(right.Memptr)
		{
			right.Memptr = nullptr;		// ムーブ元はヌルポにしておく
		}

		//! @brief デストラクタ
		~USV_PCIE7(){
			PassedLog();
			SetZeroCurrent();			// 念のためのゼロ電流指令
			munmap(Memptr, MEMMAP_SIZE);// メモリマッピング解除
		}

		//! @brief 指定した軸にq軸電流指令を設定する関数
		//! @param[in]	Current q軸電流指令 [A]
		//! @param[in]	Axis 指定の軸番号(1～7)
		void SetCurrent(const double Current, const size_t Axis){
			arcs_assert(1 <= Axis && Axis <= 7);		// 軸の範囲チェック
			UintAndFloat CurrentUnion;					// float→uint32_t変換用共用体
			CurrentUnion.FloatValue = -(float)Current;	// 倍精度から単精度へキャストしてfloatとして共用体に格納，符号をエンコーダに合わせる
			Memptr[Axis - 1] = CurrentUnion.BinaryExpr;	// uint32_tとして共用体から読み出して，q軸電流指令に対応したメモリに書き込み
		}
		
		//! @brief 指定した軸にq軸電流指令を設定する関数(配列版)
		//! @param[in]	Current q軸電流指令配列 [A]
		void SetCurrent(const std::array<double, N>& Current){
			for(size_t i = 1; i <= N; ++i) SetCurrent(Current[i - 1], i);
		}
		
		//! @brief 全軸のq軸電流指令をゼロに設定する関数
		void SetZeroCurrent(void){
			for(size_t i = 1; i <= N; ++i) SetCurrent(0, i);
		}
		
		//! @brief 指定した軸のq軸電流を取得する関数(アナログ出力設定"AOMD 4 0 1.0"に設定のこと)
		//! @param[in]	Axis 指定の軸番号(1～7)
		//! @return	q軸電流 [A]
		double GetCurrent(const size_t Axis){
			arcs_assert(1 <= Axis && Axis <= 7);		// 軸の範囲チェック
			UintAndFloat CurrentUnion;					// uint32_t→float変換用共用体
			CurrentUnion.BinaryExpr = Memptr[IDX_RDBASE + IDX_PER_AXIS*(Axis - 1) + IDX_GNRL];	// q軸電流に対応したメモリから読み込んで，uint32_tとして共用体に格納
			return (double)CurrentUnion.FloatValue;		// [A] 単精度から倍精度へキャストして返す
		}
		
		//! @brief 全軸のq軸電流指令を取得する関数(配列版) (アナログ出力設定"AOMD 4 0 1.0"に設定のこと)
		//! @param[out]	q軸電流配列 [A]
		void GetCurrent(std::array<double, N>& Current){
			for(size_t i = 1; i <= N; ++i) Current[i - 1] = GetCurrent(i);
		}
		
		//! @brief 指定した軸の速度を取得する関数(アナログ出力設定"AOMD 1 0 1.0"に設定のこと)
		//! @param[in]	Axis 指定の軸番号(1～7)
		//! @return	速度 [rad/s]
		double GetVelocity(const size_t Axis){
			arcs_assert(1 <= Axis && Axis <= 7);		// 軸の範囲チェック
			UintAndFloat VelocityUnion;					// uint32_t→float変換用共用体
			VelocityUnion.BinaryExpr = Memptr[IDX_RDBASE + IDX_PER_AXIS*(Axis - 1) + IDX_GNRL];	// 速度に対応したメモリから読み込んで，uint32_tとして共用体に格納
			return ((double)VelocityUnion.FloatValue)*RPM_TO_RAD;	// [rad/s] 単精度から倍精度へキャストしてRPM→rad/s換算して返す
		}
		
		//! @brief 全軸の速度を取得する関数(配列版) (アナログ出力設定"AOMD 1 0 1.0"に設定のこと)
		//! @param[out]	速度配列 [rad/s]
		void GetVelocity(std::array<double, N>& Velocity){
			for(size_t i = 1; i <= N; ++i) Velocity[i - 1] = GetVelocity(i);
		}
		
		//! @brief 指定した軸の位置を取得する関数(値域0～2π版)
		//! @param[in]	Axis 指定の軸番号(1～7)
		//! @return	位置 [rad]
		double GetPosition(const size_t Axis){
			arcs_assert(1 <= Axis && Axis <= 7);		// 軸の範囲チェック
			return ((double)(ENC_MASK & Memptr[IDX_RDBASE + IDX_PER_AXIS*(Axis - 1) + IDX_POSL]))*PULSE_TO_RAD;	// 位置に対応したメモリから読み込んで，マスクかけてpulse→rad換算して返す
		}
		
		//! @brief 全軸の位置を取得する関数(値域0～2π版)(配列版)
		//! @param[out]	位置配列 [rad]
		void GetPosition(std::array<double, N>& Position){
			for(size_t i = 1; i <= N; ++i) Position[i - 1] = GetPosition(i);
		}
		
		//! @brief 指定した軸の回転の数を取得する関数
		//! @param[in]	Axis 指定の軸番号(1～7)
		//! @return	回転の数 [-]
		int32_t GetNumOfRotation(const size_t Axis){
			arcs_assert(1 <= Axis && Axis <= 7);		// 軸の範囲チェック
			UintAndInt NumOfRotUnion;					// uint32_t→int32_t変換用共用体
			NumOfRotUnion.UnsignedExpr = Memptr[IDX_RDBASE + IDX_PER_AXIS*(Axis - 1) + IDX_POSH];	// 回転の数に対応したメモリから読み込む
			return NumOfRotUnion.SignedExpr;			// 符号付きで返す
		}
		
		//! @brief 指定した軸の位置を取得する関数(値域±∞版)
		//! @param[in]	Axis 指定の軸番号(1～7)
		//! @return	位置 [rad]
		double GetFullPosition(const size_t Axis){
			arcs_assert(1 <= Axis && Axis <= 7);		// 軸の範囲チェック
			return (double)GetNumOfRotation(Axis)*2.0*M_PI + GetPosition(Axis);	// 「回転の数の分のrad」と「0～2π」を足して返す
		}
		
		//! @brief 全軸の位置を取得する関数(値域±∞版)(配列版)
		//! @param[out]	位置配列 [rad]
		void GetFullPosition(std::array<double, N>& Position){
			for(size_t i = 1; i <= N; ++i) Position[i - 1] = GetFullPosition(i);
		}
		
		//! @brief 指定した軸のトルクを取得する関数
		//! @param[in]	Axis 指定の軸番号(1～7)
		//! @return	トルク [Nm]
		double GetTorque(const size_t Axis){
			arcs_assert(1 <= Axis && Axis <= 7);		// 軸の範囲チェック
			UintAndFloat TorqueUnion;					// uint32_t→float変換用共用体
			TorqueUnion.BinaryExpr = Memptr[IDX_RDBASE + IDX_PER_AXIS*(Axis - 1) + IDX_TRQ];	// トルクに対応したメモリから読み込んで，uint32_tとして共用体に格納
			return ((double)TorqueUnion.FloatValue);	// [Nm] 単精度から倍精度へキャストして返す
		}
		
		//! @brief 全軸のトルクを取得する関数(配列版)
		//! @param[out]	トルク配列 [rad]
		void GetTorque(std::array<double, N>& Torque){
			for(size_t i = 1; i <= N; ++i) Torque[i - 1] = GetTorque(i);
		}
		
		//! @brief 指定した軸の同期カウンタデータを取得する関数
		//! @param[in]	Axis 指定の軸番号(1～7)
		//! @return	同期カウンタデータ
		uint32_t GetCounter(const size_t Axis){
			arcs_assert(1 <= Axis && Axis <= 7);		// 軸の範囲チェック
			return Memptr[IDX_CNTBASE + Axis - 1];		// ユニネットデータ同期用カウント値を読み出して返す
		}

		//! @brief 全軸サーボONにする関数
		void TurnServosOn(void){
			// 書き込み値は無視され，書き込むとコマンドが全軸に同時に単純に発行される。
			Memptr[IDX_SERVO_ON] = 0x01;
			sleep(2);	// サーボON状態が落ち着くまで2秒待機
		}

		//! @brief 全軸サーボOFFにする関数
		void TurnServosOff(void){
			// 書き込み値は無視され，書き込むとコマンドが全軸に同時に単純に発行される。
			Memptr[IDX_SERVO_OFF] = 0x01;
		}

		//! @brief 制御モードを設定する関数
		void SetControlMode(ControlMode Ctrl){
			// 書き込み値は無視され，書き込むとコマンドが全軸に同時に単純に発行される。
			switch(Ctrl){
				case ControlMode::ACCELERATION_CTRL:
					Memptr[IDX_ACCCTRL] = 0x01;	// 加速度制御モード
					break;
				case ControlMode::CURRENT_CTRL:
					Memptr[IDX_CURCTRL] = 0x01;	// 電流制御モード
					break;
				default:
					break;
			}
		}

		//! @brief エラー状態をクリアする関数
		void ClearErrorStatus(void){
			// 書き込み値は無視され，書き込むとコマンドが全軸に同時に単純に発行される。
			Memptr[IDX_ERRCLR] = 0x01;
		}

		//! @brief トルクセンサのゼロ校正を実行する関数
		void CalibrateZeroTorque(void){
			// 書き込み値は無視され，書き込むとコマンドが全軸に同時に単純に発行される。
			Memptr[IDX_ZEROCAL] = 0x01;
		}

	private:
		USV_PCIE7(const USV_PCIE7&) = delete;					//!< コピーコンストラクタ使用禁止
		const USV_PCIE7& operator=(const USV_PCIE7&) = delete;	//!< 代入演算子使用禁止
		static constexpr unsigned int MEMMAP_SIZE = 0x120;		//!< メモリマップサイズ
		static constexpr size_t IDX_RDBASE = 32;				//!< 読み込み用ベース要素番号 (0x80/4bytes = 32)
		static constexpr size_t IDX_PER_AXIS = 4;				//!< 1軸当たりの要素数 (1word = 4bytes)
		static constexpr size_t IDX_POSH = 0;					//!< モータ回転数の相対要素番号 (例：1軸は0x80)
		static constexpr size_t IDX_POSL = 1;					//!< エンコーダ値の相対要素番号 (例：1軸は0x84)
		static constexpr size_t IDX_TRQ  = 2;					//!< 計測トルクの相対要素番号   (例：1軸は0x88)
		static constexpr size_t IDX_GNRL = 3;					//!< 汎用出力の相対要素番号     (例：1軸は0x8C)
		static constexpr size_t IDX_CNTBASE   = 60;				//!< 同期カウンタ用ベース要素番号 (0xF0/4bytes = 60)
		static constexpr size_t IDX_SERVO_ON  = 20;			   	//!< サーボオンコマンドの要素番号 (0x50/4bytes = 20) "SV 1"
		static constexpr size_t IDX_SERVO_OFF = 21;	 			//!< サーボオフコマンドの要素番号 (0x54/4bytes = 21) "SV 0"
		static constexpr size_t IDX_ACCCTRL   = 22;				//!< 加速度制御モード要素番号 (0x58/4bytes = 22) "CTM 11"
		static constexpr size_t IDX_CURCTRL   = 23;				//!< 電流制御モード要素番号   (0x5C/4bytes = 23) "CTM 13"
		static constexpr size_t IDX_ERRCLR    = 24;				//!< エラークリア要素番号     (0x60/4bytes = 24) "SVERCLR"
		static constexpr size_t IDX_ZEROCAL   = 25;				//!< トルクゼロ校正要素番号   (0x64/4bytes = 25) "TRQORG"
		static constexpr double RPM_TO_RAD   = 2.0*M_PI/60.0;		//!< RPMからrad/sへの換算係数
		static constexpr double PULSE_TO_RAD = 2.0*M_PI/1048576.0;	//!< RPMからrad/sへの換算係数
		static constexpr uint32_t ENC_MASK = 0x000FFFFF;		//!< エンコーダデータ不要ビットマスク
		const unsigned long ADDR_BASE;							//!< ベースアドレス
		int fd;  												//!< ファイル記述子
		uint32_t* Memptr;										//!< PCIEメモリへのポインタ
		
		//! @brief uint32_t←→float 変換用共用体
		union UintAndFloat {
			uint32_t BinaryExpr;	//!< バイナリ表現
			float FloatValue;		//!< 単精度浮動小数点での数値
		};
		
		//! @brief uint32_t←→int32_t 変換用共用体
		union UintAndInt {
			uint32_t UnsignedExpr;	//!< 符号なしバイナリ表現
			int32_t  SignedExpr;	//!< 符号ありバイナリ表現
		};
		
};
}

#endif

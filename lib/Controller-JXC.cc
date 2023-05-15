//! @file Controller-JXC.cc
//! @brief JXCシリーズコントローラでの電動グリッパの制御
//! 動作環境 PCI2826CV >> JXC5171-BC-E >> LEHF32EK2-32
//! CentOS6.10, ARCS6オフラインモードにて確認
//! @date 2022/03/18
//! @author Nagisa Kuwata
//
// Copyright (C) 2011-2020 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the BSD License.
// For details, see the License.txt file.

#include <stdint.h>
#include "Matrix.hh"
#include "Controller-JXC.hh"

using namespace ARCS;

//! @brief コンストラクタ
//! @param [in] Addr PCIボードアドレス
//! @param [in] Ts 制御周期
//! @param [in] InputPortArray コントローラ側入力ポートリスト
//! @param [in] InputBitArray コントローラ側入力ビットリスト
//! @param [in] OutputPortArray コントローラ側出力ポートリスト
//! @param [in] OutputBitArray コントローラ側出力ビットリスト
ControllerJXC::ControllerJXC(
    unsigned int Addr,
    double Ts,
    std::array<unsigned short int, 11> InputPortArray,
    std::array<unsigned short int, 11> InputBitArray,
    std::array<unsigned short int, 13> OutputPortArray,
    std::array<unsigned short int, 13> OutputBitArray)
  : ControllerJxc(Addr),
    Ts(Ts),
    InputPortArray(InputPortArray),
    InputBitArray(InputBitArray),
    OutputPortArray(OutputPortArray),
    OutputBitArray(OutputBitArray)
{
    PortIn0 = InputPortArray[0];
    PortIn1 = InputPortArray[1];
    PortIn2 = InputPortArray[2];
    PortIn3 = InputPortArray[3];
    PortIn4 = InputPortArray[4];
    PortIn5 = InputPortArray[5];
    PortSetup = InputPortArray[6];
    PortHold = InputPortArray[7];
    PortDrive = InputPortArray[8];
    PortReset = InputPortArray[9];
    PortSvon = InputPortArray[10];
    BitIn0 = InputBitArray[0];
    BitIn1 = InputBitArray[1];
    BitIn2 = InputBitArray[2];
    BitIn3 = InputBitArray[3];
    BitIn4 = InputBitArray[4];
    BitIn5 = InputBitArray[5];
    BitSetup = InputBitArray[6];
    BitHold = InputBitArray[7];
    BitDrive = InputBitArray[8];
    BitReset = InputBitArray[9];
    BitSvon = InputBitArray[10];
    PortOut0 = OutputPortArray[0];
    PortOut1 = OutputPortArray[1];
    PortOut2 = OutputPortArray[2];
    PortOut3 = OutputPortArray[3];
    PortOut4 = OutputPortArray[4];
    PortOut5 = OutputPortArray[5];
    PortBusy = OutputPortArray[6];
    PortArea = OutputPortArray[7];
    PortSeton = OutputPortArray[8];
    PortInp = OutputPortArray[9];
    PortSvre = OutputPortArray[10];
    PortEstop = OutputPortArray[11];
    PortAlarm = OutputPortArray[12];
    BitOut0 = OutputBitArray[0];
    BitOut1 = OutputBitArray[1];
    BitOut2 = OutputBitArray[2];
    BitOut3 = OutputBitArray[3];
    BitOut4 = OutputBitArray[4];
    BitOut5 = OutputBitArray[5];
    BitBusy = OutputBitArray[6];
    BitArea = OutputBitArray[7];
    BitSeton = OutputBitArray[8];
    BitInp = OutputBitArray[9];
    BitSvre = OutputBitArray[10];
    BitEstop = OutputBitArray[11];
    BitAlarm = OutputBitArray[12];
}

//! @brief 空コンストラクタ
ControllerJXC::ControllerJXC()
  : ControllerJxc(0)
{
}

//! @brief デストラクタ
ControllerJXC::~ControllerJXC(){
}

//---------------//
// サーボオフ関数
//---------------//
//! @brief サーボオン信号をオフする.
void ControllerJXC::SvonOff(void) {
    if (StateSvon) {
        ControllerJxc.SetData(0, PortSvon, BitSvon);
        StateSvon = false;
    }
}

//-------------//
// リセット関数
//-------------//
//! @brief 各出力, 及びエラーをリセットする.
//! 動作中は呼び出し続け, 戻り値を監視して運用する.
//! @param [in] Ts 制御周期[s]
//! @return リセット動作中:true, 動作終了状態:false
bool ControllerJXC::AllReset(void) {
    static double time_count = 0;
    if (!StateReset) {
        ControllerJxc.SetData(0, PortIn0, BitIn0);
        ControllerJxc.SetData(0, PortIn1, BitIn1);
        ControllerJxc.SetData(0, PortIn2, BitIn2);
        ControllerJxc.SetData(0, PortIn3, BitIn3);
        ControllerJxc.SetData(0, PortIn4, BitIn4);
        ControllerJxc.SetData(0, PortIn5, BitIn5);
        ControllerJxc.SetData(0, PortSetup, BitSetup);
        ControllerJxc.SetData(0, PortHold, BitHold);
        ControllerJxc.SetData(0, PortDrive, BitDrive);
        ControllerJxc.SetData(0, PortReset, BitReset);
        ControllerJxc.SetData(1, PortReset, BitReset);
        printf("reset:on\n");
        StateReset = true;
        ControllerJxc.SetData(1, PortSvon, BitSvon);
        printf("svon:on\n");
        StateSvon = true;
    } else {
        time_count += Ts;
        if (time_count > 0.5) {
            ControllerJxc.SetData(0, PortReset, BitReset);
            printf("reset:off\n");
            StateReset = false;
            time_count = 0;
        }
    }
    return StateReset;
}

//-------------//
// 原点復帰関数
//-------------//
//! @brief 原点復帰を行う. 
//! 動作中は呼び出し続け, 戻り値を監視して運用する.
//! @return 原点復帰動作中:true, 動作終了状態:false
bool ControllerJXC::ReturnOrg(void) {
    static double time_count = 0;
    if (!StateOrg) {
        ControllerJxc.SetData(1, PortSetup, BitSetup);
        printf("org:on\n");
        StateOrg = true;
    } else if (!ControllerJxc.GetData(PortBusy, BitBusy)) {
        time_count += Ts;
        if (time_count > 0.5) {
            ControllerJxc.SetData(0, PortSetup, BitSetup);
            printf("org:off\n");
            StateOrg = false;
            time_count = 0;
        }
    }
    return StateOrg;
}

//-----------------//
// パターン出力関数
//-----------------//
//! @brief 事前に設定された0~63の64パターンの動作を出力する.
//! 動作中は呼び出し続け, 戻り値を監視して運用する.
//! 動作中のパターン変更可
//! @param [in] PatternID 動作パターンのID(0~63)
//! @return 動作中:true, 動作終了状態:false
bool ControllerJXC::PatternOutput(const short int PatternID) {
    static double time_count = 0;
    static bool time_delay = false;
    if (StateOrg) {
        ControllerJxc.SetData(0, PortSetup, BitSetup);
        printf("org:off\n");
        StateOrg = false;
    }
    if (PrePatternID == -1 || PrePatternID != PatternID) {
        for (int i = 0; i < 6; i++) {
            bool output;
            PrePatternID = PatternID;
            if (PatternID & int(std::pow(2, i))) output = true;
            else output = false;
            switch (i) {
                case 0:
                    ControllerJxc.SetData(output, PortIn0, BitIn0);
                    break;
                case 1:
                    ControllerJxc.SetData(output, PortIn1, BitIn1);
                    break;
                case 2:
                    ControllerJxc.SetData(output, PortIn2, BitIn2);
                    break;
                case 3:
                    ControllerJxc.SetData(output, PortIn3, BitIn3);
                    break;
                case 4:
                    ControllerJxc.SetData(output, PortIn4, BitIn4);
                    break;
                case 5:
                    ControllerJxc.SetData(output, PortIn5, BitIn5);
                    break;
            }
        }
        printf("pattern%d:on\n", PatternID);
        ControllerJxc.SetData(1, PortDrive, BitDrive);
        printf("drive:on\n");
        if (!StatePattern) StatePattern = true;
    } else if (!ControllerJxc.GetData(PortBusy, BitBusy)) {
        time_count += Ts;
        if (time_count > 0.5 && !time_delay) {
            ControllerJxc.SetData(0, PortDrive, BitDrive);
            printf("drive:off\n");
            time_count = 0;
            time_delay = true;
        } else if (time_count > 0.5 && time_delay) {
            time_count = 0;      
            PrePatternID = -1;
            time_delay = false;
            StatePattern = false;
        }
    }
    return StatePattern;
}
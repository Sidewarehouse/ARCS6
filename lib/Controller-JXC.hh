//! @file Controller-JXC.hh
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

#ifndef CONTROLLER_JXC
#define CONTROLLER_JXC

#include "Matrix.hh"
#include "PCI-2826CV.hh"

namespace ARCS {
    //! @brief JXCシリーズコントローラでの電動グリッパの制御クラス
    class ControllerJXC {
        public:
            //!コンストラクタ
            ControllerJXC(
                unsigned int Addr,
                double Ts,
                std::array<unsigned short int, 11> InputPortArray,
                std::array<unsigned short int, 11> InputBitArray,
                std::array<unsigned short int, 13> OutputPortArray,
                std::array<unsigned short int, 13> OutputBitArray);
            //!空コンストラクタ
            ControllerJXC();
            //!デストラクタ
            ~ControllerJXC();
            //!リセット関数
            bool AllReset(void);
            //!原点復帰関数
            bool ReturnOrg(void);
            //!パターン出力関数
            bool PatternOutput(const short int PatternID);
        private:
            PCI2826CV ControllerJxc;
            //!コピーコンストラクタ使用禁止
            ControllerJXC(const ControllerJXC&) = delete;
            //!代入演算子使用禁止
            const ControllerJXC& operator = (const ControllerJXC&) = delete;
            //!ムーブコンストラクタ使用禁止
            ControllerJXC(ControllerJXC&& r) = delete;
            //!制御周期
            double Ts = 0;
            //!コントローラ側入力ポートリスト
            std::array<unsigned short int, 11> InputPortArray{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            //!コントローラ側入力ビットリスト
            std::array<unsigned short int, 11> InputBitArray{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            //!コントローラ側出力ポートリスト
            std::array<unsigned short int, 13> OutputPortArray{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            //!コントローラ側出力ビットリスト
            std::array<unsigned short int, 13> OutputBitArray{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
            //!状態保持変数
            bool StateDrive = false;
            bool StateOrg = false;
            bool StatePattern = false;
            bool StateReset = false;
            bool StateSvon = false;
            //!パターン出力前回値
            short int PrePatternID = -1;
            //!コントローラ側入力ポート
            unsigned short int PortIn0 = 0;
            unsigned short int PortIn1 = 0;
            unsigned short int PortIn2 = 0;
            unsigned short int PortIn3 = 0;
            unsigned short int PortIn4 = 0;
            unsigned short int PortIn5 = 0;
            unsigned short int PortSetup = 0;
            unsigned short int PortHold = 0;
            unsigned short int PortDrive = 0;
            unsigned short int PortReset = 0;
            unsigned short int PortSvon = 0;
            //!コントローラ側入力ビット
            unsigned short int BitIn0 = 0;
            unsigned short int BitIn1 = 0;
            unsigned short int BitIn2 = 0;
            unsigned short int BitIn3 = 0;
            unsigned short int BitIn4 = 0;
            unsigned short int BitIn5 = 0;
            unsigned short int BitSetup = 0;
            unsigned short int BitHold = 0;
            unsigned short int BitDrive = 0;
            unsigned short int BitReset = 0;
            unsigned short int BitSvon = 0;
            //!コントローラ側出力ポート
            unsigned short int PortOut0 = 0;
            unsigned short int PortOut1 = 0;
            unsigned short int PortOut2 = 0;
            unsigned short int PortOut3 = 0;
            unsigned short int PortOut4 = 0;
            unsigned short int PortOut5 = 0;
            unsigned short int PortBusy = 0;
            unsigned short int PortArea = 0;
            unsigned short int PortSeton = 0;
            unsigned short int PortInp = 0;
            unsigned short int PortSvre = 0;
            unsigned short int PortEstop = 0;
            unsigned short int PortAlarm = 0;
            //!コントローラ側出力ビット
            unsigned short int BitOut0 = 0;
            unsigned short int BitOut1 = 0;
            unsigned short int BitOut2 = 0;
            unsigned short int BitOut3 = 0;
            unsigned short int BitOut4 = 0;
            unsigned short int BitOut5 = 0;
            unsigned short int BitBusy = 0;
            unsigned short int BitArea = 0;
            unsigned short int BitSeton = 0;
            unsigned short int BitInp = 0;
            unsigned short int BitSvre = 0;
            unsigned short int BitEstop = 0;
            unsigned short int BitAlarm = 0;
            //!サーボオフ関数
            void SvonOff(void);
    };
}

# endif

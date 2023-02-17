//! @file PCI-2826.hh
//! @brief PCI-2826入出力クラス
//! Interface社製PCI-2826のための入出力機能を提供します。
//! @date 2020/06/09
//! @author Yuto KOBAYASHI
//
// Copyright (C) 2011-2020 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the BSD License.
// For details, see the License.txt file.

#ifndef PCI_2826CV
#define PCI_2826CV

namespace ARCS {
    //! @brief PCI-2826入出力クラス
    //! Interface社製PCI-2826のための入出力機能を提供します。
    class PCI2826CV {
        public:
            PCI2826CV(unsigned int Addr, uint32_t GetDataMask, uint32_t SetDataMask);
            PCI2826CV(unsigned int Addr); //!< コンストラクタ
            PCI2826CV();                  //!< 空コンストラクタ
            ~PCI2826CV();                 //!< デストラクタ
            

            uint32_t GetData(void);    //!< 全入力の読み取り
            uint8_t GetData(uint8_t PortNum);  //!< あるポート入力の読み取り
            bool GetData(uint8_t PortNum, uint8_t BitNum);   //!< あるポート・ビット入力の読み取り 

            void SetData(uint32_t Data);   //!< 全出力の書き込み
            void SetData(uint8_t Data, uint8_t PortNum);    //!< あるポート出力の書き込み
            void SetData(uint8_t Data, uint8_t PortNum, uint8_t BitNum);    //!< あるポート・ビットの書き込み

            uint32_t InputData;      //!< 入力データを保持（privateにすべき？)
            uint32_t OutputData;     //!< 出力データを保持（privateにすべき？)

        private:
            PCI2826CV(const PCI2826CV&) = delete;       //!< コピーコンストラクタ使用禁止
            const PCI2826CV& operator=(const PCI2826CV&) = delete;  //!< 代入演算子使用禁止

            const unsigned int ADDR_BASE;   //!< ベースアドレス
            const unsigned int ADDR_DIN0;   //!< デジタル入力 ポート0
            const unsigned int ADDR_DIN1;   //!< デジタル入力 ポート1
            const unsigned int ADDR_DIN2;   //!< デジタル入力 ポート2
            const unsigned int ADDR_DIN3;   //!< デジタル入力 ポート3
            const unsigned int ADDR_DOUT0;   //!< デジタル出力 ポート0
            const unsigned int ADDR_DOUT1;   //!< デジタル出力 ポート1
            const unsigned int ADDR_DOUT2;   //!< デジタル出力 ポート2
            const unsigned int ADDR_DOUT3;   //!< デジタル出力 ポート3
            
            const uint32_t GET_DATA_MASK;         //!< 入力用データマスク
            const uint32_t SET_DATA_MASK;         //!< 出力用データマスク

            //データ中任意のビットを更新する
            uint32_t DataUpdate32(uint32_t NewData, uint32_t OrgData, uint32_t BitArray);
            uint8_t DataUpdate8(uint8_t NewData, uint8_t OrgData, uint8_t BitArray);
    };
}

# endif

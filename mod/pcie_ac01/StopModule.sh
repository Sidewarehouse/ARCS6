#!/bin/bash
# PCIe-AC01デバイスドライバ停止用スクリプト
# 2016/11/25  Yuki YOKOKURA
#echo -n "Kernel Module Removing..."
# モジュールの停止
rmmod pcie_ac01.ko
#echo -n "..."
# デバイスアクセス用ファイルの削除
rm -f /dev/pcie_ac01
#echo "...Done"

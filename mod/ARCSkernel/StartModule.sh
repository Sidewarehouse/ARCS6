#!/bin/bash
# PCIe-AC01デバイスドライバ起動用スクリプト
# 2017/01/30  Yuki YOKOKURA
#echo -n "Kernel Module Starting..."
current_path=`dirname $0`
# 一旦削除
rmmod pcie_ac01.ko >/dev/null 2>&1
rm -f /dev/pcie_ac01 >/dev/null 2>&1
# モジュールの起動
insmod ${current_path}/pcie_ac01.ko MemAddress=$1
#echo -n "..."
# デバイスのメジャー番号の取得
MajorNum=$(awk "\$2==\"pcie_ac01\" {print \$1}" /proc/devices)
#echo -n "..."
# デバイスアクセス用ファイルの作成
mknod /dev/pcie_ac01 c $MajorNum 0
#echo "...Done"

#!/bin/bash
# PCIeデバイス駆動用FPGAのリコンフィグスクリプト
# 2016/11/16  Yuki YOKOKURA
echo -n "Please reconfigure the FPGA device for PCIe interfaces. [Enter]"
read KeyWait
sh -c "echo 1 > /sys/bus/pci/devices/0000:04:00.0/remove"
sleep 1
sh -c "echo 1 > /sys/bus/pci/rescan"
sh -c "echo 1 > /sys/bus/pci/devices/0000:04:00.0/enable"
echo "PCIe restart complete."

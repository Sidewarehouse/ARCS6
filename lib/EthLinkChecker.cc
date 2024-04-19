// イーサネットリンク確認関数
// 2014/02/27 Yuki YOKOKURA
//
// 指定したNICがリンクアップしてるかどうか確認して状態を返す関数
//
// Copyright (C) 2014 Yuki YOKOKURA
// This program is free software;
// you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3 of the License, or any later version.
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details <http://www.gnu.org/licenses/>.
// Besides, you can negotiate about other options of licenses instead of GPL.
// If you would like to get other licenses, please contact us<yuki@katsura.sd.keio.ac.jp>.

#include "EthLinkChecker.hh"

bool EthLinkChecker(const char* eth){
	// NICの名前 true: リンク確立中, false: 切断中 or 何かのエラー
	
	struct ethtool_value ethval;
	struct ifreq          ifr;
	int err;
	int fd;
	
	memset(&ifr, 0, sizeof(struct ifreq));
	strncpy(ifr.ifr_name, eth, IFNAMSIZ);
	
	ethval.cmd = ETHTOOL_GLINK;
	ifr.ifr_data = (caddr_t)&ethval;
	
	fd = socket(AF_INET, SOCK_DGRAM, 0);
	if( 0 > fd ){
		return false;
	}
	
	err = ioctl(fd, SIOCETHTOOL, &ifr);
	close(fd);
	
	if(err==0){
		if(ethval.data){
			// リンク確立中の場合は true を返す
			return true;
		}else{
			// 切断中の場合は false を返す
			return false;
		}
	}else{
		// 何らかのエラーの場合も false を返す
		return false;
	}
	
	return false;
}



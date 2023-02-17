# ARCS6 サブディレクトリアドオン用Makefile
# 2021/07/31

# 最終的に出来上がる中間ファイルの名前
OBJNAME = ARCS_ADDON.o

# 追加ライブラリ中間ファイルへのパス
LIBNAME = ./c/ARCS_ADDONC.o

# ディレクトリパス
LIBPATH = ../../../lib
SRCPATH = ../../../src
SYSPATH = ../../../sys

# ライブラリヘッダファイルへのパス
LIBHEAD = -I. -I$(LIBPATH) -I$(SRCPATH) -I$(SYSPATH) -I/usr/src/linux/include

# コンパイラとリンカの指定
CC = gcc
LD = ld

# コンパイラに渡すオプション
CFLAGS  = -Wall -Weffc++ -std=c++17 -ftree-vectorize -march=native -fPIE -DARCS_IN

# リンカに渡すオプション
LDFLAGS = -L.

# C++ファイル形式の指定
.SUFFIXES:.cc .hh .o
SRC=$(shell ls *.cc)
OBJ=$(SRC:.cc=.o)

# ディレクトリ内すべてコンパイル＆リンク
all: $(OBJNAME)

# ディレクトリ内すべてコンパイル＆リンク(オフライン計算用)
offline: $(OBJNAME)

# 依存関係ファイル読み込み
-include Makefile.depend

# ディレクトリ内すべてリンク
$(OBJNAME): $(OBJ)
	@cd ./c; $(MAKE) --no-print-directory
	@$(LD) $(LDFLAGS) -r -o $(OBJNAME) $(OBJ) $(LIBNAME)

# ディレクトリ内すべてコンパイル
.cc.o:
	$(CC) -O2 $(LIBHEAD) $(CFLAGS) -c $<

# 掃除
.PHONY: clean
clean:
	@rm -f *.o
	@rm -f *.s
	@$(MAKE) --no-print-directory clean -C ./c

# 依存関係ファイルの生成
.PHONY: dep
dep:
	@$(CC) $(LIBHEAD) $(CFLAGS) -MM $(SRC) > Makefile.depend
	@cd ./c; $(MAKE) --no-print-directory dep

# デバッグ
.PHONY: debug
debug:
	@rm -f *.o
	@rm -f *.s
	@cd ./c; $(MAKE) --no-print-directory debug
	$(CC) -ggdb3 -Og $(LIBHEAD) $(CFLAGS) -c $(SRC)
	$(LD) $(LDFLAGS) -r -o $(OBJNAME) $(OBJ) $(LIBNAME)

# すべてコンパイル (オフライン計算＋アセンブリリスト出力用)
.PHONY: asmlist
asmlist:
	@rm -f *.o
	@rm -f *.s
	@cd ./c; $(MAKE) --no-print-directory asmlist
	$(CC) -S -g $(LIBHEAD) $(CFLAGS) -c $(SRC)


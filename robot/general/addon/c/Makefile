# ARCS6 サブディレクトリC言語アドオン用Makefile
# 2021/07/31

# 最終的に出来上がる中間ファイルの名前
OBJNAME = ARCS_ADDONC.o

# 追加ライブラリ中間ファイルへのパス
LIBNAME =

# ディレクトリパス
LIBPATH = ../../../../lib
SRCPATH = ../../../../src
SYSPATH = ../../../../sys

# ライブラリヘッダファイルへのパス
LIBHEAD = -I. -I$(LIBPATH) -I$(SRCPATH) -I$(SYSPATH) -I/usr/src/linux/include

# コンパイラとリンカの指定
CC = gcc
LD = ld

# コンパイラに渡すオプション
CFLAGS  = -Wall -march=native -fPIE

# リンカに渡すオプション
LDFLAGS = -L.

# Cファイル形式の指定
.SUFFIXES:.c .h .o
SRC=$(shell ls *.c)
OBJ=$(SRC:.c=.o)

# ディレクトリ内すべてコンパイル＆リンク
all: $(OBJNAME)

# ディレクトリ内すべてコンパイル＆リンク(オフライン計算用)
offline: $(OBJNAME)

# 依存関係ファイル読み込み
-include Makefile.depend

# ディレクトリ内すべてリンク
$(OBJNAME): $(OBJ)
	@$(LD) $(LDFLAGS) -r -o $(OBJNAME) $(OBJ) $(LIBNAME)

# ディレクトリ内すべてコンパイル
.c.o:
	$(CC) -O2 $(LIBHEAD) $(CFLAGS) -c $<

# 掃除
.PHONY: clean
clean:
	@rm -f *.o
	@rm -f *.s

# 依存関係ファイルの生成
.PHONY: dep
dep:
	@$(CC) $(LIBHEAD) $(CFLAGS) -MM $(SRC) > Makefile.depend

# デバッグ
.PHONY: debug
debug:
	@rm -f *.o
	@rm -f *.s
	$(CC) -ggdb3 -Og $(LIBHEAD) $(CFLAGS) -c $(SRC)
	$(LD) $(LDFLAGS) -r -o $(OBJNAME) $(OBJ) $(LIBNAME)

# すべてコンパイル (オフライン計算＋アセンブリリスト出力用)
.PHONY: asmlist
asmlist:
	@rm -f *.o
	@rm -f *.s
	$(CC) -S -g $(LIBHEAD) $(CFLAGS) -c $(SRC)


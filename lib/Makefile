# ARCS6サブディレクトリ用Makefile
# 2021/07/31

# 最終的に出来上がる中間ファイルの名前
OBJNAME = ARCS_LIBS.o

# ライブラリヘッダファイルへのパス
LIBHEAD = -I. -I../src -I../sys -I/usr/src/linux/include

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
	@$(LD) $(LDFLAGS) -r -o $(OBJNAME) $(OBJ)

# ディレクトリ内すべてコンパイル
.cc.o:
	$(CC) -O2 -pipe $(LIBHEAD) $(CFLAGS) -c $<

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
	@$(MAKE) --no-print-directory clean
	$(CC) -ggdb3 -Og $(LIBHEAD) $(CFLAGS) -c $(SRC)
	$(LD) $(LDFLAGS) -r -o $(OBJNAME) $(OBJ)

# すべてコンパイル (オフライン計算＋アセンブリリスト出力用)
.PHONY: asmlist
asmlist:
	@$(MAKE) --no-print-directory clean
	$(CC) -S -g $(LIBHEAD) $(CFLAGS) -c $(SRC)


% MultiPlot マルチ(マルタイ)プロットクラス
% 2024/04/03 Yokokura, Yuki
%
% 本プロットの特徴
% ・余白を完全に制御下に置ける多段プロット
% 　(デフォのsubplotでは余白取りすぎて論文用の図にし難い)
% ・フォント名・サイズ、線の色・太さ、など、決まり切った設定をしなくて良い
%
% 　 サイズ定義パラメータの位置関係は↓の通り
%{
    0%                                                  100%
100% +--------------------- FigWidth ----------------------+
     |       |TopMargin
     |       +        +--------------------------+
     |       |        |       PlaneNum = 2       |
     |       |        |PlaneHeight               |
     |       |        |                          |
     |       |        +--------------------------+
     |       |        +--------------------------+ <- XaxisSubMargin
     |Fig    |All     |       PlaneNum = 1       |
     |Height |Planes  |PlaneHeight               |
     |       |Height  |                          |
     |       |        +--------------------------+
     |       |        +--------------------------+ <- XaxisSubMargin
     |       |        |       PlaneNum = 0       |
     |       |        |PlaneHeight               |
     |       |        |                          |
     |       +        +------- PlaneWidth -------+
     |       |Xaxis
     |       |Margin
  0% +----------------+                          +---------+
        YaxisMargin                              RightMargin
%}
classdef MultiPlot < handle

% 定数プロパティ
properties(Constant)

end

% クラスプロパティ
properties(SetAccess = protected)
	FigHandle;				% figureハンドル
	FigLeft    = 0;			% [px] figureの左座標
	FigBottom  = 0;			% [px] figureの下座標
	FigWidth   = 0;			% [px] figureの幅
	FigHeight  = 0;			% [px] figureの高さ
	NumPlanes  = 0;			% プロット平面の段数
	CurrPlane  = 0;			% 現在選択中のプロット平面
	PlaneHandle;			% axesハンドル
	PlaneLeft  = 0;			% [px] プロット平面の左位置
	PlaneWidth = 0;			% [px] プロット平面の横幅
	PlaneHeight = 0;		% [px] プロット平面の高さ
	AllPlanesHeight = 0;	% [px] すべてのプロット平面の高さ
	YaxisMargin = 130;		% [px] Y軸用の余白
	XaxisSubMargin = 3;		% [px] 中間段のX軸用の余白
	XaxisMargin = 70;		% [px] 最下段のX軸用の余白
	TopMargin   = 10;		% [px] 上側の余白
	RightMargin = 20;		% [px] 右側の余白
	FontSize   = 18;		% [pt] 文字サイズ
	FontName   = 'Times New Roman';	% フォント名
	AxisColor  = [0, 0, 0];			% 軸の色
	BackColor  = [1, 1, 1];			% 背景色
	ColorBlack = [      0,       0,       0];	% プロット線の色(黒系)
	ColorGray  = [127/255, 127/255, 127/255];	% プロット線の色(灰系)
	ColorRed   = [225/255,   0/255,   0/255];	% プロット線の色(赤系)
	ColorGreen = [  0/255, 255/255, 127/255];	% プロット線の色(緑系)
	ColorBlue  = [  0/255, 127/255, 255/255];	% プロット線の色(青系)
	WidthNarrow = 1.2;		% プロット線の太さ(細)
	WidthNormal = 2.0;		% プロット線の太さ(普通)
	WidthWide   = 2.8;		% プロット線の太さ(太)
	LineColors;				% プロット線の色用連想配列
	LineWidths;				% プロット線の太さ用連想配列
	LegendEdgeColor = [1, 1, 1];			% 凡例の枠色
	LegendBackColor = [0.95, 0.95, 0.95];	% 凡例の背景色
	XaxisName = 'Time [s]';	% 時間軸のラベル名
	Xmin = 0;				% X軸範囲の最小値
	Xstp = 0;				% X軸のグリッド間隔
	Xmax = 0;				% X軸範囲の最大値
end

% 静的メソッド
methods(Static)
	% CSVファイルを読み込む関数
	function varargout = LoadCsvFile(FileName)
		% CSVファイル読み込み
		FileFormat = FileName(length(FileName)-3:length(FileName));	% 拡張子読み込み
		if strcmp(FileFormat, '.csv')
			RawData = readmatrix(FileName,'FileType','text');		% CSVファイルの場合
		else
			assert(false, 'MultiPlot: 未知のファイル形式');			% 未知のファイルの場合
		end
		
		% 出力変数の数のチェック
		[h, w] = size(RawData);
		if w < nargout
			assert(false, 'MultiPlot: CSVファイルの列の数 < 出力変数の数');
		end
		
		% 出力変数として返す
		for i = 1:nargout
			varargout{i} = RawData(:, i);
		end
	end
end

% メソッド
methods
	
	% コンストラクタ
	function this = MultiPlot(Handle)
		this.FigHandle = Handle;
		clf(this.FigHandle);
		set(this.FigHandle, 'PaperPositionMode','auto');
		set(this.FigHandle, 'Position',[this.FigLeft, this.FigBottom, this.FigWidth, this.FigHeight]);
		set(this.FigHandle, 'Color',this.BackColor);
		
		% プロット線色の初期化
		ColorDefs = {this.ColorBlack, this.ColorGray, this.ColorRed, this.ColorGreen, this.ColorBlue};
		ColorNames = ["Black", "Gray", "Red", "Green", "Blue"];
		this.LineColors = dictionary(ColorNames, ColorDefs);	% ←連想配列(R2022b以降が必要)
		
		% プロット線幅の初期化
		WidthDefs = {this.WidthNarrow, this.WidthNormal, this.WidthWide};
		WidthSizes = ["Thin", "Normal", "Heavy"];
		this.LineWidths = dictionary(WidthSizes, WidthDefs);	% ←連想配列(R2022b以降が必要)
	end
	
	% グラフ位置の設定関数
	function FigurePosition(this, Left, Bottom)
		this.FigLeft = Left;
		this.FigBottom = Bottom;
		set(this.FigHandle, 'Position', [this.FigLeft, this.FigBottom, this.FigWidth, this.FigHeight]);
	end
	
	% グラフサイズの設定関数
	function FigureSize(this, Width, Height)
		this.FigWidth = Width;
		this.FigHeight = Height;
		set(this.FigHandle, 'Position', [this.FigLeft, this.FigBottom, this.FigWidth, this.FigHeight]);
	end
	
	% 余白の設定関数
	function FigureMargin(this, Left, Bottom, Right)
		this.YaxisMargin = Left;	% [px] Y軸用の余白
		this.XaxisMargin = Bottom;	% [px] 最下段のX軸用の余白
		this.RightMargin = Right;	% [px] 右側の余白
	end
	
	% フォントの設定関数
	function Font(this, Name, Size)
		this.FontName = Name;
		this.FontSize = Size;
	end
	
	% プロット平面の段数の設定関数
	function NumOfPlanes(this, NumPlanes)
		% サイズ定義パラメータの計算
		this.NumPlanes = NumPlanes;
		this.AllPlanesHeight = this.FigHeight - this.XaxisMargin - this.TopMargin + this.XaxisSubMargin;% [px] すべてのプロット平面の高さ
		this.PlaneHeight = this.AllPlanesHeight/this.NumPlanes - this.XaxisSubMargin;					% [px] プロット平面の高さ
		this.PlaneWidth = this.FigWidth - this.YaxisMargin - this.RightMargin;							% [px] プロット平面の幅
		this.PlaneLeft = this.YaxisMargin;																% [px] プロット平面の左位置
		
		% 段数分のプロット平面を生成
		for i = 0:1:(NumPlanes - 1)
			this.PlaneHandle(i + 1) = this.GeneratePlane(i);
		end
	end
	
	function FontSettings(this, Ax)
		% X軸ラベルとフォントの設定
		if this.CurrPlane ~= 0
			set(Ax,'XTickLabel',{});	% 最下段以外はX軸ラベルを消す
		else
			xlabel(this.XaxisName,'FontSize',this.FontSize,'FontName',this.FontName);	% 最下段のみX軸ラベルを付ける
		end
		
		% フォントと色の設定
		set(Ax, 'FontSize',this.FontSize, 'FontName',this.FontName);
		set(Ax, 'Color',this.BackColor, 'XColor',this.AxisColor, 'YColor',this.AxisColor);
	end
	
	% プロット平面を生成する関数
	function Ax = GeneratePlane(this, PlaneNum)
		% 座標計算
		w = 1/this.FigWidth*this.PlaneWidth;	% プロット平面の幅を換算 [px] → [%]
		h = 1/this.FigHeight*this.PlaneHeight;	% プロット平面の高さを換算 [px] → [%]
		l = 1/this.FigWidth*this.PlaneLeft;		% プロット平面の左位置を換算 [px] → [%]
		b = 1/this.FigHeight*(this.AllPlanesHeight/this.NumPlanes*PlaneNum + this.XaxisMargin);		% プロット平面の下位置を換算 [px] → [%]
		Pos = [l, b, w, h];						% [左%, 下%, 幅%, 高さ%] 
		Ax = axes('Position', Pos);				% プロット平面生成
		plot(Ax, 0, 0);							% 初期化用plot() 空プロットをしておく
		
		% X軸ラベルとフォントの設定
		if PlaneNum ~= 0
			set(Ax,'XTickLabel',{});	% 最下段以外はX軸ラベルを消す
		else
			xlabel(this.XaxisName,'FontSize',this.FontSize,'FontName',this.FontName);	% 最下段のみX軸ラベルを付ける
		end
		
		% フォントと色の設定
		ylabel(Ax, 'y', 'FontSize',this.FontSize, 'FontName',this.FontName);
		set(Ax, 'FontSize',this.FontSize, 'FontName',this.FontName);
		set(Ax, 'Color',this.BackColor, 'XColor',this.AxisColor, 'YColor',this.AxisColor);
	end
	
	% プロット平面を選択する関数(何段目にプロットするかの設定)
	function SelectPlane(this, PlaneNum)
		this.CurrPlane = this.NumPlanes - PlaneNum;
		axes(this.PlaneHandle( this.NumPlanes - PlaneNum + 1 ))
	end
	
	% X軸範囲の設定関数
	function XaxisRange(this, Xmin, Xstp, Xmax)
		this.Xmin = Xmin;
		this.Xstp = Xstp;
		this.Xmax = Xmax;
		
		% 段数分のプロット平面に対してX軸範囲を設定
		for i = 1:1:this.NumPlanes
			xlim(this.PlaneHandle(i), [Xmin, Xmax]);
			set(this.PlaneHandle(i),'XTickMode','manual','XTick', Xmin:Xstp:Xmax);
		end
	end
	
	% X軸ラベル名の設定関数
	function XaxisLabel(this, Xname)
		this.XaxisName = Xname;
		xlabel(this.PlaneHandle(1), this.XaxisName,'FontSize',this.FontSize,'FontName',this.FontName);	% 最下段のX軸ラベル
	end
	
	% X/Y軸グリッドを自動設定する関数
	function AutoGrid(this, Ax)
		Ymax = max(Ax.YTick);
		Ymin = min(Ax.YTick);
		Ystp = (Ymax - Ymin)/(length(Ax.YTick) - 1);
		Xmax = max(Ax.XTick);
		Xmin = min(Ax.XTick);
		Xstp = (Xmax - Xmin)/(length(Ax.XTick) - 1);
		grid on;
		axis([Xmin, Xmax, Ymin - Ystp/2, Ymax + Ystp/2]);
		set(Ax,'XTickMode','manual','YTickMode','manual','XTick', Xmin:Xstp:Xmax,'YTick', Ymin:Ystp:Ymax);
		Ax.YAxis.Exponent = 0;
	end
	
	% Y軸グリッドを手動設定する関数
	function ManualGrid(this, Ymin, Ystp, Ymax)
		grid on;
		ylim([Ymin - Ystp/2, Ymax + Ystp/2]);
		set(gca,'YTickMode','manual','YTick', Ymin:Ystp:Ymax);
	end
	
	% 縦軸ラベルの設定関数
	function Label(this, Name)
		ylabel(gca, Name, 'FontSize',this.FontSize, 'FontName',this.FontName);
	end
	
	% 凡例の設定関数
	function Legend(this, Name, Position, Direction)
		nll = {char};		% 空文字
		Name = [nll, Name];	% 初期化用plot()の分の凡例は非表示
		h = legend(Name, 'Location',Position, 'Orientation',Direction);
		set(h, 'TextColor',this.AxisColor);
		set(h, 'FontSize',this.FontSize, 'FontName',this.FontName);
		set(h, 'EdgeColor', this.LegendEdgeColor, 'Color', this.LegendBackColor);
	end
	
	% データをダウンサンプリングして間引く関数
	function xout = DownSampling(this, xin, dsmpl)
		len = length(xin(:,1));
		xout = xin(1:dsmpl:len, :);
	end
	
	% プロット関数
	function Plot(this, varargin)
		% 引数を読み込み
		PlotFunc = varargin{1};	% 関数ハンドル
		Argin = varargin{2};	% 引数セル配列
		x = Argin{1};			% X軸データ
		y = Argin{2};			% Y軸データ
		Color = Argin{3};		% 線の色
		Width = Argin{4};		% 線の太さ
		if size(varargin) == [1, 3]
			% プロット点の種類が指定された場合
			PlotType = varargin{3};
		else
			% プロット点の種類が指定されなかった場合
			PlotType = ' ';
		end
		
		% 引数の数で動作を変更
		if length(Argin) == 4
			% ダウンサンプリングしない場合
			DownSmplRate = 1;
		else
			% ダウンサンプリングする場合
			DownSmplRate = Argin{5};
			x = this.DownSampling(x, DownSmplRate);
			y = this.DownSampling(y, DownSmplRate);
		end
		
		% 実際のプロット動作
		hold on;
		h = PlotFunc(gca, x, y, PlotType);
		c = this.LineColors(Color);
		w = this.LineWidths(Width);
		set(h, 'LineWidth', w{1}, 'Color', c{1});
		hold off;
		this.AutoGrid(gca);
	end
	
	% 階段プロット関数
	function StairsPlot(this, varargin)
		this.Plot(@stairs, varargin);
	end
	
	% 線形プロット関数
	function LinePlot(this, varargin)
		this.Plot(@plot, varargin);
	end
	
	% 点プロット関数
	function DotPlot(this, varargin)
		this.Plot(@plot, varargin, '.');
	end
	
	% ×プロット関数
	function CrossPlot(this, varargin)
		this.Plot(@plot, varargin, 'x');
	end
	
	% ○プロット関数
	function CirclePlot(this, varargin)
		this.Plot(@plot, varargin, 'o');
	end
	
	% ＋プロット関数
	function PlusPlot(this, varargin)
		this.Plot(@plot, varargin, '+');
	end
	
	% PNG画像とEPSファイル生成
	function SavePNGandEPS(this, FileName)
		PngFileName = strcat(FileName(1:length(FileName)-4), '.png');	% 拡張子をPNGに変換
		EpsFileName = strcat(FileName(1:length(FileName)-4), '.eps');	% 拡張子をEPSに変換
		exportgraphics(this.FigHandle, PngFileName,'ContentType','image', 'Resolution',150);
		exportgraphics(this.FigHandle, EpsFileName,'ContentType','vector','BackgroundColor','none');
	end
	
	% PNG画像とPDFファイル生成
	function SavePNGandPDF(this, FileName)
		PngFileName = strcat(FileName(1:length(FileName)-4), '.png');	% 拡張子をPNGに変換
		PdfFileName = strcat(FileName(1:length(FileName)-4), '.pdf');	% 拡張子をPDFに変換
		exportgraphics(this.FigHandle, PngFileName,'ContentType','image', 'Resolution',150);
		exportgraphics(this.FigHandle, PdfFileName,'ContentType','vector');
	end
end
end


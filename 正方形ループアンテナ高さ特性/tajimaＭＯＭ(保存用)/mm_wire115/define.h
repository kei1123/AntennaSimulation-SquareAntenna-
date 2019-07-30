//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//
//			変数、定数、配列宣言ヘッダー　　define.h
//
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
//				プロトタイプ宣言
//////////////////////////////////////////////////////////////////////////////////
//--初期化
void Initialization(void);						//配列と係数の初期化
//--配列
void MakeMatAll(void);							//配列確保
void DelMatAll(void);							//配列開放
//--計算
void Calculation(void);							//計算順序
//--インピーダンス行列関係
void MakeZmn(void);								//インピーダンス行列作成
void MakeGaussPara(void);						//ガウス積分のパラメータ作成
complex grmn(int ,int ,int ,int);				//gρmnの計算 (m ,n ,i ,j)
complex hij(int ,int ,int ,int ,int ,double);	//hji (<1 or 2 or 3>,m ,n ,i ,j ,xi)
//--電流分布関係
void MakeCurrent(void);							//電流分布
void MakePhase(void);							//電流位相計算
void Pibot(int);								//ピボットの選択(Zmn[i][i]≠0の判定)
void Gauss(void);								//ガウスの消去法
complex Idx(int ,int , double);					//電流計算 (ワイヤーNo 
												//		   ,ワイヤー中のセグメントNo
												//         ,位置)
//--放射界関係
void Radiation(void);							//放射界計算
complex FarFieldT(double ,double);				//遠方界θ成分 (θ，φ)方向
complex FarFieldF(double ,double);				//遠方界φ成分 (θ，φ)方向
//--入出力関係
void ConfInput(void);							//入力
void OutputCurrent(void);						//電流出力
void OutputConf(void);							//形状出力
void OutputRad(void);							//放射界　放射パターン出力
void OutputChara(void);							//特性データ出力
void FileOpen(void);							//ファイル開く
void FileClose(void);							//ファイル閉じる
void StartTime(void);							//タイマースタート
void LapsedTime(void);							//経過時間
void ConfCheck(void);							//形状検査

//////////////////////////////////////////////////////////////////////////////////
//			定数
//////////////////////////////////////////////////////////////////////////////////
//--基本物理定数
const double PI = 3.1415926535897932385;		//π
const double C  = 2.99792458*pow(10.0,8.0);			//光速
const double e0 = 8.85418780*pow(10.0,-12.0);		//真空中の透磁率ε0
const double u0 = 4.0*PI*pow(10.0,-7.0);			//真空中の誘電率μ0
const double R  = 1.0;							//1.0 固定(遠方界計算距離)
//--式簡略用
complex J=Complex(0.0,1.0);                     //虚数単位　J=√-1

//--積分用
const int GaussTenNor = 4 ;						//ガウス積分分点  [ 4固定]（通常用）
const int GaussTenSpe = 40;						//ガウス積分分点数[40固定]（特異点用）

//////////////////////////////////////////////////////////////////////////////////
//			配列の名称
//////////////////////////////////////////////////////////////////////////////////
//--形状配列
double *RX , *RY , *RZ;							//各セグメントの始点　	RX[ ]  RY[ ]  RZ[ ]
double *SX , *SY , *SZ;							//単位ベクトルの成分　	SX[ ]  SY[ ]  SZ[ ]
double *SEGL , *RA;								//各セグメントの長さ 半径	SEGL[ ]	RA[ ]　		
int *SEGN , *FEDP;								//各ワイヤーのセグメント数　SEGN[ ]
complex *FEDV;									//各給電電圧			FEDV[ ]
int *RAD_MODE;									//各セグメントの放射モード		RAD_MODE[ ]
int *LOADP;										//各インピーダンス装荷位置
complex *LOADZ;									//各インピーダンス装荷値

//--電流配列
complex **Zmn;									//インピーダンス行列	Zmn[ ][ ]
complex *Im;									//電流分布				Im[ ]
double *PhaseIm;								//電流の位相			PhaseIm[ ]

//--入力インピーダンス
complex *Zin;									//各給電点の入力インピーダンス
double *VSWR_50 , *VSWR_75;						//各給電点の50Ωと75Ωに対するVSWR

//--放射界配列
complex *TRAD ,*FRAD;							//EθとEφの放射界 TRAD[ ] FRAD[ ]
complex *RRAD ,*LRAD;							//ER とEL の放射界 RRAD[ ] LRAD[ ]
double *TPHASE ,*FPHASE;						//EθとEφの位相 TPHASE[ ] FPHASE[ ]
double *RPHASE ,*LPHASE;						//ER とEL の位相 RPHASE[ ] LPHASE[ ]
double *TGAI ,*FGAI;							//EθとEφの利得 TGAI[ ] FGAI[ ]
double *TFGAI;									//EθとEφの合計の利得 TFGAI[ ]
double *RLGAI;									//ER とEL の合計の利得 RLGAI[ ]
double *RGAI ,*LGAI;							//ER とEL の利得 RGAI[ ] LGAI[ ]
double *TPAT ,*FPAT;							//EθとEφのパターン TPAT[ ] FPAT[ ]
double *RPAT ,*LPAT;							//ER とEL のパターン RPAT[ ] LPAT[ ]
double *AR;										//軸比 AR[ ]

//--ガウス積分のパラメータ
double *GaussWeightNor , *GaussBuntenNor;		//重み 分点
double *GaussWeightSpe , *GaussBuntenSpe;		//重み 分点 特異点用

//////////////////////////////////////////////////////////////////////////////////
//			ファイルポインタ
//////////////////////////////////////////////////////////////////////////////////
FILE *fp_conf;									//形状出力用
FILE *fp_curr;									//電流分布出力用
FILE *fp_radf;									//放射界出力用
FILE *fp_radp;									//放射パターン出力用
FILE *fp_char1;									//特性出力用
FILE *fp_char2;									//特性出力用
FILE *fp_free;									//任意データ出力用

//////////////////////////////////////////////////////////////////////////////////
//			変数
//////////////////////////////////////////////////////////////////////////////////
//--形状 給電
double FREQ0;									//設計周波数f0
double RAMDA0;									//設計周波数の自由空間波長
double USEF;									//給電周波数
int NWIR;										//全ワイヤー数
int NSEG;										//全セグメント数
int NFED;										//全給電点数
int NSEG0;										//全電流計算点数
int NLOAD;										//抵抗等装荷数

//--計算用
double k0;										//k0
complex uair;									//u(air)
double Beta;									//β

//--放射界計算
static int    AXMODE   ;						//(0=φ固定, 1=θ固定)	
double DegDelta ;								//刻み幅("0.0"禁止)
double DegStart ;								//初期角
double DegWidth ;								//範囲
double FixAngle ;								//固定軸角度

//--時間測定
int		la_min,la_sec;							//分 秒
time_t	start, finish;							//開始、終了

//--計算回し用パラメータ
double PARA1,PARA2,PARA3;	//出力に使うだけ
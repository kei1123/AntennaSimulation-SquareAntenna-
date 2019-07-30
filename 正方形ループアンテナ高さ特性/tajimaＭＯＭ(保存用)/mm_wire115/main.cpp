//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//
//			mm_wire	ver 1.15
//			Method of Moments
//			Arbitrarily shaped Wire Antennas
//																	Ryoichi Tajima
//																	    2006/02/08
//
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <conio.h>
#include <time.h>
#include <direct.h>
#include "complex.h"					//複素数計算
#include "define.h"						//宣言とか ("complex.h" の後に呼び出す)
#include "makemat.h"					//配列確保 ("define.h"  の後に呼び出す)
#include "antenna.h"					//入力部分 ("makemat.h" の後に呼び出す)

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
void OutputFree(void);							//任意データ出力
void FileOpen(void);							//ファイル開く
void FileClose(void);							//ファイル閉じる
void StartTime(void);							//タイマースタート
void LapsedTime(void);							//経過時間
void ConfCheck(void);							//形状検査

//////////////////////////////////////////////////////////////////////////////////
//				メイン
//////////////////////////////////////////////////////////////////////////////////
void main(void)
{
	FileOpen();					//ファイル開く
	StartTime();				//タイマースタート
	Calculation();				//"antenna.h"の計算
	FileClose();				//ファイル閉じる
	getch();					//待った
}

//////////////////////////////////////////////////////////////////////////////////
//				形状検査
//////////////////////////////////////////////////////////////////////////////////
void ConfCheck(void)					//形状検査
{
	int i;						//カウンタ
	int f = 0;					//間違えがあると[f = 1]になる
	//----------------------------------------------------
	// 給電点数１以上？
	if(NFED == 0){
		printf("Make Conf Error !!  NFED = 0\n");
		f = 1;					//間違え
	}
	//----------------------------------------------------

	//----------------------------------------------------
	// 第１給電の電圧？
	//		第２給電がある場合でも，第２給電電圧=(0.0+j0.0)の可能性が
	//		あるので，検査は第１給電のみ
	if(Abs(FEDV[0]) == 0.0){
		printf("Make Conf Error !!  FEDV[0] = 0.0 + j0.0 [V]\n");
		f = 1;					//間違え
	}
	//----------------------------------------------------

	//----------------------------------------------------
	// セグメント長が入力されているか？
	for (i = 0; i < NSEG; ++i){
		if(SEGL[i] == 0.0){
			printf("Make Conf Error !!  SEGL[%d] = 0.0\n",i);
			f = 1;				//間違え
			i = NSEG;			//この項目の検査終了
		}
	}
	//----------------------------------------------------

	//----------------------------------------------------
	// ワイヤー半径は入力されているか？
	for (i = 0; i < NSEG; ++i){
		if(RA[i] == 0.0){
			printf("Make Conf Error !!  RA[%d] = 0.0\n",i);
			f = 1;				//間違え
			i = NSEG;			//この項目の検査終了
		}
	}
	//----------------------------------------------------

	//----------------------------------------------------
	// ワイヤー毎のセグメント数が入力されているか？
	for (i = 0; i < NWIR; ++i){
		if(SEGN[i] == 0){
			printf("Make Conf Error !!  SEGN[%d] = 0\n",i);
			f = 1;				//間違え
			i = NWIR;			//この項目の検査終了
		}
	}
	//----------------------------------------------------

	//----------------------------------------------------
	// 検査でエラーがあれば終了
	if(f==1){
		printf("\n\n");
		OutputConf();				//形状出力
		getch();					//待った
		exit( 0 );					//終わり
	}
	//----------------------------------------------------
}

//////////////////////////////////////////////////////////////////////////////////
//				初期化等
//////////////////////////////////////////////////////////////////////////////////
void Initialization(void)
{
	int i;		//カウンタ
	//----------------------------------------------------
	// 計算用の変数代入
	k0 = 2.0*PI*USEF * sqrt(e0*u0);
	uair = -1.0 * J * 2.0*PI*USEF *u0/(4 * PI * pow(k0,2));
	Beta = 2.0*PI*USEF*sqrt(e0*u0);
	//----------------------------------------------------
	
	//----------------------------------------------------	
	// 放射モードを全放射に設定
	for (i = 0; i < NSEG; ++i){
		RAD_MODE[i] = 1;
	}
	//----------------------------------------------------

	//----------------------------------------------------
	// 電流(電圧)行列初期化
	for(i = 0 ;i < NSEG0 ;++i){
		Im[i] = Complex(0.0,0.0);
	}
	//----------------------------------------------------

	//----------------------------------------------------
	// 抵抗装荷
	for(i = 0 ;i < NLOAD ;++i){
		LOADP[i] = 1;
		LOADZ[i] = Complex(0.0,0.0);
	}
	//----------------------------------------------------

	//----------------------------------------------------
	// ガウス積分のパラメータ作成
	MakeGaussPara();
	//----------------------------------------------------
}

//////////////////////////////////////////////////////////////////////////////////
//				インピーダンス行列作成
//////////////////////////////////////////////////////////////////////////////////
void MakeZmn(void)
{
	int m , n;					//セグメントNo
	int wm , wn;				//ワイヤーNo
	int cwirem , cwiren;		//ワイヤーNo
	int tempm , tempn;			//次のワイヤーの終点No

	wm = 0;
	cwirem = 0;
	tempm = 0;

	//----------------------------------------------------
	// 画面出力
	printf("Impedance Matrix                 0/%6d ",NSEG0);
	//----------------------------------------------------
	for(m = 1;m <= NSEG; ++m){						//波源側のセグメントNo
		if(m == tempm + SEGN[cwirem]){				//ｍがワイヤー(エレメント)の終端かを判断
			tempm = tempm + SEGN[cwirem];			//次のワイヤー(エレメント)の終端を計算
			cwirem = cwirem + 1;					//つぎのワイヤーNo
		}
		else{
			//----------------------------------------------------
			// 画面出力
			printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b");	//前回の画面出力を消去
			printf("%6d /%6d",wm+1,NSEG0);			//計算している行を画面出力
			//----------------------------------------------------
			wm = wm +1;								//配列行番号
			wn = 0;
			cwiren = 0;
			tempn = 0;
			for(n = 1; n <= NSEG; ++n){				//観測点のセグメントNo
				if(n == tempn + SEGN[cwiren]){		//ｎがワイヤー(エレメント)の終端かを判断
					tempn = tempn + SEGN[cwiren];	//次のワイヤー(エレメント)の終端を計算
					cwiren = cwiren + 1;			//次ののワイヤーNo
				}
				else{
					wn = wn +1;
					Zmn[wm-1][wn-1] = grmn(m, n, m-1, n-1) + grmn(m, n, m-1, n-0)
								    + grmn(m, n, m-0, n-1) + grmn(m, n, m-0, n-0);
				}
			}
		}
	}

	//----------------------------------------------------
	// インピーダンス行列操作(抵抗付加)
	//     antenna.hで指定した位置(LOADP[])にインピーダンスの値(LOADZ[])を加える
	//     加える場所は，インピーダンス行列の対角成分(自己インピーダンス成分)
	for(m=0;m<NLOAD;++m){
		Zmn[LOADP[m]][LOADP[m]].re = Zmn[LOADP[m]][LOADP[m]].re + LOADZ[m].re;
		Zmn[LOADP[m]][LOADP[m]].im = Zmn[LOADP[m]][LOADP[m]].im + LOADZ[m].im;
	}
	//----------------------------------------------------
	printf("\n");
}

complex grmn(int gm ,int gn ,int gi ,int gj)
{
	double nmi , xmi;	//ηm(i) , xm(i)
	double xi_xj;		//xi・xj  (xiとxjの内積結果)
	complex answ_h13 = Complex(0.0 , 0.0); //結果代入用
	complex answ_h2  = Complex(0.0 , 0.0); //結果代入用
	complex answer   = Complex(0.0 , 0.0); //結果代入用
	int gaussi;			//積分で使うカウンタ
	int GaussTen ;		//ガウス積分分点数代入用
	double gausst;		//一次変換結果の代入
	double temp;		//一時保管用(double)
	complex temp_a;		//一時保管用(complex)

	//----------------------------------------------------
	// 係数計算
	if(gi == gm-1)	  nmi =  1.0;								//ηm(i)
	else if(gi == gm) nmi = -1.0;								//ηm(i)
	xi_xj = SX[gi]*SX[gj]  +  SY[gi]*SY[gj]  +  SZ[gi]*SZ[gj];	//内積
	//----------------------------------------------------

	//----------------------------------------------------
	// 数値積分はノーマルで行う
	//    grmnの積分はノーマル(積分分点数=4)で行う
	GaussTen = GaussTenNor;
	//----------------------------------------------------

	//----------------------------------------------------
	// h<1>とh<3>の計算
	//    hij<3>はεr=1.0の場合0.0になるので，ここでは「0.0」を足している．
	//    hij<3>=0.0の意味なので，消さずに残す
	for(gaussi = 0; gaussi < GaussTen; ++gaussi){
		//--積分位置
		gausst = GaussBuntenNor[gaussi] * SEGL[gi]/2.0 + SEGL[gi]/2.0;
		//--係数計算
		if(gi == gm-1)	  xmi = k0 * gausst;
		else if(gi == gm) xmi = k0 * (SEGL[gm] - gausst);
		//--積分関数
		temp_a = ( hij(1,gm,gn,gi,gj,gausst) + 0.0 );	//hji<3> はεr=1.0の場合「0.0」
		temp = (cos(xmi) * (GaussWeightNor[gaussi] * SEGL[gi]/2.0));
		answ_h13.re = answ_h13.re + temp * temp_a.re;
		answ_h13.im = answ_h13.im + temp * temp_a.im;
	}
	//----------------------------------------------------

	//----------------------------------------------------
	// h<2>の計算
	for(gaussi = 0; gaussi < GaussTen; ++gaussi){
		//--積分位置
		gausst = GaussBuntenNor[gaussi] * SEGL[gi]/2.0 + SEGL[gi]/2.0;
		//--係数計算
		if(gi == gm-1)	  xmi = k0 * gausst;
		else if(gi == gm) xmi = k0 * (SEGL[gm] - gausst);
		//--積分関数
		temp_a = hij(2,gm,gn,gi,gj,gausst);
		temp = (sin(xmi) * (GaussWeightNor[gaussi] * SEGL[gi]/2.0));
		answ_h2.re = answ_h2.re + temp * temp_a.re;
		answ_h2.im = answ_h2.im + temp * temp_a.im;
	}
	//----------------------------------------------------

	answer.re = nmi * answ_h13.re  +  xi_xj * answ_h2.re;
	answer.im = nmi * answ_h13.im  +  xi_xj * answ_h2.im;

	//----------------------------------------------------
	// 論文では２重積分の外にある処理( k0/sin(d*k0) )
	//     観測点側のセグメント長SEGL[gi]を使用するため
	//     論文ではこの部分でtempの値が2乗にするが，観測点側と波源側のセグメント長
	//     に自由度を持たせるために，こんな処理(２乗をせずに，SEGL[gi]はhijの計算
	//     の最後に移す)をしている
	temp = k0/sin(k0*SEGL[gi]);
	answer.re = answer.re * temp;
	answer.im = answer.im * temp;
	//----------------------------------------------------
	return (answer);
}

complex hij(int hmode ,int gm ,int gn ,int gi,int gj,double gxi)
{
	double nnj , xnj;						//ηn(j) , xn(j)
	double rspx , rspy , rspz;				//距離代入用
	double rowij;							//波源と観測点の距離
	double sin_cos;							//係数(sin or cos)
	int gaussi;								//積分で使うカウンタ
	int GaussTen ;							//ガウス積分分点数代入用
	double seki_nor_spe;					//積分の特異点の境目(セグメント長との比)
	double gausst;							//一次変換結果の代入
	complex answer = Complex(0.0 , 0.0);	//結果代入用
	double temp ;							//一時保管用
	double temp_x , temp_y , temp_z;		//一時保管用 x,y,z
	complex temp_a , temp_b ,temp_c;		//一時保管用 複素数

	//----------------------------------------------------
	// 積分でこの値(×セグメント長)より近いと特異点
	//     この閾値(1.01)に決めたのは，同じ長さのセグメントが２つ隣以内の場合特異点
	//     処理を行うため．(近い位置のセグメント計算は特異点)
	//     この値が大きすぎると，特異点処理が多くなり，無駄計算になる．
	//     この値が小さすぎると，近いセグメント同士の計算で誤差の原因になる．
	//     セグメント長よりも少し大きいくらいがいい．だから(1.01)を使用．
	seki_nor_spe = 1.01;
	//----------------------------------------------------

	//----------------------------------------------------
	// 係数計算
	if(hmode == 1){							//h<1>を計算している場合
		if(gj == gn-1)    nnj =  1.0;
		else if(gj == gn) nnj = -1.0;
		nnj = nnj * (-1);
	}
	else if(hmode == 2)	  nnj =  1.0;		//h<2>を計算している場合
	//----------------------------------------------------

	//----------------------------------------------------
	// 数値積分の分点数を決定
	//     ２つのセグメントの端同士で近い距離が[観測点側のセグメント長×seki_nor_spe]より
	//     近い場合には特異点として処理をする．(分点数をGaussTenSpeにする)
	//     GaussTenSpe：積分の分点数を40点に増加させ，積分精度を上げている
	GaussTen = GaussTenNor;								//最初はノーマルに設定
	//--始点 vs. 始点の距離計算
	rspx = (RX[gi]) - (RX[gj]);							//x軸
	rspy = (RY[gi]) - (RY[gj]);							//y軸
	rspz = (RZ[gi]) - (RZ[gj]);							//z軸
	temp = sqrt(rspx*rspx + rspy*rspy + rspz*rspz);		//絶対値
	if(temp<=(SEGL[gi]+SEGL[gj])){
		if(temp<=SEGL[gi]*seki_nor_spe) GaussTen = GaussTenSpe;
		else if(GaussTen == GaussTenNor){
			//--始点 vs. 終点の距離計算
			rspx = (RX[gi]) - (RX[gj]+SX[gj]*SEGL[gj]);
			rspy = (RY[gi]) - (RY[gj]+SY[gj]*SEGL[gj]);
			rspz = (RZ[gi]) - (RZ[gj]+SZ[gj]*SEGL[gj]);
			temp = sqrt(rspx*rspx + rspy*rspy + rspz*rspz);
			if(temp<=SEGL[gi]*seki_nor_spe) GaussTen = GaussTenSpe;
		}
		else if(GaussTen == GaussTenNor){
			//--終点 vs. 始点の距離計算
			rspx = (RX[gi]+SX[gi]*SEGL[gi]) - (RX[gj]);
			rspy = (RY[gi]+SY[gi]*SEGL[gi]) - (RY[gj]);
			rspz = (RZ[gi]+SZ[gi]*SEGL[gi]) - (RZ[gj]);
			temp = sqrt(rspx*rspx + rspy*rspy + rspz*rspz);
			if(temp<=SEGL[gi]*seki_nor_spe) GaussTen = GaussTenSpe;
		}
		else if(GaussTen == GaussTenNor){
			//--終点 vs. 終点の距離計算
			rspx = (RX[gi]+SX[gi]*SEGL[gi]) - (RX[gj]+SX[gj]*SEGL[gj]);
			rspy = (RY[gi]+SY[gi]*SEGL[gi]) - (RY[gj]+SY[gj]*SEGL[gj]);
			rspz = (RZ[gi]+SZ[gi]*SEGL[gi]) - (RZ[gj]+SZ[gj]*SEGL[gj]);
			temp = sqrt(rspx*rspx + rspy*rspy + rspz*rspz);
			if(temp<=SEGL[gi]*seki_nor_spe) GaussTen = GaussTenSpe;
		}
	}
	//----------------------------------------------------
	
	//----------------------------------------------------
	// 数値積分
	for(gaussi = 0; gaussi < GaussTen; ++gaussi){
		if(GaussTen == GaussTenNor){
			gausst = GaussBuntenNor[gaussi] * SEGL[gj]/2.0 + SEGL[gj]/2.0;
			//--係数計算
			if(gj == gn-1)		xnj = k0 * gausst;
			else if(gj == gn)	xnj = k0 * (SEGL[gn] - gausst);
			if(hmode == 1)		sin_cos = cos(xnj);
			else if(hmode == 2) sin_cos = sin(xnj);
			//--距離計算
			temp_x = (RX[gi]+SX[gi]*gxi) - (RX[gj]+SX[gj]*gausst);
			temp_y = (RY[gi]+SY[gi]*gxi) - (RY[gj]+SY[gj]*gausst);
			temp_z = (RZ[gi]+SZ[gi]*gxi) - (RZ[gj]+SZ[gj]*gausst);
			rowij = sqrt( temp_x*temp_x + temp_y*temp_y + temp_z*temp_z + RA[gi]*RA[gi] );
			//--積分関数			
			temp =  ( sin_cos * GaussWeightNor[gaussi] * SEGL[gj]/(2.0 * rowij) ) ;
		}
		else if(GaussTen == GaussTenSpe){
			gausst = GaussBuntenSpe[gaussi] * SEGL[gj]/2.0 + SEGL[gj]/2.0;
			//--係数計算
			if(gj == gn-1)		xnj = k0 * gausst;
			else if(gj == gn)	xnj = k0 * (SEGL[gn] - gausst);
			if(hmode == 1)		sin_cos = cos(xnj);
			else if(hmode == 2) sin_cos = sin(xnj);
			//--距離計算
			temp_x = (RX[gi]+SX[gi]*gxi) - (RX[gj]+SX[gj]*gausst);
			temp_y = (RY[gi]+SY[gi]*gxi) - (RY[gj]+SY[gj]*gausst);
			temp_z = (RZ[gi]+SZ[gi]*gxi) - (RZ[gj]+SZ[gj]*gausst);
			rowij = sqrt( temp_x*temp_x + temp_y*temp_y + temp_z*temp_z + RA[gi]*RA[gi] );
			//--積分関数
			temp =  ( sin_cos * GaussWeightSpe[gaussi] * SEGL[gj]/(2.0 * rowij) ) ;
		}
		//--係数計算
		temp_c.re = 0.0;
		temp_c.im = -1.0*k0*rowij;
		temp_b.re = exp(temp_c.re)*cos(temp_c.im);
		temp_b.im = exp(temp_c.re)*sin(temp_c.im);
		temp_a.re = uair.re * temp_b.re - uair.im * temp_b.im;
		temp_a.im = uair.re * temp_b.im + uair.im * temp_b.re;
		answer.re = answer.re + temp * temp_a.re;
		answer.im = answer.im + temp * temp_a.im;
	}
	//----------------------------------------------------

	answer.re = answer.re * nnj;
	answer.im = answer.im * nnj;
	
	//----------------------------------------------------
	// 論文では２重積分の外にある処理( k0/sin(d*k0) )
	// 波源側のセグメント長SEGL[gi]を使用するため
	temp = k0/sin(k0*SEGL[gj]);
	answer.re = answer.re * temp;
	answer.im = answer.im * temp;
	//----------------------------------------------------
	return (answer);	
}


//////////////////////////////////////////////////////////////////////////////////
//				電流分布計算
//////////////////////////////////////////////////////////////////////////////////
void MakeCurrent(void)				//電流分布
{
	int i;							//カウンタ
	double fedv1_abs;				//第１給電電圧の絶対値
	complex hansya;					//反射係数
	//----------------------------------------------------
	// 画面出力
	printf("Current Distribution");
	printf("                      ");
	//----------------------------------------------------

	//----------------------------------------------------
	// 給電電圧代入	
	//	  ※電圧条件
	//         ①第１給電の電圧はイメージ法を用いる場合   -2×(1+j0)[V]
	//                          イメージ法を用いない場合 -1×(1+j0)[V]
	//         ②第１給電がイメージ法の場合，他の給電点もイメージ法を適用している．
	//           (振幅，位相にかかわらず-2倍)
	//         ③給電点数を複数にして，電圧値制御等を行う場合，第１給電は条件①
	//           のままで，第１給電以外の電圧で制御を行う．（条件①は絶対）
	//    電流行列に電圧値を入れているが，行列式計算後にIm[ ]に入っている値は電流
	for(i = 0 ;i < NFED ;++i){
		Im[ FEDP[i] ] = FEDV[i];
	}
	//----------------------------------------------------
	
	//----------------------------------------------------
	// 行列式計算
	Gauss();
	//----------------------------------------------------
	
	//----------------------------------------------------
	// 入力インピーダンス計算
	//    fedv1_absについて
	//	  	  イメージ法を用いている場合[fedv1_abs=2.0]がとなり，
	//		  イメージ法を用いていない場合[fedv1_abs=1.0]とする．
	//	  インピーダンスを計算する時，「NFEDの全て」の給電電圧を[fedv1_abs]
	//    で割ることで，イメージ法による電圧２倍に対応する．
	fedv1_abs = Abs(FEDV[0]);
	for(i = 0; i < NFED ;++i){
		Zin[i] = (-1.0*FEDV[i]/fedv1_abs) / Im[FEDP[i]];		//給電位相任意
	}
	//----------------------------------------------------

	//----------------------------------------------------
	//　各給電点の50Ωに対するVSWR
	for(i = 0; i < NFED ;++i){
		hansya = (Zin[i]-50.0)/(Zin[i]+50.0);
		VSWR_50[i] = (1.0+Abs(hansya))/(1.0-Abs(hansya));
	}
	//----------------------------------------------------

	//----------------------------------------------------
	//　各給電点の75Ωに対するVSWR
	for(i = 0; i < NFED ;++i){
		hansya = (Zin[i]-75.0)/(Zin[i]+75.0);
		VSWR_75[i] = (1.0+Abs(hansya))/(1.0-Abs(hansya));
	}
	//----------------------------------------------------
}

//////////////////////////////////////////////////////////////////////////////////
//				電流位相計算
//////////////////////////////////////////////////////////////////////////////////
void MakePhase(void)
{
	int i,w,s,ns;					//カウンタ
	int flag;						//位相計算用フラグ
	double temp1,temp2;				//位相計算用

	ns = 0;
	for(w = 0; w < NWIR; ++w){				//ワイヤーNo.
		for(s = 0; s <SEGN[w]-1; ++s){		//セグメントNo.
			PhaseIm[ns] =180.0/PI * atan2(Im[ns].im,Im[ns].re);		//位相計算
			if(s != 0){
				flag = 0;					//とりあえず初期化
				//----------------------------------------------------
				// 隣り合うセグメント(ワイヤーの切れ目以外)で位相
				// が360度以上変わる場合の補正
				// この処理を行うと，進行波の位相が進行波のように見える
				// ※±360度を行うだけなので，本質的な変化は無い
				for(i = 1 ; flag <= 0 ; ++i){		//「flag==0」の間は何度でも
					flag = 1;						//とりあえず終了フラグ
					temp1 = sqrt(pow(PhaseIm[ns]-PhaseIm[ns-1],2));
					temp2 = sqrt(pow(PhaseIm[ns]+360-PhaseIm[ns-1],2));
					if( temp1 > temp2 ){
						PhaseIm[ns] = PhaseIm[ns] + 360.0;
						flag = 0;					//もう一度
					}
					temp2 = sqrt(pow(PhaseIm[ns]-360-PhaseIm[ns-1],2));
					if( temp1 > temp2 ){
						PhaseIm[ns] = PhaseIm[ns] - 360.0;
						flag = 0;					//もう一度
					}
				}
				//----------------------------------------------------
			}
			ns = ns + 1 ;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////
//				ガウスの消去法
//////////////////////////////////////////////////////////////////////////////////
void Gauss(void)										//ガウスの消去法
{
	int m,n,j,k;										//カウンタ
	//----------------------------------------------------
	// 前進消去
	complex temp,u;
	for(m = 0 ;m < NSEG0 ;++m){
		//----------------------------------------------------
		// 画面出力
		printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b");			//前回の画面出力を消去
		printf("%6d /%6d",m+1,NSEG0);					//計算している行を画面出力
		//----------------------------------------------------
		if( (Zmn[m][m].re == 0.0) && (Zmn[m][m].im == 0.0) ){
			Pibot( m );									//ピボットの確認
		}
		u.re =  1.0*Zmn[m][m].re/(Zmn[m][m].re*Zmn[m][m].re + Zmn[m][m].im*Zmn[m][m].im);
		u.im = -1.0*Zmn[m][m].im/(Zmn[m][m].re*Zmn[m][m].re + Zmn[m][m].im*Zmn[m][m].im);
		for(n = 0 ;n < NSEG0 ;++n){
			temp = Zmn[m][n];
			Zmn[m][n].re = temp.re * u.re - temp.im * u .im;
			Zmn[m][n].im = temp.re * u.im + temp.im * u .re;
		}
		temp = Im[m];
		Im[m].re = temp.re * u.re - temp.im * u .im;
		Im[m].im = temp.re * u.im + temp.im * u .re;
		for(k = m+1 ;k < NSEG0 ;++k){
			temp = Zmn[k][m] ;							//掃出すための係数の取得
			for(j=m+1 ; j < NSEG0 ; ++j){				//対角成分以下の掃出し
				Zmn[k][j].re = Zmn[k][j].re - (Zmn[m][j].re*temp.re-Zmn[m][j].im*temp.im);
				Zmn[k][j].im = Zmn[k][j].im - (Zmn[m][j].re*temp.im+Zmn[m][j].im*temp.re);
			}
			Im[k].re = Im[k].re - (Im[m].re*temp.re - Im[m].im*temp.im);
			Im[k].im = Im[k].im - (Im[m].re*temp.im + Im[m].im*temp.re);
		}
	}
	printf("\n");
	//----------------------------------------------------

	//----------------------------------------------------
	// 後退代入
	for(m = 0 ;m < NSEG0 ;++m){
		temp = Complex(0.0,0.0);
		for(int k = NSEG0-m ;k < NSEG0 ;++k){			//m＝0(1回目)の時には入らない
			temp = temp + Zmn[NSEG0-1-m][k]*Im[k];
		}
		Im[NSEG0-1-m] = Im[NSEG0-1-m] - temp;
	}
	//----------------------------------------------------
}

void Pibot(int i)										//ピボットの選択(Zmn[i][i]≠0の判定)
{
	//----------------------------------------------------
	// 数値を｢0｣で割らないようにする
	// アンテナの形状入力等が正常ならば，必要ない部分
	// ガウスの消去法で，分母が0にならないようにしている
	complex temp;
	if( Abs(Zmn[i][i]) == 0.0 ){						//Zmn[i-1][i-1]=0の時入る
		for(int x=i+1 ;Abs(Zmn[i][i]) == 0.0 ;++x){
			//----------------------------------------------------
			// xが最後の行になると入る
			if(x>=NSEG0){								//計算続行不可能
				OutputConf();							//形状出力
				printf("\n\nDivision by Zero !!\n\n");	//分母=0になる場合(エラー)
				getch();
				exit ( 0 );								//終了
			}
			//----------------------------------------------------

			//----------------------------------------------------
			// 行の入れ替え
			else{
				for(int y = 0 ;y < NSEG0 ; ++y){
					temp = Zmn[i][y];
					Zmn[i][y] = Zmn[x][y];
					Zmn[x][y] = temp;	
				}
			}
			//----------------------------------------------------
		}
	}
	//----------------------------------------------------
}

//////////////////////////////////////////////////////////////////////////////////
//				セグメント上の任意位置の電流計算
//////////////////////////////////////////////////////////////////////////////////
complex Idx(int wn,int sn, double dx)		//電流計算 (ワイヤーNo 
											//		   ,ワイヤー中のセグメントNo
											//         ,位置)
{
	//----------------------------------------------------
	//   放射界を計算するとき，セグメントの端の電流だけで遠方界を計算すると，誤差の源
	//   プログラム中では，セグメント上のどの位置でも，電流値を算出できるようにしている．
	//----------------------------------------------------
	complex answer1 = Complex(0.0 , 0.0);	//結果代入用
	complex answer2 = Complex(0.0 , 0.0);	//結果代入用
	int INo;								//セグメント始点の電流No (最初は "0" )
	int SNo;								//計算対象のセグメントNo
	int i;									//カウンタ

	//----------------------------------------------------
	// 電流No , セグメントNo　算出
	SNo = 0;
	for(i = 0 ;i < wn; ++i){
		SNo = SNo + SEGN[i] ;
	}
	SNo = SNo + sn;
	INo = SNo - wn;
	//---------------------------------------------------

	//----------------------------------------------------
	// ワイヤーの始点側のセグメント用電流
	if(sn != SEGN[wn]-1){							//終点の時は計算しない
		answer1 = Im[INo] * sin(k0*(dx)) / sin(k0*(SEGL[SNo]));
	}
	//---------------------------------------------------

	//----------------------------------------------------
	// ワイヤーの終点側のセグメント用電流
	if(sn != 0){									//始点の時は計算しない
		answer2 = Im[INo-1] * sin(k0*(SEGL[SNo]-dx)) / sin(k0*(SEGL[SNo]));
	}
	//---------------------------------------------------

	//---------------------------------------------------
	// 位置毎に電流値を返す
	if(sn == 0)					return (answer1);	//ワイヤーの始点のセグメント用電流
	else if(sn == SEGN[wn]-1)	return (answer2);	//ワイヤーの終点のセグメント用電流
	else 	return (answer1 + answer2);				//その他の位置のセグメント用電流
	//---------------------------------------------------
}

//////////////////////////////////////////////////////////////////////////////////
//				放射界計算
//////////////////////////////////////////////////////////////////////////////////
void Radiation(void)
{
	double temp;
	double i;		//計算角度
	int n;			//カウンタ
	complex Zi0;	//第1給電の入力インピーダンス
	complex Ii0;	//第1給電部の電流

	//----------------------------------------------------
	// 利得計算用係数
	//   第1給電のみ利得計算に使用
	//     ※複数給電は利得補正必要
	//       自動で利得補正を行う事も出来るが，自動化すると汎用性をそこなう可能性が有る
	//       ので，このプログラムでは，利得計算時に第１給電点の電力のみ用い，後に利得補
	//		 正を行う
	Ii0 = Im[ FEDP[0] ];				//第１給電点の電流
	Zi0 = Zin[0];						//第１給電点の入力インピーダンス
	//----------------------------------------------------

	i = 0.0;
	n = 0;
	
	//----------------------------------------------------
	// 画面出力
	printf("Radiation Pattern               %10.2f",i);
	//----------------------------------------------------

	for(i = DegStart ;i <= DegStart+DegWidth ;i = i + DegDelta){
		//----------------------------------------------------
		// 画面表示
		printf("\b\b\b\b\b\b\b\b\b\b");
		printf("%10.2f",i);
		//----------------------------------------------------

		//----------------------------------------------------
		// 直線偏波の遠方界
		if(AXMODE == 0){				//φ固定の場合
			TRAD[n] = FarFieldT(i , FixAngle);
			FRAD[n] = FarFieldF(i , FixAngle);
		}
		else if(AXMODE == 1){			//θ固定の場合
			TRAD[n] = FarFieldT(FixAngle , i);
			FRAD[n] = FarFieldF(FixAngle , i);
		}
		FPHASE[n] =180.0/PI * atan2(FRAD[n].im,FRAD[n].re);		//位相
		TPHASE[n] =180.0/PI * atan2(TRAD[n].im,TRAD[n].re);		//位相
		//----------------------------------------------------

		//----------------------------------------------------
		// 円偏波の遠方界
		LRAD[n] = 0.5 * (TRAD[n]-FRAD[n]*J) ;					//左旋円偏波
		RRAD[n] = 0.5 * (TRAD[n]+FRAD[n]*J) ;					//右旋円偏波
		RPHASE[n] =180.0/PI * atan2(RRAD[n].im,RRAD[n].re);		//位相
		LPHASE[n] =180.0/PI * atan2(LRAD[n].im,LRAD[n].re);		//位相
		//----------------------------------------------------

		//----------------------------------------------------
		// 軸比の計算
		temp = (Abs(RRAD[n])+Abs(LRAD[n]))/ (Abs(RRAD[n])-Abs(LRAD[n]));
		temp = sqrt(temp*temp);									//絶対値(逆旋回に対応)
		if(temp >= 1000) temp = 1000;							//軸比の最大値は60dB
		AR[n] = 20*log10(temp);		
		//----------------------------------------------------

		//----------------------------------------------------
		// 直線偏波θの利得計算
		temp = pow(Abs(TRAD[n]),2)*pow(R,2) / (30.0*pow(Abs(Ii0),2)*Zi0.re) ;
		if(temp <= 0.000001) temp = 0.000001;					//利得最小値は-60dB
		TGAI[n] = 10*log10( temp );
		//----------------------------------------------------
		
		//----------------------------------------------------
		// 直線偏波φの利得計算
		temp = pow(Abs(FRAD[n]),2)*pow(R,2) / (30.0*pow(Abs(Ii0),2)*Zi0.re) ;
		if(temp <= 0.000001) temp = 0.000001;					//利得最小値は-60dB
		FGAI[n] = 10*log10( temp );
		//----------------------------------------------------
		
		//----------------------------------------------------
		// 直線偏波θとφの利得計算
		temp = (pow(Abs(TRAD[n]),2)+pow(Abs(FRAD[n]),2))*pow(R,2)
			 / (30.0*pow(Abs(Ii0),2)*Zi0.re) ;
		if(temp <= 0.000001) temp = 0.000001;					//利得最小値は-60dB
		TFGAI[n] = 10*log10( temp );
		//----------------------------------------------------

		//----------------------------------------------------
		// 円偏波Rの利得計算
		temp = 2.0*pow(Abs(RRAD[n]),2)*pow(R,2) / (30.0*pow(Abs(Ii0),2)*Zi0.re) ;
		if(temp <= 0.000001) temp = 0.000001;					//利得最小値は-60dB
		RGAI[n] = 10*log10( temp );
		//----------------------------------------------------
		
		//----------------------------------------------------
		// 円偏波Lの利得計算
		temp = 2.0*pow(Abs(LRAD[n]),2)*pow(R,2) / (30.0*pow(Abs(Ii0),2)*Zi0.re) ;
		if(temp <= 0.000001) temp = 0.000001;					//利得最小値は-60dB
		LGAI[n] = 10*log10( temp );
		//----------------------------------------------------
		
		//----------------------------------------------------
		// 円偏波RとLの利得計算
		temp = (2.0*pow(Abs(RRAD[n]),2)+2.0*pow(Abs(LRAD[n]),2))*pow(R,2)
			 / (30.0*pow(Abs(Ii0),2)*Zi0.re) ;
		if(temp <= 0.000001) temp = 0.000001;					//利得最小値は-60dB
		RLGAI[n] = 10*log10( temp );
		//----------------------------------------------------
		n = n + 1;
	}
	//----------------------------------------------------
	// 直線偏波の放射パターン計算
	n = 0;
	temp = TGAI[0];
	for(i = DegStart ;i <= DegStart+DegWidth ;i = i + DegDelta){
		if(temp < TGAI[n])	temp = TGAI[n];		//最大利得値探し
		if(temp < FGAI[n])	temp = FGAI[n];		//最大利得値探し
		n = n + 1;	
	}
	n = 0;
	for(i = DegStart ;i <= DegStart+DegWidth ;i = i + DegDelta){
		TPAT[n] = TGAI[n] - temp;				//[利得]-[最大利得]
		FPAT[n] = FGAI[n] - temp;				//[利得]-[最大利得]
		if(TPAT[n] < -30.0) TPAT[n] = -30.0;	//パターン最小値は-30dB
		if(FPAT[n] < -30.0) FPAT[n] = -30.0;	//パターン最小値は-30dB
		n = n + 1;	
	}
	//----------------------------------------------------

	//----------------------------------------------------
	// 円偏波の放射パターン計算
	n = 0;
	temp = RGAI[0];
	for(i = DegStart ;i <= DegStart+DegWidth ;i = i + DegDelta){
		if(temp < RGAI[n])	temp = RGAI[n];		//最大利得値探し
		if(temp < LGAI[n])	temp = LGAI[n];		//最大利得値探し
		n = n + 1;	
	}
	n = 0;
	for(i = DegStart ;i <= DegStart+DegWidth ;i = i + DegDelta){
		RPAT[n] = RGAI[n] - temp;				//[利得]-[最大利得]
		LPAT[n] = LGAI[n] - temp;				//[利得]-[最大利得]
		if(RPAT[n] < -30.0) RPAT[n] = -30.0;	//パターン最小値は-30dB
		if(LPAT[n] < -30.0) LPAT[n] = -30.0;	//パターン最小値は-30dB
		n = n + 1;	
	}
	//----------------------------------------------------
	printf("\n");
}

complex FarFieldT(double theta,double fai)		//遠方界θ成分 (θ，φ)方向
{
	int w,s;				//ワイヤー　セグメントのカウンタ
	int ns;					//セグメントナンバー
	complex Ith;			//Iθ
	double delR;			//ΔR
	double radth , radfa;	//θとφのラジアン値
	int gaussi;				//積分で使うカウンタ
	double gausst;			//一次変換結果の代入
	complex answer;			//結果代入用
	complex answer_s;		//各セグメントの結果代入用

	//----------------------------------------------------
	// deg -> rad の変換
	radth = theta * PI / 180.0;
	radfa =   fai * PI / 180.0;
	//----------------------------------------------------

	answer = Complex(0.0 , 0.0);	//結果代入用
	ns = 0;
	for(w = 0; w < NWIR; ++w){
		for(s = 0; s < SEGN[w]; ++s){
			answer_s = Complex(0.0 , 0.0);
			for(gaussi = 0; gaussi < GaussTenNor; ++gaussi){
				gausst = GaussBuntenNor[gaussi] * SEGL[ns]/2.0 + SEGL[ns]/2.0;
				//--Iθ (電流のθ方向成分)
				Ith = Idx(w,s,gausst)
						*(  SX[ns] * cos(radth) * cos(radfa)
						  + SY[ns] * cos(radth) * sin(radfa)
						  - SZ[ns] * sin(radth)   );
				//--ΔR
				delR = (RX[ns]+SX[ns]*gausst) * sin(radth) * cos(radfa)
				     + (RY[ns]+SY[ns]*gausst) * sin(radth) * sin(radfa)
					 + (RZ[ns]+SZ[ns]*gausst) * cos(radth);
				//--非積分関数
				answer_s = answer_s
					     + Ith * Exp(J*(Beta*delR)) 
						       * (GaussWeightNor[gaussi] * SEGL[ns]/2.0);
			}
			
			//----------------------------------------------------
			// 放射モードを適用
			//    RAD_MODE[ ]が「0」の場合，ここで放射を止める
			//    RAD_MODE[ ]をdouble型で宣言すれば，セグメント毎に放射率を指定する事も可能
			//    アンテナの解析で使わないのでRAD_MODE[ ]をint型にしている(0％または100％)
			answer_s = answer_s * RAD_MODE[ns];
			//----------------------------------------------------
			answer = answer + answer_s;
			ns = ns + 1 ;
		}
	}
	answer = answer * ((-1.0)*2.0*PI*USEF*u0/(4.0*PI*R))*(Exp((-1.0*Beta*R)*J)*J);
	return(answer);
}

complex FarFieldF(double theta,double fai)				//遠方界φ成分 (θ，φ)方向
{
	int w,s;				//ワイヤー　セグメントのカウンタ
	int ns;					//セグメントナンバー
	complex Ifa;			//Iφ
	double delR;			//ΔR
	double radth , radfa;	//ラジアン

	int gaussi;				//積分で使うカウンタ
	double gausst;			//一次変換結果の代入
	complex answer;			//結果代入用
	complex answer_s;		//各セグメントの結果代入用

	//----------------------------------------------------
	// deg -> rad の変換
	radth = theta * PI / 180.0;
	radfa =   fai * PI / 180.0;
	//----------------------------------------------------

	answer = Complex(0.0 , 0.0);	//結果代入用
	ns = 0;
	for(w = 0; w < NWIR; ++w){
		for(s = 0; s <SEGN[w]; ++s){
			answer_s = Complex(0.0 , 0.0);
			for(gaussi = 0; gaussi < GaussTenNor; ++gaussi){
				gausst = GaussBuntenNor[gaussi] * SEGL[ns]/2.0 + SEGL[ns]/2.0;
				//--Iφ (電流のφ方向成分)
				Ifa = Idx(w,s,gausst)
						*(  SX[ns] * sin(radfa) * (-1.0)
						  + SY[ns] * cos(radfa)   );
				//--ΔR
				delR = (RX[ns]+SX[ns]*gausst) * sin(radth) * cos(radfa)
				     + (RY[ns]+SY[ns]*gausst) * sin(radth) * sin(radfa)
					 + (RZ[ns]+SZ[ns]*gausst) * cos(radth);
				//--非積分関数
				answer_s = answer_s
					     + Ifa * Exp(J*(Beta*delR)) 
						       * (GaussWeightNor[gaussi] * SEGL[ns]/2.0);
			}
			//----------------------------------------------------
			// 放射モードを適用
			//    RAD_MODE[ ]が「0」の場合，ここで放射を止める
			//    RAD_MODE[ ]をdouble型で宣言すれば，セグメント毎に放射率を指定する事も可能
			//    アンテナの解析で使わないのでRAD_MODE[ ]をint型にしている(0％または100％)
			answer_s = answer_s * RAD_MODE[ns];
			//----------------------------------------------------
			answer = answer + answer_s;
			ns = ns + 1 ;
		}
	}
	answer = answer * ((-1.0)*2.0*PI*USEF*u0/(4.0*PI*R))*(Exp((-1.0*Beta*R)*J)*J);
	return(answer);
}

//////////////////////////////////////////////////////////////////////////////////
//				ファイル開く
//////////////////////////////////////////////////////////////////////////////////
void FileOpen(void)
{
	_mkdir("AntennaData") ;											//データの出力ファイル
	//----------------------------------------------------
	// 形状用ファイルを開く ※書き換え禁止
	//    ファイル名や中身を書き換えると，形状を見るプログ
	//    ラム[mm_conf.exe]が使えなくなるので，書き換え禁止
	if((fp_conf=fopen("AntennaData\\antenna_conf.dat","w"))==NULL){		//エラーの場合｢NULL｣
		printf("write open error [ antenna_conf.dat ] !!\n\n");			//オープンエラー
		getch();
		exit(0);														//エラーの場合止める
	}
	fprintf(fp_conf," [antenna_conf.dat]\n\n");
	//----------------------------------------------------

	//----------------------------------------------------
	// 電流用ファイルを開く
	if((fp_curr=fopen("AntennaData\\antenna_current.dat","w"))==NULL){	//エラーの場合｢NULL｣
		printf("write open error [ antenna_current.dat ] !!\n\n");		//オープンエラー
		getch();
		exit(0);														//エラーの場合止める
	}
	fprintf(fp_curr," [antenna_current.dat]\n\n");
	//----------------------------------------------------

	//----------------------------------------------------
	// 放射界用ファイルを開く
	if((fp_radf=fopen("AntennaData\\antenna_field.dat","w"))==NULL){	//エラーの場合｢NULL｣
		printf("write open error [ antenna_field.dat ] !!\n\n");		//オープンエラー
		getch();
		exit(0);														//エラーの場合止める
	}
	fprintf(fp_radf," [antenna_field.dat]\n\n");
	//----------------------------------------------------

	//----------------------------------------------------
	// 放射パターン用ファイルを開く
	if((fp_radp=fopen("AntennaData\\antenna_pattern.dat","w"))==NULL){	//エラーの場合｢NULL｣
		printf("write open error [ antenna_pattern.dat ] !!\n\n");		//オープンエラー
		getch();
		exit(0);														//エラーの場合止める
	}
	fprintf(fp_radp," [antenna_pattern.dat]\n\n");
	//----------------------------------------------------

	//----------------------------------------------------
	// 特性データ出力1用ファイルを開く
	if((fp_char1=fopen("AntennaData\\antenna_chara1.dat","w"))==NULL){	//エラーの場合｢NULL｣
		printf("write open error [ antenna_chara1.dat ] !!\n\n");		//オープンエラー
		getch();
		exit(0);														//エラーの場合止める
	}
	fprintf(fp_char1," [antenna_chara1.dat]\n\n");
	fprintf(fp_char1," >>---------PARA ---------<< **");
	fprintf(fp_char1," >>------- Impedance -------- ------<< **");
	fprintf(fp_char1," >>------ Max-Direction Gain(dB) ----<< **");
	fprintf(fp_char1," >>---- 0(deg)-Gain(dB) ----<<\n");

	fprintf(fp_char1,"         para1         para2         para3 **");
	fprintf(fp_char1,"    Zin.re    Zin.im  VSWR-50  VSWR-75 **");
	fprintf(fp_char1,"    Angle       Eθ       Eφ   Eθ+Eφ **");
	fprintf(fp_char1,"       Eθ       Eφ   Eθ+Eφ\n");
	//----------------------------------------------------

	//----------------------------------------------------
	// 特性データ出力2用ファイルを開く
	if((fp_char2=fopen("AntennaData\\antenna_chara2.dat","w"))==NULL){	//エラーの場合｢NULL｣
		printf("write open error [ antenna_chara2.dat ] !!\n\n");		//オープンエラー
		getch();
		exit(0);														//エラーの場合止める
	}
	fprintf(fp_char2," [antenna_chara2.dat]\n\n");
	fprintf(fp_char2," >>------------------PARA ---------------<< **");
	fprintf(fp_char2," >>------- Impedance -------- ------<< **");
	fprintf(fp_char2," >>-------    Max Direction Gain&AR(dB) -----<< **");
	fprintf(fp_char2," >>------  0(deg) Gain&AR(dB)  -----<<\n");

	fprintf(fp_char2,"         para1         para2         para3 **");
	fprintf(fp_char2,"    Zin.re    Zin.im  VSWR-50  VSWR-75 **");
	fprintf(fp_char2,"      Angle        Er        El       Er+El      AR **");
	fprintf(fp_char2,"         Er         El        Er+El       AR\n");
	//----------------------------------------------------

	//----------------------------------------------------
	// 任意データ出力用ファイルを開く
	//   何を出力しても問題ないが，出力作業は"antenna.h"で行う
	if((fp_free=fopen("AntennaData\\antenna_free.dat","w"))==NULL){		//エラーの場合｢NULL｣
		printf("write open error [ antenna_free.dat ] !!\n\n");			//オープンエラー

		getch();
		exit(0);														//エラーの場合止める
	}
	fprintf(fp_free," [antenna_free.dat]\n\n");
	//----------------------------------------------------
}

//////////////////////////////////////////////////////////////////////////////////
//				ファイル閉じる
//////////////////////////////////////////////////////////////////////////////////
void FileClose(void)
{
	fclose(fp_curr);		//電流分布
	fclose(fp_conf);		//形状
	fclose(fp_radf);		//放射界
	fclose(fp_radp);		//パターン
	fclose(fp_char1);		//特性１
	fclose(fp_char2);		//特性２
	fclose(fp_free);		//適当出力
}

//////////////////////////////////////////////////////////////////////////////////
//			特性データ出力	antenna_chara1.dat  antenna_chara2.dat
//////////////////////////////////////////////////////////////////////////////////
void OutputChara(void)							//特性データ出力
{
	double gaimax;		//最大利得代入用
	double i;			//計算角度
	int n;				//配列No
	int maxang,zang;	//最大方向No，Z軸方向No

	//----------------------------------------------------
	//	直線偏波成分 antenna_chara1
	//--最大放射方向探し

	n = 0;//初期化
	maxang = 0;//初期化
	gaimax = TGAI[maxang];
	for(i = DegStart ;i <= DegStart+DegWidth ;i = i + DegDelta){
		if(TGAI[n]>=gaimax){
			gaimax = TGAI[n];
			maxang = n;
		}
		if(FGAI[n]>=gaimax){
			gaimax = FGAI[n];
			maxang = n;
		}
		n = n + 1;
	}

	//--Z軸方向探し(0度方向探し)θ固定の場合X軸方向探し
	n = 0;
	for(i = DegStart ;i <= DegStart+DegWidth ;i = i + DegDelta){
		if( (int)(i*10000000)==0.0){
			zang = n;
		}
		n = n + 1;
	}

	//--特性出力
	//    PARA1,PARA2
	//      特性出力で出力しているPARA1,PARA2は周波数特性などを計算した場合に表
	//      やグラフを作りやすくするためのモノ計算中では何も使用していないので，
	//      出力しなくても問題ない場合は消さない
	//    表示桁数
	//      アンテナの研究を行う上で，必要以上の少数以下の桁数を出力していない．
	//      これ以上の桁数が必要ならば，任意出力を用いる
	fprintf(fp_char1," %13.6f %13.6f %13.6f **",PARA1,PARA2,PARA3);
	fprintf(fp_char1," %9.2f %9.2f",Zin[0].re,Zin[0].im);
	fprintf(fp_char1," %8.3f %8.3f **",VSWR_50[0],VSWR_75[0]);
	fprintf(fp_char1," %8.3f",DegStart+DegDelta*maxang);
	fprintf(fp_char1," %9.3f %9.3f",TGAI[maxang],FGAI[maxang]);
	fprintf(fp_char1," %9.3f **",TFGAI[maxang]);
	fprintf(fp_char1," %9.3f %9.3f",TGAI[zang],FGAI[zang]);
	fprintf(fp_char1," %9.3f\n",TFGAI[zang]);
	//----------------------------------------------------

	//----------------------------------------------------
	//	円偏波成分 antenna_chara2
	//--最大放射方向探し
	n = 0;
	maxang = 0;
	gaimax = RGAI[maxang];
	for(i = DegStart ;i <= DegStart+DegWidth ;i = i + DegDelta){
		if(RGAI[n]>=gaimax){
			gaimax = RGAI[n];
			maxang = n;
		}
		if(LGAI[n]>=gaimax){
			gaimax = LGAI[n];
			maxang = n;
		}
		n = n + 1;
	}

	//--Z軸方向探し(0度方向探し)θ固定の場合X軸方向探し
	n = 0;
	for(i = DegStart ;i <= DegStart+DegWidth ;i = i + DegDelta){
		if( (int)(i*10000000)==0.0){
			zang = n;
		}
		n = n + 1;
	}

	//--特性出力
	//    PARA1,PARA2
	//      特性出力で出力しているPARA1,PARA2は周波数特性などを計算した場合に表
	//      やグラフを作りやすくするためのモノ計算中では何も使用していないので，
	//      出力しなくても問題ない場合は消さない
	//    表示桁数
	//      アンテナの研究を行う上で，必要以上の少数以下の桁数を出力していない．
	//      これ以上の桁数が必要ならば，任意出力を用いる
	fprintf(fp_char2," %13.6f %13.6f %13.6f **",PARA1,PARA2,PARA3);
	fprintf(fp_char2," %9.2f %9.2f",Zin[0].re,Zin[0].im);
	fprintf(fp_char2," %8.3f %8.3f **",VSWR_50[0],VSWR_75[0]);
	fprintf(fp_char2," %8.3f",DegStart+DegDelta*maxang);
	fprintf(fp_char2," %9.3f %9.3f",RGAI[maxang],LGAI[maxang]);
	fprintf(fp_char2," %9.3f %7.3f **",RLGAI[maxang],AR[maxang]);
	fprintf(fp_char2," %9.3f %9.3f",RGAI[zang],LGAI[zang]);
	fprintf(fp_char2," %9.3f %7.3f\n",RLGAI[zang],AR[zang]);
	//----------------------------------------------------
}

//////////////////////////////////////////////////////////////////////////////////
//			形状出力	antenna_conf.dat
//			********  書き換え不可  ********
//				形状を見るプログラム[mm_conf.exe]を使わない場合は書き換え可
//////////////////////////////////////////////////////////////////////////////////
void OutputConf(void)
{
	int i;			//カウンタ
	//----------------------------------------------------
	// 初期情報
	fprintf(fp_conf,"      NSEG      NWIR      NFED     NLOAD");
	fprintf(fp_conf,"           FREQ0          RAMDA0            USEF\n");
	fprintf(fp_conf," %9d",NSEG);
	fprintf(fp_conf," %9d",NWIR);
	fprintf(fp_conf," %9d",NFED);
	fprintf(fp_conf," %9d",NLOAD);
	fprintf(fp_conf," %15.0f",FREQ0);
	fprintf(fp_conf," %15.8f",RAMDA0);
	fprintf(fp_conf," %15.0f\n",USEF);
	//----------------------------------------------------
	
	//----------------------------------------------------
	// セグメント毎のデータ
	//     位置や長さを波長単位で出力する
	//       アンテナの設計が波長単位だから
	//       形状を見るときに，周波数に依存しないで見やすい
	fprintf(fp_conf,"\nSegment\n");
	fprintf(fp_conf,"   NO.");
	fprintf(fp_conf,"           SX           SY           SZ");
	fprintf(fp_conf,"             RX             RY             RZ");
	fprintf(fp_conf,"           SEGL             RA     RAD_MODE\n");
	for(i = 0; i < NSEG; ++i){
		fprintf(fp_conf," %5d",i);
		fprintf(fp_conf," %12.7f %12.7f %12.7f",SX[i],SY[i],SZ[i]);
		fprintf(fp_conf," %14.6f %14.6f %14.6f",RX[i]/RAMDA0,RY[i]/RAMDA0,RZ[i]/RAMDA0);
		fprintf(fp_conf," %14.6f",SEGL[i]/RAMDA0);
		fprintf(fp_conf," %14.6f %12d",RA[i]/RAMDA0,RAD_MODE[i]);
		fprintf(fp_conf,"\n");
	}
	//----------------------------------------------------

	//----------------------------------------------------
	// ワイヤー毎のデータ
	fprintf(fp_conf,"\nWire\n");
	fprintf(fp_conf,"   NO.");
	fprintf(fp_conf,"     SEGN\n");
	for(i = 0; i < NWIR; ++i){
		fprintf(fp_conf," %5d",i);
		fprintf(fp_conf," %8d",SEGN[i]);
		fprintf(fp_conf,"\n");
	}
	//----------------------------------------------------

	//----------------------------------------------------
	// 給電毎のデータ
	fprintf(fp_conf,"\nFeed\n");
	fprintf(fp_conf,"   NO.");
	fprintf(fp_conf,"     FEDP      FEDV.re      FEDV.im\n");
	for(i = 0; i < NFED; ++i){
		fprintf(fp_conf," %5d",i);
		fprintf(fp_conf," %8d %12.5f %12.5f",FEDP[i],FEDV[i].re,FEDV[i].im);
		fprintf(fp_conf,"\n");
	}
	//----------------------------------------------------

	//----------------------------------------------------
	// 抵抗毎のデータ
	fprintf(fp_conf,"\nLoad\n");
	fprintf(fp_conf,"   NO.");
	fprintf(fp_conf,"     LOADP      LOADZ.re      LOADZ.im\n");
	for(i = 0; i < NLOAD; ++i){
		fprintf(fp_conf," %5d",i);
		fprintf(fp_conf," %8d %12.5f %12.5f",LOADP[i],LOADZ[i].re,LOADZ[i].im);
		fprintf(fp_conf,"\n");
	}
	fprintf(fp_conf,"\n\n\n");
	//----------------------------------------------------
}

//////////////////////////////////////////////////////////////////////////////////
//			電流分布Iの出力		antenna_current.dat
//////////////////////////////////////////////////////////////////////////////////
void OutputCurrent(void)
{
	int i;			//カウンタ
	fprintf(fp_curr,"     i");
	fprintf(fp_curr,"     Im[i].re     Im[i].im");
	fprintf(fp_curr,"      |Im[i]|   PhaseIm[i] \n");
	for(i = 0; i < NSEG0; ++i){
		fprintf(fp_curr," %5d",i+1);
		fprintf(fp_curr," %12.8f %12.8f",Im[i].re ,Im[i].im);
		fprintf(fp_curr," %12.8f %12.2f",Abs(Im[i]) ,PhaseIm[i]);
		fprintf(fp_curr,"\n");
	}
	fprintf(fp_curr,"\n\n\n");
}

//////////////////////////////////////////////////////////////////////////////////
//			放射界　放射パターン出力	antenna_field.dat   antenna_pattern.dat
//////////////////////////////////////////////////////////////////////////////////
void OutputRad(void)
{
	double i,delp;			//角度，位相差
	int n;					//カウンタ

	//========================================================
	//			放射界		antenna_field.dat
	//========================================================

	//----------------------------------------------------
	// 直線偏波関係
	if(AXMODE == 0)fprintf(fp_radf,"  φ = %.3f °plane\n",FixAngle);
	if(AXMODE == 1)fprintf(fp_radf,"  θ = %.3f °plane\n",FixAngle);
	fprintf(fp_radf,"   --    ");
	fprintf(fp_radf,">>-------   -------   --Eφ--   -------   -------<<");
	fprintf(fp_radf," ** ");
	fprintf(fp_radf,">>-------   -------   --Eθ--   -------   -------<<");
	fprintf(fp_radf," ** ");
	fprintf(fp_radf,">>-------   --------<<\n");
	if(AXMODE == 0)fprintf(fp_radf," θ(deg)");
	if(AXMODE == 1)fprintf(fp_radf," φ(deg)");
	fprintf(fp_radf,"   FRAD.re   FRAD.im    |FRAD|");
	fprintf(fp_radf,"    FPHASE  FGAI(dB)");
	fprintf(fp_radf,"   **");
	fprintf(fp_radf,"   TRAD.re   TRAD.im    |TRAD|");
	fprintf(fp_radf,"    TPHASE  TGAI(dB)");
	fprintf(fp_radf,"   **");
	fprintf(fp_radf," TFGAI(dB)    TPHASE\n");

	n = 0;
	for(i = DegStart ;i <= DegStart+DegWidth ;i = i + DegDelta){
		delp = FPHASE[n]-TPHASE[n];
		fprintf(fp_radf," %7.3f",i);
		fprintf(fp_radf," %9.5f %9.5f %9.5f",FRAD[n].re,FRAD[n].im,Abs(FRAD[n]) );
		fprintf(fp_radf," %9.3f %9.4f",FPHASE[n],FGAI[n]);
		fprintf(fp_radf,"   **");
		fprintf(fp_radf," %9.5f %9.5f %9.5f",TRAD[n].re,TRAD[n].im,Abs(TRAD[n]) );
		fprintf(fp_radf," %9.3f %9.4f",TPHASE[n],TGAI[n]);
		fprintf(fp_radf,"   **");
		fprintf(fp_radf," %9.4f  %9.4f",TFGAI[n],delp);
		fprintf(fp_radf," \n");
		n = n + 1;
	}
	//----------------------------------------------------

	//----------------------------------------------------
	// 円偏波関係
	fprintf(fp_radf,"\n\n   --    ");
	fprintf(fp_radf,">>-------   -------   --ER --   -------   -------<<");
	fprintf(fp_radf," ** ");
	fprintf(fp_radf,">>-------   -------   --EL --   -------   -------<<");
	fprintf(fp_radf," ** ");
	fprintf(fp_radf,">>-------   --------<<\n");
	if(AXMODE == 0)fprintf(fp_radf," θ(deg)");
	if(AXMODE == 1)fprintf(fp_radf," φ(deg)");
	fprintf(fp_radf,"   RRAD.re   RRAD.im    |RRAD|");
	fprintf(fp_radf,"    RPHASE  RGAI(dB)");
	fprintf(fp_radf,"   **");
	fprintf(fp_radf,"   LRAD.re   LRAD.im    |LRAD|");
	fprintf(fp_radf,"    LPHASE  LGAI(dB)");
	fprintf(fp_radf,"   **");
	fprintf(fp_radf," TFGAI(dB)     AR(dB)\n");

	n = 0;
	for(i = DegStart ;i <= DegStart+DegWidth ;i = i + DegDelta){
		fprintf(fp_radf," %7.3f",i);
		fprintf(fp_radf," %9.5f %9.5f %9.5f",RRAD[n].re,RRAD[n].im,Abs(RRAD[n]) );
		fprintf(fp_radf," %9.3f %9.4f",RPHASE[n],RGAI[n]);
		fprintf(fp_radf,"   **");
		fprintf(fp_radf," %9.5f %9.5f %9.5f",LRAD[n].re,LRAD[n].im,Abs(LRAD[n]) );
		fprintf(fp_radf," %9.3f %9.4f",LPHASE[n],LGAI[n]);
		fprintf(fp_radf,"   **");
		fprintf(fp_radf," %9.4f  %9.4f",RLGAI[n],AR[n] );
		fprintf(fp_radf," \n");
		n = n + 1;
	}
	fprintf(fp_radf,"\n\n\n\n");
	//----------------------------------------------------

	//========================================================
	//			放射パターン	antenna_pattern.dat
	//========================================================
	if(AXMODE == 0)fprintf(fp_radp,"  φ = %.3f °plane\n",FixAngle);
	if(AXMODE == 1)fprintf(fp_radp,"  θ = %.3f °plane\n",FixAngle);
	if(AXMODE == 0)fprintf(fp_radp,"   θ(deg)");
	if(AXMODE == 1)fprintf(fp_radp,"   φ(deg)");
	fprintf(fp_radp,"       Eφ       Eθ       ER        EL \n");
	n = 0;
	for(i = DegStart ;i <= DegStart+DegWidth ;i = i + DegDelta){
		fprintf(fp_radp," %9.3f",i);
		fprintf(fp_radp," %9.4f %9.4f %9.4f %9.4f",FPAT[n],TPAT[n],RPAT[n],LPAT[n]);
		fprintf(fp_radp,"\n");
		n = n + 1;
	}
	fprintf(fp_radp,"\n\n\n\n");
}

//////////////////////////////////////////////////////////////////////////////////
//				計算時間測定
//////////////////////////////////////////////////////////////////////////////////
void StartTime(void)						//タイマースタート
{
	time( &start );							//タイマースタート
}

void   LapsedTime(void)                     //経過時間
{
	double d_time;
	time( &finish );						//時間取得
	d_time = difftime( finish, start );		//時間差計算
	la_sec=(int)(d_time);					//秒単位で代入
	la_min=(int)(la_sec/60);				//分単位で代入
	la_sec=la_sec%60;						//[秒単位の時間÷60]の余り=秒
	printf("Lapsed Time                         %3d:%2d\n",la_min,la_sec);
}

//////////////////////////////////////////////////////////////////////////////////
//				配列の確保
//					"makemat.h"の中にある関数を呼び出し配列の確保を行う．
//////////////////////////////////////////////////////////////////////////////////
void MakeMatAll(void)						//配列確保
{
	int NN;									//放射界の計算点数代入用

	//--形状配列
	RX = MakeMat1D_double( NSEG );			//double RX[NSEG]
	RY = MakeMat1D_double( NSEG );			//double RY[NSEG]
	RZ = MakeMat1D_double( NSEG );			//double RZ[NSEG]
	SX = MakeMat1D_double( NSEG );			//double SX[NSEG]
	SY = MakeMat1D_double( NSEG );			//double SY[NSEG]
	SZ = MakeMat1D_double( NSEG );			//double SZ[NSEG]
	SEGL = MakeMat1D_double( NSEG );		//double SEGL[NSEG]
	RA = MakeMat1D_double( NSEG );			//double RA[NSEG]	
	SEGN = MakeMat1D_int( NWIR );			//int SEGN[NWIR]
	FEDP = MakeMat1D_int( NFED );			//int FEDP[NFED]
	FEDV = MakeMat1D_complex( NFED );		//complex FEDV[NFED]
	RAD_MODE = MakeMat1D_int( NSEG );		//int RAD_MODE[ ]
	LOADP = MakeMat1D_int( NLOAD );			//各インピーダンス装荷位置
	LOADZ = MakeMat1D_complex( NLOAD );		//各インピーダンス装荷値

	//--電流配列
	Zmn = MakeMat2D_complex( NSEG0 );		//complex Zmn[NSEG0][NSEG0]
	Im = MakeMat1D_complex( NSEG0 );		//complex Im[NSEG0]
	PhaseIm = MakeMat1D_double( NSEG0 );	//double PhaseIm[NSEG0]
	Zin = MakeMat1D_complex( NFED );		//各給電点の入力インピーダンス
	VSWR_50 = MakeMat1D_double( NFED );		//各給電点の50Ωに対するVSWR
	VSWR_75 = MakeMat1D_double( NFED );		//各給電点の75Ωに対するVSWR

	//--放射界計算用
	NN = int(DegWidth/DegDelta) + 1;		//放射界の計算点数(+1は予備)
	TRAD = MakeMat1D_complex( NN );			//Eθの放射界 TRAD[  ]
	FRAD = MakeMat1D_complex( NN );			//Eφの放射界 FRAD[  ]
	RRAD = MakeMat1D_complex( NN );			//ER の放射界 RRAD[  ]
	LRAD = MakeMat1D_complex( NN );			//EL の放射界 LRAD[  ]
	TPHASE = MakeMat1D_double( NN );		//Eθの位相 TPHASE[  ]
	FPHASE = MakeMat1D_double( NN );		//Eφの位相 FPHASE[  ]
	RPHASE = MakeMat1D_double( NN );		//ER の位相 RPHASE[  ]
	LPHASE = MakeMat1D_double( NN );		//EL の位相 LPHASE[  ]
	TGAI = MakeMat1D_double( NN );			//Eθの利得 TGAI[  ]
	FGAI = MakeMat1D_double( NN );			//Eφの利得 FGAI[  ]
	TFGAI = MakeMat1D_double( NN );			//EθとEφの合計の利得 TFGAI[  ]
	RLGAI = MakeMat1D_double( NN );			//ER とEL の合計の利得 RLGAI[  ]
	RGAI = MakeMat1D_double( NN );			//ER の利得 RGAI[  ]
	LGAI = MakeMat1D_double( NN );			//EL の利得 LGAI[  ]
	TPAT = MakeMat1D_double( NN );			//Eθのパターン TPAT[  ]
	FPAT = MakeMat1D_double( NN );			//Eφのパターン FPAT[  ]
	RPAT = MakeMat1D_double( NN );			//ER のパターン RPAT[  ]
	LPAT = MakeMat1D_double( NN );			//EL のパターン LPAT[  ]
	AR = MakeMat1D_double( NN );			//軸比 AR[ ]

	//--ガウス積分のパラメータ
	GaussWeightNor = MakeMat1D_double( GaussTenNor );	//通常用
	GaussBuntenNor = MakeMat1D_double( GaussTenNor );	//通常用
	GaussWeightSpe = MakeMat1D_double( GaussTenSpe );	//特異点用
	GaussBuntenSpe = MakeMat1D_double( GaussTenSpe );	//特異点用
}

//////////////////////////////////////////////////////////////////////////////////
//				配列の開放
//////////////////////////////////////////////////////////////////////////////////
void DelMatAll(void)			//配列開放
{
	//--形状配列
	DelMat1D_double(RX);			//double RX[NSEG]
	DelMat1D_double(RY);			//double RY[NSEG]
	DelMat1D_double(RZ);			//double RZ[NSEG]
	DelMat1D_double(SX);			//double SX[NSEG]
	DelMat1D_double(SY);			//double SY[NSEG]
	DelMat1D_double(SZ);			//double SZ[NSEG]
	DelMat1D_double(SEGL);			//double SEGL[NSEG]
	DelMat1D_double(RA);			//double RA[NWIR]
	DelMat1D_int(SEGN);				//int SEGN[NWIR]
	DelMat1D_int(FEDP);				//int FEDP[NFED]
	DelMat1D_complex(FEDV);			//complex FEDV[NFED]
	DelMat1D_int(RAD_MODE);			//int RAD_MODE[NSEG]
	DelMat1D_int( LOADP );			//各インピーダンス装荷位置
	DelMat1D_complex( LOADZ );		//各インピーダンス装荷値

	//--電流配列
	DelMat2D_complex(Zmn,NSEG0);	//complex Zmn[NSEG0][NSEG0]
	DelMat1D_complex(Im);			//complex Im[NSEG0]
	DelMat1D_double(PhaseIm);		//double PhaseIm[NSEG0]
	DelMat1D_complex(Zin);			//double Zin[NFED]
	DelMat1D_double(VSWR_50);		//double VSWR_50[NFED]
	DelMat1D_double(VSWR_75);		//double VSWR_75[NFED]

	//--放射界計算用
	DelMat1D_complex(TRAD);			//Eθの放射界 TRAD[  ]
	DelMat1D_complex(FRAD);			//Eφの放射界 FRAD[  ]
	DelMat1D_complex(RRAD);			//ER の放射界 RRAD[  ]
	DelMat1D_complex(LRAD);			//EL の放射界 LRAD[  ]
	DelMat1D_double(TPHASE);		//Eθの位相 TPHASE[  ]
	DelMat1D_double(FPHASE);		//Eφの位相 FPHASE[  ]
	DelMat1D_double(RPHASE);		//ER の位相 RPHASE[  ]
	DelMat1D_double(LPHASE);		//EL の位相 LPHASE[  ]
	DelMat1D_double(TGAI);			//Eθの利得 TGAI[  ]
	DelMat1D_double(FGAI);			//Eφの利得 FGAI[  ]
	DelMat1D_double(TFGAI);			//EθとEφの合計の利得 TFGAI[  ]
	DelMat1D_double(RLGAI);			//ER とEL の合計の利得 RLGAI[  ]
	DelMat1D_double(RGAI);			//ER の利得 RGAI[  ]
	DelMat1D_double(LGAI);			//EL の利得 LGAI[  ]
	DelMat1D_double(TPAT);			//Eθのパターン TPAT[  ]
	DelMat1D_double(FPAT);			//Eφのパターン FPAT[  ]
	DelMat1D_double(RPAT);			//ER のパターン RPAT[  ]
	DelMat1D_double(LPAT);			//EL のパターン LPAT[  ]
	DelMat1D_double(AR);			//軸比 AR[ ]

	//--ガウス積分のパラメータ
	DelMat1D_double(GaussWeightNor);
	DelMat1D_double(GaussBuntenNor);
	DelMat1D_double(GaussWeightSpe);
	DelMat1D_double(GaussBuntenSpe);
}

//////////////////////////////////////////////////////////////////////////////////
//				ガウス積分のパラメータ作成
//////////////////////////////////////////////////////////////////////////////////
void MakeGaussPara(void)
{
	//----------------------------------------------------
	// 通常用
	GaussWeightNor[ 0] =  0.347854845137454;
	GaussWeightNor[ 1] =  0.652145154862546;
	GaussWeightNor[ 2] =  0.652145154862546;
	GaussWeightNor[ 3] =  0.347854845137454;

	GaussBuntenNor[ 0] = -0.861136311594053;
	GaussBuntenNor[ 1] = -0.339981043584856;
	GaussBuntenNor[ 2] =  0.339981043584856;
	GaussBuntenNor[ 3] =  0.861136311594053;
	//----------------------------------------------------

	//----------------------------------------------------
	// 特異点用
	GaussWeightSpe[ 0] =  0.001282051282051;
	GaussWeightSpe[ 1] =  0.007891011588601;
	GaussWeightSpe[ 2] =  0.014159307549920;
	GaussWeightSpe[ 3] =  0.020334759063387;
	GaussWeightSpe[ 4] =  0.026381190653141;
	GaussWeightSpe[ 5] =  0.032260717927117;
	GaussWeightSpe[ 6] =  0.037936243700708;
	GaussWeightSpe[ 7] =  0.043371908194758;
	GaussWeightSpe[ 8] =  0.048533353845914;
	GaussWeightSpe[ 9] =  0.053387951971494;
	GaussWeightSpe[10] =  0.057905011981786;
	GaussWeightSpe[11] =  0.062055976475710;
	GaussWeightSpe[12] =  0.065814602222896;
	GaussWeightSpe[13] =  0.069157126276081;
	GaussWeightSpe[14] =  0.072062416302054;
	GaussWeightSpe[15] =  0.074512104235389;
	GaussWeightSpe[16] =  0.076490702433397;
	GaussWeightSpe[17] =  0.077985701608681;
	GaussWeightSpe[18] =  0.078987649925364;
	GaussWeightSpe[19] =  0.079490212761550;
	GaussWeightSpe[20] =  0.079490212761550;
	GaussWeightSpe[21] =  0.078987649925364;
	GaussWeightSpe[22] =  0.077985701608681;
	GaussWeightSpe[23] =  0.076490702433397;
	GaussWeightSpe[24] =  0.074512104235389;
	GaussWeightSpe[25] =  0.072062416302054;
	GaussWeightSpe[26] =  0.069157126276081;
	GaussWeightSpe[27] =  0.065814602222896;
	GaussWeightSpe[28] =  0.062055976475710;
	GaussWeightSpe[29] =  0.057905011981786;
	GaussWeightSpe[30] =  0.053387951971494;
	GaussWeightSpe[31] =  0.048533353845914;
	GaussWeightSpe[32] =  0.043371908194758;
	GaussWeightSpe[33] =  0.037936243700708;
	GaussWeightSpe[34] =  0.032260717927117;
	GaussWeightSpe[35] =  0.026381190653141;
	GaussWeightSpe[36] =  0.020334759063387;
	GaussWeightSpe[37] =  0.014159307549920;
	GaussWeightSpe[38] =  0.007891011588601;
	GaussWeightSpe[39] =  0.001282051282051;

	GaussBuntenSpe[ 0] = -1.000000000000000;
	GaussBuntenSpe[ 1] = -0.995297929244349;
	GaussBuntenSpe[ 2] = -0.984266280717503;
	GaussBuntenSpe[ 3] = -0.967010076487989;
	GaussBuntenSpe[ 4] = -0.943639764943602;
	GaussBuntenSpe[ 5] = -0.914303339690209;
	GaussBuntenSpe[ 6] = -0.879186343479340;
	GaussBuntenSpe[ 7] = -0.838510822778106;
	GaussBuntenSpe[ 8] = -0.792533952601552;
	GaussBuntenSpe[ 9] = -0.741546419147384;
	GaussBuntenSpe[10] = -0.685870585084314;
	GaussBuntenSpe[11] = -0.625858452755258;
	GaussBuntenSpe[12] = -0.561889439294723;
	GaussBuntenSpe[13] = -0.494367978125254;
	GaussBuntenSpe[14] = -0.423720962155551;
	GaussBuntenSpe[15] = -0.350395044914181;
	GaussBuntenSpe[16] = -0.274853816714324;
	GaussBuntenSpe[17] = -0.197574873718911;
	GaussBuntenSpe[18] = -0.119046798444971;
	GaussBuntenSpe[19] = -0.039766070802182;
	GaussBuntenSpe[20] =  0.039766070802182;
	GaussBuntenSpe[21] =  0.119046798444971;
	GaussBuntenSpe[22] =  0.197574873718911;
	GaussBuntenSpe[23] =  0.274853816714324;
	GaussBuntenSpe[24] =  0.350395044914181;
	GaussBuntenSpe[25] =  0.423720962155551;
	GaussBuntenSpe[26] =  0.494367978125254;
	GaussBuntenSpe[27] =  0.561889439294723;
	GaussBuntenSpe[28] =  0.625858452755258;
	GaussBuntenSpe[29] =  0.685870585084314;
	GaussBuntenSpe[30] =  0.741546419147384;
	GaussBuntenSpe[31] =  0.792533952601552;
	GaussBuntenSpe[32] =  0.838510822778106;
	GaussBuntenSpe[33] =  0.879186343479340;
	GaussBuntenSpe[34] =  0.914303339690209;
	GaussBuntenSpe[35] =  0.943639764943602;
	GaussBuntenSpe[36] =  0.967010076487989;
	GaussBuntenSpe[37] =  0.984266280717503;
	GaussBuntenSpe[38] =  0.995297929244349;
	GaussBuntenSpe[39] =  1.000000000000000;
	//----------------------------------------------------
}

//============================================================================
//============================================================================
//		履歴
//
//		ver 0_00	2Dを3D用に書き換え
//		ver 0_01	特異点を40点で積分
//		ver 0_02	ρの計算を見直し
//		ver 0_03	sin(dk0)の計算を2重積分の中に移動
//		ver 0_04	初期設定関数を追加
//		ver 0_05	分離パターンを追加
//		ver 0_06	if文を訂正計算時間を短くする
//					timerを追加
//					pow関数を取り除き計算時間を短くする
//		ver 0_07	形状出力に給電点と放射モードを追加
//		ver 0_08	複素数計算を展開して時間短縮
//		ver 0_09	出力ファイル名を変える
//					antenna_chara2.datを追加
//		ver 0_10	周波数特性をPARA1で行う
//		ver 0_11	ファイル名をmomawaに変更
//		ver 0_12	hij()内の式をif文で切り替えに変更
//		ver 0_13	sin(dk0)の計算を観測点と波源のそれぞれに分けて使用
//		ver 0_14	出力ファイル(char1,2)を追加
//
//		ver 1_00	一部宣言をヘッダー化
//					計算順序関数Calculation( )を追加
//		ver 1_01	積分無駄計算取り除き
//		ver 1_02	antenna.hのみで形状入力
//		ver 1_03	距離の閾値計算見直し
//		ver 1_04	形状入力方法を変更
//		ver 1_05	形状入力をミスしても形状を出力させる
//		ver 1_06	出力書式を計算回しに対応
//		ver 1_07	入力インピーダンス計算の給電位相考慮
//		ver 1_08	放射界計算のワイヤー端(配列外参照)修正
//					形状検査に給電関係追加
//		ver 1_09	任意出力をCalculation( )の中に移動
//		ver 1_10	入力インピーダンス計算修正
//		ver 1_11	インピーダンス行列操作(抵抗ロード)を修正
//		ver 1_12	形状出力に抵抗ロードを追加
//		ver 1_13	ガウスの消去法 無駄計算取り除き
//		ver 1_14	入力インピーダンス計算修正(第２給電点以降にイメージ法適用)
//					MakeGaussPara()の呼び出し位置を変更
//					Betaの計算を初期値計算に追加
//					直線偏波の位相差の無駄な部分削除
//		ver 1_15	[Division by Zero]でも形状出力
//					コメント追加
//					利得計算とインピーダンスの付加に最近少し疑問を感じる??
//
//============================================================================
//============================================================================





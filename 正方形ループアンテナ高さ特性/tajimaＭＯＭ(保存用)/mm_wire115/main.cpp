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
#include "complex.h"					//���f���v�Z
#include "define.h"						//�錾�Ƃ� ("complex.h" �̌�ɌĂяo��)
#include "makemat.h"					//�z��m�� ("define.h"  �̌�ɌĂяo��)
#include "antenna.h"					//���͕��� ("makemat.h" �̌�ɌĂяo��)

//////////////////////////////////////////////////////////////////////////////////
//				�v���g�^�C�v�錾
//////////////////////////////////////////////////////////////////////////////////
//--������
void Initialization(void);						//�z��ƌW���̏�����
//--�z��
void MakeMatAll(void);							//�z��m��
void DelMatAll(void);							//�z��J��
//--�v�Z
void Calculation(void);							//�v�Z����
//--�C���s�[�_���X�s��֌W
void MakeZmn(void);								//�C���s�[�_���X�s��쐬
void MakeGaussPara(void);						//�K�E�X�ϕ��̃p�����[�^�쐬
complex grmn(int ,int ,int ,int);				//g��mn�̌v�Z (m ,n ,i ,j)
complex hij(int ,int ,int ,int ,int ,double);	//hji (<1 or 2 or 3>,m ,n ,i ,j ,xi)
//--�d�����z�֌W
void MakeCurrent(void);							//�d�����z
void MakePhase(void);							//�d���ʑ��v�Z
void Pibot(int);								//�s�{�b�g�̑I��(Zmn[i][i]��0�̔���)
void Gauss(void);								//�K�E�X�̏����@
complex Idx(int ,int , double);					//�d���v�Z (���C���[No 
												//		   ,���C���[���̃Z�O�����gNo
												//         ,�ʒu)
//--���ˊE�֌W
void Radiation(void);							//���ˊE�v�Z
complex FarFieldT(double ,double);				//�����E�Ɛ��� (�ƁC��)����
complex FarFieldF(double ,double);				//�����E�Ӑ��� (�ƁC��)����
//--���o�͊֌W
void ConfInput(void);							//����
void OutputCurrent(void);						//�d���o��
void OutputConf(void);							//�`��o��
void OutputRad(void);							//���ˊE�@���˃p�^�[���o��
void OutputChara(void);							//�����f�[�^�o��
void OutputFree(void);							//�C�Ӄf�[�^�o��
void FileOpen(void);							//�t�@�C���J��
void FileClose(void);							//�t�@�C������
void StartTime(void);							//�^�C�}�[�X�^�[�g
void LapsedTime(void);							//�o�ߎ���
void ConfCheck(void);							//�`�󌟍�

//////////////////////////////////////////////////////////////////////////////////
//				���C��
//////////////////////////////////////////////////////////////////////////////////
void main(void)
{
	FileOpen();					//�t�@�C���J��
	StartTime();				//�^�C�}�[�X�^�[�g
	Calculation();				//"antenna.h"�̌v�Z
	FileClose();				//�t�@�C������
	getch();					//�҂���
}

//////////////////////////////////////////////////////////////////////////////////
//				�`�󌟍�
//////////////////////////////////////////////////////////////////////////////////
void ConfCheck(void)					//�`�󌟍�
{
	int i;						//�J�E���^
	int f = 0;					//�ԈႦ�������[f = 1]�ɂȂ�
	//----------------------------------------------------
	// ���d�_���P�ȏ�H
	if(NFED == 0){
		printf("Make Conf Error !!  NFED = 0\n");
		f = 1;					//�ԈႦ
	}
	//----------------------------------------------------

	//----------------------------------------------------
	// ��P���d�̓d���H
	//		��Q���d������ꍇ�ł��C��Q���d�d��=(0.0+j0.0)�̉\����
	//		����̂ŁC�����͑�P���d�̂�
	if(Abs(FEDV[0]) == 0.0){
		printf("Make Conf Error !!  FEDV[0] = 0.0 + j0.0 [V]\n");
		f = 1;					//�ԈႦ
	}
	//----------------------------------------------------

	//----------------------------------------------------
	// �Z�O�����g�������͂���Ă��邩�H
	for (i = 0; i < NSEG; ++i){
		if(SEGL[i] == 0.0){
			printf("Make Conf Error !!  SEGL[%d] = 0.0\n",i);
			f = 1;				//�ԈႦ
			i = NSEG;			//���̍��ڂ̌����I��
		}
	}
	//----------------------------------------------------

	//----------------------------------------------------
	// ���C���[���a�͓��͂���Ă��邩�H
	for (i = 0; i < NSEG; ++i){
		if(RA[i] == 0.0){
			printf("Make Conf Error !!  RA[%d] = 0.0\n",i);
			f = 1;				//�ԈႦ
			i = NSEG;			//���̍��ڂ̌����I��
		}
	}
	//----------------------------------------------------

	//----------------------------------------------------
	// ���C���[���̃Z�O�����g�������͂���Ă��邩�H
	for (i = 0; i < NWIR; ++i){
		if(SEGN[i] == 0){
			printf("Make Conf Error !!  SEGN[%d] = 0\n",i);
			f = 1;				//�ԈႦ
			i = NWIR;			//���̍��ڂ̌����I��
		}
	}
	//----------------------------------------------------

	//----------------------------------------------------
	// �����ŃG���[������ΏI��
	if(f==1){
		printf("\n\n");
		OutputConf();				//�`��o��
		getch();					//�҂���
		exit( 0 );					//�I���
	}
	//----------------------------------------------------
}

//////////////////////////////////////////////////////////////////////////////////
//				��������
//////////////////////////////////////////////////////////////////////////////////
void Initialization(void)
{
	int i;		//�J�E���^
	//----------------------------------------------------
	// �v�Z�p�̕ϐ����
	k0 = 2.0*PI*USEF * sqrt(e0*u0);
	uair = -1.0 * J * 2.0*PI*USEF *u0/(4 * PI * pow(k0,2));
	Beta = 2.0*PI*USEF*sqrt(e0*u0);
	//----------------------------------------------------
	
	//----------------------------------------------------	
	// ���˃��[�h��S���˂ɐݒ�
	for (i = 0; i < NSEG; ++i){
		RAD_MODE[i] = 1;
	}
	//----------------------------------------------------

	//----------------------------------------------------
	// �d��(�d��)�s�񏉊���
	for(i = 0 ;i < NSEG0 ;++i){
		Im[i] = Complex(0.0,0.0);
	}
	//----------------------------------------------------

	//----------------------------------------------------
	// ��R����
	for(i = 0 ;i < NLOAD ;++i){
		LOADP[i] = 1;
		LOADZ[i] = Complex(0.0,0.0);
	}
	//----------------------------------------------------

	//----------------------------------------------------
	// �K�E�X�ϕ��̃p�����[�^�쐬
	MakeGaussPara();
	//----------------------------------------------------
}

//////////////////////////////////////////////////////////////////////////////////
//				�C���s�[�_���X�s��쐬
//////////////////////////////////////////////////////////////////////////////////
void MakeZmn(void)
{
	int m , n;					//�Z�O�����gNo
	int wm , wn;				//���C���[No
	int cwirem , cwiren;		//���C���[No
	int tempm , tempn;			//���̃��C���[�̏I�_No

	wm = 0;
	cwirem = 0;
	tempm = 0;

	//----------------------------------------------------
	// ��ʏo��
	printf("Impedance Matrix                 0/%6d ",NSEG0);
	//----------------------------------------------------
	for(m = 1;m <= NSEG; ++m){						//�g�����̃Z�O�����gNo
		if(m == tempm + SEGN[cwirem]){				//�������C���[(�G�������g)�̏I�[���𔻒f
			tempm = tempm + SEGN[cwirem];			//���̃��C���[(�G�������g)�̏I�[���v�Z
			cwirem = cwirem + 1;					//���̃��C���[No
		}
		else{
			//----------------------------------------------------
			// ��ʏo��
			printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b");	//�O��̉�ʏo�͂�����
			printf("%6d /%6d",wm+1,NSEG0);			//�v�Z���Ă���s����ʏo��
			//----------------------------------------------------
			wm = wm +1;								//�z��s�ԍ�
			wn = 0;
			cwiren = 0;
			tempn = 0;
			for(n = 1; n <= NSEG; ++n){				//�ϑ��_�̃Z�O�����gNo
				if(n == tempn + SEGN[cwiren]){		//�������C���[(�G�������g)�̏I�[���𔻒f
					tempn = tempn + SEGN[cwiren];	//���̃��C���[(�G�������g)�̏I�[���v�Z
					cwiren = cwiren + 1;			//���̂̃��C���[No
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
	// �C���s�[�_���X�s�񑀍�(��R�t��)
	//     antenna.h�Ŏw�肵���ʒu(LOADP[])�ɃC���s�[�_���X�̒l(LOADZ[])��������
	//     ������ꏊ�́C�C���s�[�_���X�s��̑Ίp����(���ȃC���s�[�_���X����)
	for(m=0;m<NLOAD;++m){
		Zmn[LOADP[m]][LOADP[m]].re = Zmn[LOADP[m]][LOADP[m]].re + LOADZ[m].re;
		Zmn[LOADP[m]][LOADP[m]].im = Zmn[LOADP[m]][LOADP[m]].im + LOADZ[m].im;
	}
	//----------------------------------------------------
	printf("\n");
}

complex grmn(int gm ,int gn ,int gi ,int gj)
{
	double nmi , xmi;	//��m(i) , xm(i)
	double xi_xj;		//xi�Exj  (xi��xj�̓��ό���)
	complex answ_h13 = Complex(0.0 , 0.0); //���ʑ���p
	complex answ_h2  = Complex(0.0 , 0.0); //���ʑ���p
	complex answer   = Complex(0.0 , 0.0); //���ʑ���p
	int gaussi;			//�ϕ��Ŏg���J�E���^
	int GaussTen ;		//�K�E�X�ϕ����_������p
	double gausst;		//�ꎟ�ϊ����ʂ̑��
	double temp;		//�ꎞ�ۊǗp(double)
	complex temp_a;		//�ꎞ�ۊǗp(complex)

	//----------------------------------------------------
	// �W���v�Z
	if(gi == gm-1)	  nmi =  1.0;								//��m(i)
	else if(gi == gm) nmi = -1.0;								//��m(i)
	xi_xj = SX[gi]*SX[gj]  +  SY[gi]*SY[gj]  +  SZ[gi]*SZ[gj];	//����
	//----------------------------------------------------

	//----------------------------------------------------
	// ���l�ϕ��̓m�[�}���ōs��
	//    grmn�̐ϕ��̓m�[�}��(�ϕ����_��=4)�ōs��
	GaussTen = GaussTenNor;
	//----------------------------------------------------

	//----------------------------------------------------
	// h<1>��h<3>�̌v�Z
	//    hij<3>�̓�r=1.0�̏ꍇ0.0�ɂȂ�̂ŁC�����ł́u0.0�v�𑫂��Ă���D
	//    hij<3>=0.0�̈Ӗ��Ȃ̂ŁC�������Ɏc��
	for(gaussi = 0; gaussi < GaussTen; ++gaussi){
		//--�ϕ��ʒu
		gausst = GaussBuntenNor[gaussi] * SEGL[gi]/2.0 + SEGL[gi]/2.0;
		//--�W���v�Z
		if(gi == gm-1)	  xmi = k0 * gausst;
		else if(gi == gm) xmi = k0 * (SEGL[gm] - gausst);
		//--�ϕ��֐�
		temp_a = ( hij(1,gm,gn,gi,gj,gausst) + 0.0 );	//hji<3> �̓�r=1.0�̏ꍇ�u0.0�v
		temp = (cos(xmi) * (GaussWeightNor[gaussi] * SEGL[gi]/2.0));
		answ_h13.re = answ_h13.re + temp * temp_a.re;
		answ_h13.im = answ_h13.im + temp * temp_a.im;
	}
	//----------------------------------------------------

	//----------------------------------------------------
	// h<2>�̌v�Z
	for(gaussi = 0; gaussi < GaussTen; ++gaussi){
		//--�ϕ��ʒu
		gausst = GaussBuntenNor[gaussi] * SEGL[gi]/2.0 + SEGL[gi]/2.0;
		//--�W���v�Z
		if(gi == gm-1)	  xmi = k0 * gausst;
		else if(gi == gm) xmi = k0 * (SEGL[gm] - gausst);
		//--�ϕ��֐�
		temp_a = hij(2,gm,gn,gi,gj,gausst);
		temp = (sin(xmi) * (GaussWeightNor[gaussi] * SEGL[gi]/2.0));
		answ_h2.re = answ_h2.re + temp * temp_a.re;
		answ_h2.im = answ_h2.im + temp * temp_a.im;
	}
	//----------------------------------------------------

	answer.re = nmi * answ_h13.re  +  xi_xj * answ_h2.re;
	answer.im = nmi * answ_h13.im  +  xi_xj * answ_h2.im;

	//----------------------------------------------------
	// �_���ł͂Q�d�ϕ��̊O�ɂ��鏈��( k0/sin(d*k0) )
	//     �ϑ��_���̃Z�O�����g��SEGL[gi]���g�p���邽��
	//     �_���ł͂��̕�����temp�̒l��2��ɂ��邪�C�ϑ��_���Ɣg�����̃Z�O�����g��
	//     �Ɏ��R�x���������邽�߂ɁC����ȏ���(�Q��������ɁCSEGL[gi]��hij�̌v�Z
	//     �̍Ō�Ɉڂ�)�����Ă���
	temp = k0/sin(k0*SEGL[gi]);
	answer.re = answer.re * temp;
	answer.im = answer.im * temp;
	//----------------------------------------------------
	return (answer);
}

complex hij(int hmode ,int gm ,int gn ,int gi,int gj,double gxi)
{
	double nnj , xnj;						//��n(j) , xn(j)
	double rspx , rspy , rspz;				//��������p
	double rowij;							//�g���Ɗϑ��_�̋���
	double sin_cos;							//�W��(sin or cos)
	int gaussi;								//�ϕ��Ŏg���J�E���^
	int GaussTen ;							//�K�E�X�ϕ����_������p
	double seki_nor_spe;					//�ϕ��̓��ٓ_�̋���(�Z�O�����g���Ƃ̔�)
	double gausst;							//�ꎟ�ϊ����ʂ̑��
	complex answer = Complex(0.0 , 0.0);	//���ʑ���p
	double temp ;							//�ꎞ�ۊǗp
	double temp_x , temp_y , temp_z;		//�ꎞ�ۊǗp x,y,z
	complex temp_a , temp_b ,temp_c;		//�ꎞ�ۊǗp ���f��

	//----------------------------------------------------
	// �ϕ��ł��̒l(�~�Z�O�����g��)���߂��Ɠ��ٓ_
	//     ����臒l(1.01)�Ɍ��߂��̂́C���������̃Z�O�����g���Q�׈ȓ��̏ꍇ���ٓ_
	//     �������s�����߁D(�߂��ʒu�̃Z�O�����g�v�Z�͓��ٓ_)
	//     ���̒l���傫������ƁC���ٓ_�����������Ȃ�C���ʌv�Z�ɂȂ�D
	//     ���̒l������������ƁC�߂��Z�O�����g���m�̌v�Z�Ō덷�̌����ɂȂ�D
	//     �Z�O�����g�����������傫�����炢�������D������(1.01)���g�p�D
	seki_nor_spe = 1.01;
	//----------------------------------------------------

	//----------------------------------------------------
	// �W���v�Z
	if(hmode == 1){							//h<1>���v�Z���Ă���ꍇ
		if(gj == gn-1)    nnj =  1.0;
		else if(gj == gn) nnj = -1.0;
		nnj = nnj * (-1);
	}
	else if(hmode == 2)	  nnj =  1.0;		//h<2>���v�Z���Ă���ꍇ
	//----------------------------------------------------

	//----------------------------------------------------
	// ���l�ϕ��̕��_��������
	//     �Q�̃Z�O�����g�̒[���m�ŋ߂�������[�ϑ��_���̃Z�O�����g���~seki_nor_spe]���
	//     �߂��ꍇ�ɂ͓��ٓ_�Ƃ��ď���������D(���_����GaussTenSpe�ɂ���)
	//     GaussTenSpe�F�ϕ��̕��_����40�_�ɑ��������C�ϕ����x���グ�Ă���
	GaussTen = GaussTenNor;								//�ŏ��̓m�[�}���ɐݒ�
	//--�n�_ vs. �n�_�̋����v�Z
	rspx = (RX[gi]) - (RX[gj]);							//x��
	rspy = (RY[gi]) - (RY[gj]);							//y��
	rspz = (RZ[gi]) - (RZ[gj]);							//z��
	temp = sqrt(rspx*rspx + rspy*rspy + rspz*rspz);		//��Βl
	if(temp<=(SEGL[gi]+SEGL[gj])){
		if(temp<=SEGL[gi]*seki_nor_spe) GaussTen = GaussTenSpe;
		else if(GaussTen == GaussTenNor){
			//--�n�_ vs. �I�_�̋����v�Z
			rspx = (RX[gi]) - (RX[gj]+SX[gj]*SEGL[gj]);
			rspy = (RY[gi]) - (RY[gj]+SY[gj]*SEGL[gj]);
			rspz = (RZ[gi]) - (RZ[gj]+SZ[gj]*SEGL[gj]);
			temp = sqrt(rspx*rspx + rspy*rspy + rspz*rspz);
			if(temp<=SEGL[gi]*seki_nor_spe) GaussTen = GaussTenSpe;
		}
		else if(GaussTen == GaussTenNor){
			//--�I�_ vs. �n�_�̋����v�Z
			rspx = (RX[gi]+SX[gi]*SEGL[gi]) - (RX[gj]);
			rspy = (RY[gi]+SY[gi]*SEGL[gi]) - (RY[gj]);
			rspz = (RZ[gi]+SZ[gi]*SEGL[gi]) - (RZ[gj]);
			temp = sqrt(rspx*rspx + rspy*rspy + rspz*rspz);
			if(temp<=SEGL[gi]*seki_nor_spe) GaussTen = GaussTenSpe;
		}
		else if(GaussTen == GaussTenNor){
			//--�I�_ vs. �I�_�̋����v�Z
			rspx = (RX[gi]+SX[gi]*SEGL[gi]) - (RX[gj]+SX[gj]*SEGL[gj]);
			rspy = (RY[gi]+SY[gi]*SEGL[gi]) - (RY[gj]+SY[gj]*SEGL[gj]);
			rspz = (RZ[gi]+SZ[gi]*SEGL[gi]) - (RZ[gj]+SZ[gj]*SEGL[gj]);
			temp = sqrt(rspx*rspx + rspy*rspy + rspz*rspz);
			if(temp<=SEGL[gi]*seki_nor_spe) GaussTen = GaussTenSpe;
		}
	}
	//----------------------------------------------------
	
	//----------------------------------------------------
	// ���l�ϕ�
	for(gaussi = 0; gaussi < GaussTen; ++gaussi){
		if(GaussTen == GaussTenNor){
			gausst = GaussBuntenNor[gaussi] * SEGL[gj]/2.0 + SEGL[gj]/2.0;
			//--�W���v�Z
			if(gj == gn-1)		xnj = k0 * gausst;
			else if(gj == gn)	xnj = k0 * (SEGL[gn] - gausst);
			if(hmode == 1)		sin_cos = cos(xnj);
			else if(hmode == 2) sin_cos = sin(xnj);
			//--�����v�Z
			temp_x = (RX[gi]+SX[gi]*gxi) - (RX[gj]+SX[gj]*gausst);
			temp_y = (RY[gi]+SY[gi]*gxi) - (RY[gj]+SY[gj]*gausst);
			temp_z = (RZ[gi]+SZ[gi]*gxi) - (RZ[gj]+SZ[gj]*gausst);
			rowij = sqrt( temp_x*temp_x + temp_y*temp_y + temp_z*temp_z + RA[gi]*RA[gi] );
			//--�ϕ��֐�			
			temp =  ( sin_cos * GaussWeightNor[gaussi] * SEGL[gj]/(2.0 * rowij) ) ;
		}
		else if(GaussTen == GaussTenSpe){
			gausst = GaussBuntenSpe[gaussi] * SEGL[gj]/2.0 + SEGL[gj]/2.0;
			//--�W���v�Z
			if(gj == gn-1)		xnj = k0 * gausst;
			else if(gj == gn)	xnj = k0 * (SEGL[gn] - gausst);
			if(hmode == 1)		sin_cos = cos(xnj);
			else if(hmode == 2) sin_cos = sin(xnj);
			//--�����v�Z
			temp_x = (RX[gi]+SX[gi]*gxi) - (RX[gj]+SX[gj]*gausst);
			temp_y = (RY[gi]+SY[gi]*gxi) - (RY[gj]+SY[gj]*gausst);
			temp_z = (RZ[gi]+SZ[gi]*gxi) - (RZ[gj]+SZ[gj]*gausst);
			rowij = sqrt( temp_x*temp_x + temp_y*temp_y + temp_z*temp_z + RA[gi]*RA[gi] );
			//--�ϕ��֐�
			temp =  ( sin_cos * GaussWeightSpe[gaussi] * SEGL[gj]/(2.0 * rowij) ) ;
		}
		//--�W���v�Z
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
	// �_���ł͂Q�d�ϕ��̊O�ɂ��鏈��( k0/sin(d*k0) )
	// �g�����̃Z�O�����g��SEGL[gi]���g�p���邽��
	temp = k0/sin(k0*SEGL[gj]);
	answer.re = answer.re * temp;
	answer.im = answer.im * temp;
	//----------------------------------------------------
	return (answer);	
}


//////////////////////////////////////////////////////////////////////////////////
//				�d�����z�v�Z
//////////////////////////////////////////////////////////////////////////////////
void MakeCurrent(void)				//�d�����z
{
	int i;							//�J�E���^
	double fedv1_abs;				//��P���d�d���̐�Βl
	complex hansya;					//���ˌW��
	//----------------------------------------------------
	// ��ʏo��
	printf("Current Distribution");
	printf("                      ");
	//----------------------------------------------------

	//----------------------------------------------------
	// ���d�d�����	
	//	  ���d������
	//         �@��P���d�̓d���̓C���[�W�@��p����ꍇ   -2�~(1+j0)[V]
	//                          �C���[�W�@��p���Ȃ��ꍇ -1�~(1+j0)[V]
	//         �A��P���d���C���[�W�@�̏ꍇ�C���̋��d�_���C���[�W�@��K�p���Ă���D
	//           (�U���C�ʑ��ɂ�����炸-2�{)
	//         �B���d�_���𕡐��ɂ��āC�d���l���䓙���s���ꍇ�C��P���d�͏����@
	//           �̂܂܂ŁC��P���d�ȊO�̓d���Ő�����s���D�i�����@�͐�΁j
	//    �d���s��ɓd���l�����Ă��邪�C�s�񎮌v�Z���Im[ ]�ɓ����Ă���l�͓d��
	for(i = 0 ;i < NFED ;++i){
		Im[ FEDP[i] ] = FEDV[i];
	}
	//----------------------------------------------------
	
	//----------------------------------------------------
	// �s�񎮌v�Z
	Gauss();
	//----------------------------------------------------
	
	//----------------------------------------------------
	// ���̓C���s�[�_���X�v�Z
	//    fedv1_abs�ɂ���
	//	  	  �C���[�W�@��p���Ă���ꍇ[fedv1_abs=2.0]���ƂȂ�C
	//		  �C���[�W�@��p���Ă��Ȃ��ꍇ[fedv1_abs=1.0]�Ƃ���D
	//	  �C���s�[�_���X���v�Z���鎞�C�uNFED�̑S�āv�̋��d�d����[fedv1_abs]
	//    �Ŋ��邱�ƂŁC�C���[�W�@�ɂ��d���Q�{�ɑΉ�����D
	fedv1_abs = Abs(FEDV[0]);
	for(i = 0; i < NFED ;++i){
		Zin[i] = (-1.0*FEDV[i]/fedv1_abs) / Im[FEDP[i]];		//���d�ʑ��C��
	}
	//----------------------------------------------------

	//----------------------------------------------------
	//�@�e���d�_��50���ɑ΂���VSWR
	for(i = 0; i < NFED ;++i){
		hansya = (Zin[i]-50.0)/(Zin[i]+50.0);
		VSWR_50[i] = (1.0+Abs(hansya))/(1.0-Abs(hansya));
	}
	//----------------------------------------------------

	//----------------------------------------------------
	//�@�e���d�_��75���ɑ΂���VSWR
	for(i = 0; i < NFED ;++i){
		hansya = (Zin[i]-75.0)/(Zin[i]+75.0);
		VSWR_75[i] = (1.0+Abs(hansya))/(1.0-Abs(hansya));
	}
	//----------------------------------------------------
}

//////////////////////////////////////////////////////////////////////////////////
//				�d���ʑ��v�Z
//////////////////////////////////////////////////////////////////////////////////
void MakePhase(void)
{
	int i,w,s,ns;					//�J�E���^
	int flag;						//�ʑ��v�Z�p�t���O
	double temp1,temp2;				//�ʑ��v�Z�p

	ns = 0;
	for(w = 0; w < NWIR; ++w){				//���C���[No.
		for(s = 0; s <SEGN[w]-1; ++s){		//�Z�O�����gNo.
			PhaseIm[ns] =180.0/PI * atan2(Im[ns].im,Im[ns].re);		//�ʑ��v�Z
			if(s != 0){
				flag = 0;					//�Ƃ肠����������
				//----------------------------------------------------
				// �ׂ荇���Z�O�����g(���C���[�̐؂�ڈȊO)�ňʑ�
				// ��360�x�ȏ�ς��ꍇ�̕␳
				// ���̏������s���ƁC�i�s�g�̈ʑ����i�s�g�̂悤�Ɍ�����
				// ���}360�x���s�������Ȃ̂ŁC�{���I�ȕω��͖���
				for(i = 1 ; flag <= 0 ; ++i){		//�uflag==0�v�̊Ԃ͉��x�ł�
					flag = 1;						//�Ƃ肠�����I���t���O
					temp1 = sqrt(pow(PhaseIm[ns]-PhaseIm[ns-1],2));
					temp2 = sqrt(pow(PhaseIm[ns]+360-PhaseIm[ns-1],2));
					if( temp1 > temp2 ){
						PhaseIm[ns] = PhaseIm[ns] + 360.0;
						flag = 0;					//������x
					}
					temp2 = sqrt(pow(PhaseIm[ns]-360-PhaseIm[ns-1],2));
					if( temp1 > temp2 ){
						PhaseIm[ns] = PhaseIm[ns] - 360.0;
						flag = 0;					//������x
					}
				}
				//----------------------------------------------------
			}
			ns = ns + 1 ;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////
//				�K�E�X�̏����@
//////////////////////////////////////////////////////////////////////////////////
void Gauss(void)										//�K�E�X�̏����@
{
	int m,n,j,k;										//�J�E���^
	//----------------------------------------------------
	// �O�i����
	complex temp,u;
	for(m = 0 ;m < NSEG0 ;++m){
		//----------------------------------------------------
		// ��ʏo��
		printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b");			//�O��̉�ʏo�͂�����
		printf("%6d /%6d",m+1,NSEG0);					//�v�Z���Ă���s����ʏo��
		//----------------------------------------------------
		if( (Zmn[m][m].re == 0.0) && (Zmn[m][m].im == 0.0) ){
			Pibot( m );									//�s�{�b�g�̊m�F
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
			temp = Zmn[k][m] ;							//�|�o�����߂̌W���̎擾
			for(j=m+1 ; j < NSEG0 ; ++j){				//�Ίp�����ȉ��̑|�o��
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
	// ��ޑ��
	for(m = 0 ;m < NSEG0 ;++m){
		temp = Complex(0.0,0.0);
		for(int k = NSEG0-m ;k < NSEG0 ;++k){			//m��0(1���)�̎��ɂ͓���Ȃ�
			temp = temp + Zmn[NSEG0-1-m][k]*Im[k];
		}
		Im[NSEG0-1-m] = Im[NSEG0-1-m] - temp;
	}
	//----------------------------------------------------
}

void Pibot(int i)										//�s�{�b�g�̑I��(Zmn[i][i]��0�̔���)
{
	//----------------------------------------------------
	// ���l��0��Ŋ���Ȃ��悤�ɂ���
	// �A���e�i�̌`����͓�������Ȃ�΁C�K�v�Ȃ�����
	// �K�E�X�̏����@�ŁC���ꂪ0�ɂȂ�Ȃ��悤�ɂ��Ă���
	complex temp;
	if( Abs(Zmn[i][i]) == 0.0 ){						//Zmn[i-1][i-1]=0�̎�����
		for(int x=i+1 ;Abs(Zmn[i][i]) == 0.0 ;++x){
			//----------------------------------------------------
			// x���Ō�̍s�ɂȂ�Ɠ���
			if(x>=NSEG0){								//�v�Z���s�s�\
				OutputConf();							//�`��o��
				printf("\n\nDivision by Zero !!\n\n");	//����=0�ɂȂ�ꍇ(�G���[)
				getch();
				exit ( 0 );								//�I��
			}
			//----------------------------------------------------

			//----------------------------------------------------
			// �s�̓���ւ�
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
//				�Z�O�����g��̔C�ӈʒu�̓d���v�Z
//////////////////////////////////////////////////////////////////////////////////
complex Idx(int wn,int sn, double dx)		//�d���v�Z (���C���[No 
											//		   ,���C���[���̃Z�O�����gNo
											//         ,�ʒu)
{
	//----------------------------------------------------
	//   ���ˊE���v�Z����Ƃ��C�Z�O�����g�̒[�̓d�������ŉ����E���v�Z����ƁC�덷�̌�
	//   �v���O�������ł́C�Z�O�����g��̂ǂ̈ʒu�ł��C�d���l���Z�o�ł���悤�ɂ��Ă���D
	//----------------------------------------------------
	complex answer1 = Complex(0.0 , 0.0);	//���ʑ���p
	complex answer2 = Complex(0.0 , 0.0);	//���ʑ���p
	int INo;								//�Z�O�����g�n�_�̓d��No (�ŏ��� "0" )
	int SNo;								//�v�Z�Ώۂ̃Z�O�����gNo
	int i;									//�J�E���^

	//----------------------------------------------------
	// �d��No , �Z�O�����gNo�@�Z�o
	SNo = 0;
	for(i = 0 ;i < wn; ++i){
		SNo = SNo + SEGN[i] ;
	}
	SNo = SNo + sn;
	INo = SNo - wn;
	//---------------------------------------------------

	//----------------------------------------------------
	// ���C���[�̎n�_���̃Z�O�����g�p�d��
	if(sn != SEGN[wn]-1){							//�I�_�̎��͌v�Z���Ȃ�
		answer1 = Im[INo] * sin(k0*(dx)) / sin(k0*(SEGL[SNo]));
	}
	//---------------------------------------------------

	//----------------------------------------------------
	// ���C���[�̏I�_���̃Z�O�����g�p�d��
	if(sn != 0){									//�n�_�̎��͌v�Z���Ȃ�
		answer2 = Im[INo-1] * sin(k0*(SEGL[SNo]-dx)) / sin(k0*(SEGL[SNo]));
	}
	//---------------------------------------------------

	//---------------------------------------------------
	// �ʒu���ɓd���l��Ԃ�
	if(sn == 0)					return (answer1);	//���C���[�̎n�_�̃Z�O�����g�p�d��
	else if(sn == SEGN[wn]-1)	return (answer2);	//���C���[�̏I�_�̃Z�O�����g�p�d��
	else 	return (answer1 + answer2);				//���̑��̈ʒu�̃Z�O�����g�p�d��
	//---------------------------------------------------
}

//////////////////////////////////////////////////////////////////////////////////
//				���ˊE�v�Z
//////////////////////////////////////////////////////////////////////////////////
void Radiation(void)
{
	double temp;
	double i;		//�v�Z�p�x
	int n;			//�J�E���^
	complex Zi0;	//��1���d�̓��̓C���s�[�_���X
	complex Ii0;	//��1���d���̓d��

	//----------------------------------------------------
	// �����v�Z�p�W��
	//   ��1���d�̂ݗ����v�Z�Ɏg�p
	//     ���������d�͗����␳�K�v
	//       �����ŗ����␳���s�������o���邪�C����������Ɣėp���������Ȃ��\�����L��
	//       �̂ŁC���̃v���O�����ł́C�����v�Z���ɑ�P���d�_�̓d�͂̂ݗp���C��ɗ�����
	//		 �����s��
	Ii0 = Im[ FEDP[0] ];				//��P���d�_�̓d��
	Zi0 = Zin[0];						//��P���d�_�̓��̓C���s�[�_���X
	//----------------------------------------------------

	i = 0.0;
	n = 0;
	
	//----------------------------------------------------
	// ��ʏo��
	printf("Radiation Pattern               %10.2f",i);
	//----------------------------------------------------

	for(i = DegStart ;i <= DegStart+DegWidth ;i = i + DegDelta){
		//----------------------------------------------------
		// ��ʕ\��
		printf("\b\b\b\b\b\b\b\b\b\b");
		printf("%10.2f",i);
		//----------------------------------------------------

		//----------------------------------------------------
		// �����Δg�̉����E
		if(AXMODE == 0){				//�ӌŒ�̏ꍇ
			TRAD[n] = FarFieldT(i , FixAngle);
			FRAD[n] = FarFieldF(i , FixAngle);
		}
		else if(AXMODE == 1){			//�ƌŒ�̏ꍇ
			TRAD[n] = FarFieldT(FixAngle , i);
			FRAD[n] = FarFieldF(FixAngle , i);
		}
		FPHASE[n] =180.0/PI * atan2(FRAD[n].im,FRAD[n].re);		//�ʑ�
		TPHASE[n] =180.0/PI * atan2(TRAD[n].im,TRAD[n].re);		//�ʑ�
		//----------------------------------------------------

		//----------------------------------------------------
		// �~�Δg�̉����E
		LRAD[n] = 0.5 * (TRAD[n]-FRAD[n]*J) ;					//�����~�Δg
		RRAD[n] = 0.5 * (TRAD[n]+FRAD[n]*J) ;					//�E���~�Δg
		RPHASE[n] =180.0/PI * atan2(RRAD[n].im,RRAD[n].re);		//�ʑ�
		LPHASE[n] =180.0/PI * atan2(LRAD[n].im,LRAD[n].re);		//�ʑ�
		//----------------------------------------------------

		//----------------------------------------------------
		// ����̌v�Z
		temp = (Abs(RRAD[n])+Abs(LRAD[n]))/ (Abs(RRAD[n])-Abs(LRAD[n]));
		temp = sqrt(temp*temp);									//��Βl(�t����ɑΉ�)
		if(temp >= 1000) temp = 1000;							//����̍ő�l��60dB
		AR[n] = 20*log10(temp);		
		//----------------------------------------------------

		//----------------------------------------------------
		// �����Δg�Ƃ̗����v�Z
		temp = pow(Abs(TRAD[n]),2)*pow(R,2) / (30.0*pow(Abs(Ii0),2)*Zi0.re) ;
		if(temp <= 0.000001) temp = 0.000001;					//�����ŏ��l��-60dB
		TGAI[n] = 10*log10( temp );
		//----------------------------------------------------
		
		//----------------------------------------------------
		// �����Δg�ӂ̗����v�Z
		temp = pow(Abs(FRAD[n]),2)*pow(R,2) / (30.0*pow(Abs(Ii0),2)*Zi0.re) ;
		if(temp <= 0.000001) temp = 0.000001;					//�����ŏ��l��-60dB
		FGAI[n] = 10*log10( temp );
		//----------------------------------------------------
		
		//----------------------------------------------------
		// �����Δg�Ƃƃӂ̗����v�Z
		temp = (pow(Abs(TRAD[n]),2)+pow(Abs(FRAD[n]),2))*pow(R,2)
			 / (30.0*pow(Abs(Ii0),2)*Zi0.re) ;
		if(temp <= 0.000001) temp = 0.000001;					//�����ŏ��l��-60dB
		TFGAI[n] = 10*log10( temp );
		//----------------------------------------------------

		//----------------------------------------------------
		// �~�ΔgR�̗����v�Z
		temp = 2.0*pow(Abs(RRAD[n]),2)*pow(R,2) / (30.0*pow(Abs(Ii0),2)*Zi0.re) ;
		if(temp <= 0.000001) temp = 0.000001;					//�����ŏ��l��-60dB
		RGAI[n] = 10*log10( temp );
		//----------------------------------------------------
		
		//----------------------------------------------------
		// �~�ΔgL�̗����v�Z
		temp = 2.0*pow(Abs(LRAD[n]),2)*pow(R,2) / (30.0*pow(Abs(Ii0),2)*Zi0.re) ;
		if(temp <= 0.000001) temp = 0.000001;					//�����ŏ��l��-60dB
		LGAI[n] = 10*log10( temp );
		//----------------------------------------------------
		
		//----------------------------------------------------
		// �~�ΔgR��L�̗����v�Z
		temp = (2.0*pow(Abs(RRAD[n]),2)+2.0*pow(Abs(LRAD[n]),2))*pow(R,2)
			 / (30.0*pow(Abs(Ii0),2)*Zi0.re) ;
		if(temp <= 0.000001) temp = 0.000001;					//�����ŏ��l��-60dB
		RLGAI[n] = 10*log10( temp );
		//----------------------------------------------------
		n = n + 1;
	}
	//----------------------------------------------------
	// �����Δg�̕��˃p�^�[���v�Z
	n = 0;
	temp = TGAI[0];
	for(i = DegStart ;i <= DegStart+DegWidth ;i = i + DegDelta){
		if(temp < TGAI[n])	temp = TGAI[n];		//�ő嗘���l�T��
		if(temp < FGAI[n])	temp = FGAI[n];		//�ő嗘���l�T��
		n = n + 1;	
	}
	n = 0;
	for(i = DegStart ;i <= DegStart+DegWidth ;i = i + DegDelta){
		TPAT[n] = TGAI[n] - temp;				//[����]-[�ő嗘��]
		FPAT[n] = FGAI[n] - temp;				//[����]-[�ő嗘��]
		if(TPAT[n] < -30.0) TPAT[n] = -30.0;	//�p�^�[���ŏ��l��-30dB
		if(FPAT[n] < -30.0) FPAT[n] = -30.0;	//�p�^�[���ŏ��l��-30dB
		n = n + 1;	
	}
	//----------------------------------------------------

	//----------------------------------------------------
	// �~�Δg�̕��˃p�^�[���v�Z
	n = 0;
	temp = RGAI[0];
	for(i = DegStart ;i <= DegStart+DegWidth ;i = i + DegDelta){
		if(temp < RGAI[n])	temp = RGAI[n];		//�ő嗘���l�T��
		if(temp < LGAI[n])	temp = LGAI[n];		//�ő嗘���l�T��
		n = n + 1;	
	}
	n = 0;
	for(i = DegStart ;i <= DegStart+DegWidth ;i = i + DegDelta){
		RPAT[n] = RGAI[n] - temp;				//[����]-[�ő嗘��]
		LPAT[n] = LGAI[n] - temp;				//[����]-[�ő嗘��]
		if(RPAT[n] < -30.0) RPAT[n] = -30.0;	//�p�^�[���ŏ��l��-30dB
		if(LPAT[n] < -30.0) LPAT[n] = -30.0;	//�p�^�[���ŏ��l��-30dB
		n = n + 1;	
	}
	//----------------------------------------------------
	printf("\n");
}

complex FarFieldT(double theta,double fai)		//�����E�Ɛ��� (�ƁC��)����
{
	int w,s;				//���C���[�@�Z�O�����g�̃J�E���^
	int ns;					//�Z�O�����g�i���o�[
	complex Ith;			//I��
	double delR;			//��R
	double radth , radfa;	//�Ƃƃӂ̃��W�A���l
	int gaussi;				//�ϕ��Ŏg���J�E���^
	double gausst;			//�ꎟ�ϊ����ʂ̑��
	complex answer;			//���ʑ���p
	complex answer_s;		//�e�Z�O�����g�̌��ʑ���p

	//----------------------------------------------------
	// deg -> rad �̕ϊ�
	radth = theta * PI / 180.0;
	radfa =   fai * PI / 180.0;
	//----------------------------------------------------

	answer = Complex(0.0 , 0.0);	//���ʑ���p
	ns = 0;
	for(w = 0; w < NWIR; ++w){
		for(s = 0; s < SEGN[w]; ++s){
			answer_s = Complex(0.0 , 0.0);
			for(gaussi = 0; gaussi < GaussTenNor; ++gaussi){
				gausst = GaussBuntenNor[gaussi] * SEGL[ns]/2.0 + SEGL[ns]/2.0;
				//--I�� (�d���̃ƕ�������)
				Ith = Idx(w,s,gausst)
						*(  SX[ns] * cos(radth) * cos(radfa)
						  + SY[ns] * cos(radth) * sin(radfa)
						  - SZ[ns] * sin(radth)   );
				//--��R
				delR = (RX[ns]+SX[ns]*gausst) * sin(radth) * cos(radfa)
				     + (RY[ns]+SY[ns]*gausst) * sin(radth) * sin(radfa)
					 + (RZ[ns]+SZ[ns]*gausst) * cos(radth);
				//--��ϕ��֐�
				answer_s = answer_s
					     + Ith * Exp(J*(Beta*delR)) 
						       * (GaussWeightNor[gaussi] * SEGL[ns]/2.0);
			}
			
			//----------------------------------------------------
			// ���˃��[�h��K�p
			//    RAD_MODE[ ]���u0�v�̏ꍇ�C�����ŕ��˂��~�߂�
			//    RAD_MODE[ ]��double�^�Ő錾����΁C�Z�O�����g���ɕ��˗����w�肷�鎖���\
			//    �A���e�i�̉�͂Ŏg��Ȃ��̂�RAD_MODE[ ]��int�^�ɂ��Ă���(0���܂���100��)
			answer_s = answer_s * RAD_MODE[ns];
			//----------------------------------------------------
			answer = answer + answer_s;
			ns = ns + 1 ;
		}
	}
	answer = answer * ((-1.0)*2.0*PI*USEF*u0/(4.0*PI*R))*(Exp((-1.0*Beta*R)*J)*J);
	return(answer);
}

complex FarFieldF(double theta,double fai)				//�����E�Ӑ��� (�ƁC��)����
{
	int w,s;				//���C���[�@�Z�O�����g�̃J�E���^
	int ns;					//�Z�O�����g�i���o�[
	complex Ifa;			//I��
	double delR;			//��R
	double radth , radfa;	//���W�A��

	int gaussi;				//�ϕ��Ŏg���J�E���^
	double gausst;			//�ꎟ�ϊ����ʂ̑��
	complex answer;			//���ʑ���p
	complex answer_s;		//�e�Z�O�����g�̌��ʑ���p

	//----------------------------------------------------
	// deg -> rad �̕ϊ�
	radth = theta * PI / 180.0;
	radfa =   fai * PI / 180.0;
	//----------------------------------------------------

	answer = Complex(0.0 , 0.0);	//���ʑ���p
	ns = 0;
	for(w = 0; w < NWIR; ++w){
		for(s = 0; s <SEGN[w]; ++s){
			answer_s = Complex(0.0 , 0.0);
			for(gaussi = 0; gaussi < GaussTenNor; ++gaussi){
				gausst = GaussBuntenNor[gaussi] * SEGL[ns]/2.0 + SEGL[ns]/2.0;
				//--I�� (�d���̃ӕ�������)
				Ifa = Idx(w,s,gausst)
						*(  SX[ns] * sin(radfa) * (-1.0)
						  + SY[ns] * cos(radfa)   );
				//--��R
				delR = (RX[ns]+SX[ns]*gausst) * sin(radth) * cos(radfa)
				     + (RY[ns]+SY[ns]*gausst) * sin(radth) * sin(radfa)
					 + (RZ[ns]+SZ[ns]*gausst) * cos(radth);
				//--��ϕ��֐�
				answer_s = answer_s
					     + Ifa * Exp(J*(Beta*delR)) 
						       * (GaussWeightNor[gaussi] * SEGL[ns]/2.0);
			}
			//----------------------------------------------------
			// ���˃��[�h��K�p
			//    RAD_MODE[ ]���u0�v�̏ꍇ�C�����ŕ��˂��~�߂�
			//    RAD_MODE[ ]��double�^�Ő錾����΁C�Z�O�����g���ɕ��˗����w�肷�鎖���\
			//    �A���e�i�̉�͂Ŏg��Ȃ��̂�RAD_MODE[ ]��int�^�ɂ��Ă���(0���܂���100��)
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
//				�t�@�C���J��
//////////////////////////////////////////////////////////////////////////////////
void FileOpen(void)
{
	_mkdir("AntennaData") ;											//�f�[�^�̏o�̓t�@�C��
	//----------------------------------------------------
	// �`��p�t�@�C�����J�� �����������֎~
	//    �t�@�C�����⒆�g������������ƁC�`�������v���O
	//    ����[mm_conf.exe]���g���Ȃ��Ȃ�̂ŁC���������֎~
	if((fp_conf=fopen("AntennaData\\antenna_conf.dat","w"))==NULL){		//�G���[�̏ꍇ�NULL�
		printf("write open error [ antenna_conf.dat ] !!\n\n");			//�I�[�v���G���[
		getch();
		exit(0);														//�G���[�̏ꍇ�~�߂�
	}
	fprintf(fp_conf," [antenna_conf.dat]\n\n");
	//----------------------------------------------------

	//----------------------------------------------------
	// �d���p�t�@�C�����J��
	if((fp_curr=fopen("AntennaData\\antenna_current.dat","w"))==NULL){	//�G���[�̏ꍇ�NULL�
		printf("write open error [ antenna_current.dat ] !!\n\n");		//�I�[�v���G���[
		getch();
		exit(0);														//�G���[�̏ꍇ�~�߂�
	}
	fprintf(fp_curr," [antenna_current.dat]\n\n");
	//----------------------------------------------------

	//----------------------------------------------------
	// ���ˊE�p�t�@�C�����J��
	if((fp_radf=fopen("AntennaData\\antenna_field.dat","w"))==NULL){	//�G���[�̏ꍇ�NULL�
		printf("write open error [ antenna_field.dat ] !!\n\n");		//�I�[�v���G���[
		getch();
		exit(0);														//�G���[�̏ꍇ�~�߂�
	}
	fprintf(fp_radf," [antenna_field.dat]\n\n");
	//----------------------------------------------------

	//----------------------------------------------------
	// ���˃p�^�[���p�t�@�C�����J��
	if((fp_radp=fopen("AntennaData\\antenna_pattern.dat","w"))==NULL){	//�G���[�̏ꍇ�NULL�
		printf("write open error [ antenna_pattern.dat ] !!\n\n");		//�I�[�v���G���[
		getch();
		exit(0);														//�G���[�̏ꍇ�~�߂�
	}
	fprintf(fp_radp," [antenna_pattern.dat]\n\n");
	//----------------------------------------------------

	//----------------------------------------------------
	// �����f�[�^�o��1�p�t�@�C�����J��
	if((fp_char1=fopen("AntennaData\\antenna_chara1.dat","w"))==NULL){	//�G���[�̏ꍇ�NULL�
		printf("write open error [ antenna_chara1.dat ] !!\n\n");		//�I�[�v���G���[
		getch();
		exit(0);														//�G���[�̏ꍇ�~�߂�
	}
	fprintf(fp_char1," [antenna_chara1.dat]\n\n");
	fprintf(fp_char1," >>---------PARA ---------<< **");
	fprintf(fp_char1," >>------- Impedance -------- ------<< **");
	fprintf(fp_char1," >>------ Max-Direction Gain(dB) ----<< **");
	fprintf(fp_char1," >>---- 0(deg)-Gain(dB) ----<<\n");

	fprintf(fp_char1,"         para1         para2         para3 **");
	fprintf(fp_char1,"    Zin.re    Zin.im  VSWR-50  VSWR-75 **");
	fprintf(fp_char1,"    Angle       E��       E��   E��+E�� **");
	fprintf(fp_char1,"       E��       E��   E��+E��\n");
	//----------------------------------------------------

	//----------------------------------------------------
	// �����f�[�^�o��2�p�t�@�C�����J��
	if((fp_char2=fopen("AntennaData\\antenna_chara2.dat","w"))==NULL){	//�G���[�̏ꍇ�NULL�
		printf("write open error [ antenna_chara2.dat ] !!\n\n");		//�I�[�v���G���[
		getch();
		exit(0);														//�G���[�̏ꍇ�~�߂�
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
	// �C�Ӄf�[�^�o�͗p�t�@�C�����J��
	//   �����o�͂��Ă����Ȃ����C�o�͍�Ƃ�"antenna.h"�ōs��
	if((fp_free=fopen("AntennaData\\antenna_free.dat","w"))==NULL){		//�G���[�̏ꍇ�NULL�
		printf("write open error [ antenna_free.dat ] !!\n\n");			//�I�[�v���G���[

		getch();
		exit(0);														//�G���[�̏ꍇ�~�߂�
	}
	fprintf(fp_free," [antenna_free.dat]\n\n");
	//----------------------------------------------------
}

//////////////////////////////////////////////////////////////////////////////////
//				�t�@�C������
//////////////////////////////////////////////////////////////////////////////////
void FileClose(void)
{
	fclose(fp_curr);		//�d�����z
	fclose(fp_conf);		//�`��
	fclose(fp_radf);		//���ˊE
	fclose(fp_radp);		//�p�^�[��
	fclose(fp_char1);		//�����P
	fclose(fp_char2);		//�����Q
	fclose(fp_free);		//�K���o��
}

//////////////////////////////////////////////////////////////////////////////////
//			�����f�[�^�o��	antenna_chara1.dat  antenna_chara2.dat
//////////////////////////////////////////////////////////////////////////////////
void OutputChara(void)							//�����f�[�^�o��
{
	double gaimax;		//�ő嗘������p
	double i;			//�v�Z�p�x
	int n;				//�z��No
	int maxang,zang;	//�ő����No�CZ������No

	//----------------------------------------------------
	//	�����Δg���� antenna_chara1
	//--�ő���˕����T��

	n = 0;//������
	maxang = 0;//������
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

	//--Z�������T��(0�x�����T��)�ƌŒ�̏ꍇX�������T��
	n = 0;
	for(i = DegStart ;i <= DegStart+DegWidth ;i = i + DegDelta){
		if( (int)(i*10000000)==0.0){
			zang = n;
		}
		n = n + 1;
	}

	//--�����o��
	//    PARA1,PARA2
	//      �����o�͂ŏo�͂��Ă���PARA1,PARA2�͎��g�������Ȃǂ��v�Z�����ꍇ�ɕ\
	//      ��O���t�����₷�����邽�߂̃��m�v�Z���ł͉����g�p���Ă��Ȃ��̂ŁC
	//      �o�͂��Ȃ��Ă����Ȃ��ꍇ�͏����Ȃ�
	//    �\������
	//      �A���e�i�̌������s����ŁC�K�v�ȏ�̏����ȉ��̌������o�͂��Ă��Ȃ��D
	//      ����ȏ�̌������K�v�Ȃ�΁C�C�ӏo�͂�p����
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
	//	�~�Δg���� antenna_chara2
	//--�ő���˕����T��
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

	//--Z�������T��(0�x�����T��)�ƌŒ�̏ꍇX�������T��
	n = 0;
	for(i = DegStart ;i <= DegStart+DegWidth ;i = i + DegDelta){
		if( (int)(i*10000000)==0.0){
			zang = n;
		}
		n = n + 1;
	}

	//--�����o��
	//    PARA1,PARA2
	//      �����o�͂ŏo�͂��Ă���PARA1,PARA2�͎��g�������Ȃǂ��v�Z�����ꍇ�ɕ\
	//      ��O���t�����₷�����邽�߂̃��m�v�Z���ł͉����g�p���Ă��Ȃ��̂ŁC
	//      �o�͂��Ȃ��Ă����Ȃ��ꍇ�͏����Ȃ�
	//    �\������
	//      �A���e�i�̌������s����ŁC�K�v�ȏ�̏����ȉ��̌������o�͂��Ă��Ȃ��D
	//      ����ȏ�̌������K�v�Ȃ�΁C�C�ӏo�͂�p����
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
//			�`��o��	antenna_conf.dat
//			********  ���������s��  ********
//				�`�������v���O����[mm_conf.exe]���g��Ȃ��ꍇ�͏���������
//////////////////////////////////////////////////////////////////////////////////
void OutputConf(void)
{
	int i;			//�J�E���^
	//----------------------------------------------------
	// �������
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
	// �Z�O�����g���̃f�[�^
	//     �ʒu�Ⓑ����g���P�ʂŏo�͂���
	//       �A���e�i�̐݌v���g���P�ʂ�����
	//       �`�������Ƃ��ɁC���g���Ɉˑ����Ȃ��Ō��₷��
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
	// ���C���[���̃f�[�^
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
	// ���d���̃f�[�^
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
	// ��R���̃f�[�^
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
//			�d�����zI�̏o��		antenna_current.dat
//////////////////////////////////////////////////////////////////////////////////
void OutputCurrent(void)
{
	int i;			//�J�E���^
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
//			���ˊE�@���˃p�^�[���o��	antenna_field.dat   antenna_pattern.dat
//////////////////////////////////////////////////////////////////////////////////
void OutputRad(void)
{
	double i,delp;			//�p�x�C�ʑ���
	int n;					//�J�E���^

	//========================================================
	//			���ˊE		antenna_field.dat
	//========================================================

	//----------------------------------------------------
	// �����Δg�֌W
	if(AXMODE == 0)fprintf(fp_radf,"  �� = %.3f ��plane\n",FixAngle);
	if(AXMODE == 1)fprintf(fp_radf,"  �� = %.3f ��plane\n",FixAngle);
	fprintf(fp_radf,"   --    ");
	fprintf(fp_radf,">>-------   -------   --E��--   -------   -------<<");
	fprintf(fp_radf," ** ");
	fprintf(fp_radf,">>-------   -------   --E��--   -------   -------<<");
	fprintf(fp_radf," ** ");
	fprintf(fp_radf,">>-------   --------<<\n");
	if(AXMODE == 0)fprintf(fp_radf," ��(deg)");
	if(AXMODE == 1)fprintf(fp_radf," ��(deg)");
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
	// �~�Δg�֌W
	fprintf(fp_radf,"\n\n   --    ");
	fprintf(fp_radf,">>-------   -------   --ER --   -------   -------<<");
	fprintf(fp_radf," ** ");
	fprintf(fp_radf,">>-------   -------   --EL --   -------   -------<<");
	fprintf(fp_radf," ** ");
	fprintf(fp_radf,">>-------   --------<<\n");
	if(AXMODE == 0)fprintf(fp_radf," ��(deg)");
	if(AXMODE == 1)fprintf(fp_radf," ��(deg)");
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
	//			���˃p�^�[��	antenna_pattern.dat
	//========================================================
	if(AXMODE == 0)fprintf(fp_radp,"  �� = %.3f ��plane\n",FixAngle);
	if(AXMODE == 1)fprintf(fp_radp,"  �� = %.3f ��plane\n",FixAngle);
	if(AXMODE == 0)fprintf(fp_radp,"   ��(deg)");
	if(AXMODE == 1)fprintf(fp_radp,"   ��(deg)");
	fprintf(fp_radp,"       E��       E��       ER        EL \n");
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
//				�v�Z���ԑ���
//////////////////////////////////////////////////////////////////////////////////
void StartTime(void)						//�^�C�}�[�X�^�[�g
{
	time( &start );							//�^�C�}�[�X�^�[�g
}

void   LapsedTime(void)                     //�o�ߎ���
{
	double d_time;
	time( &finish );						//���Ԏ擾
	d_time = difftime( finish, start );		//���ԍ��v�Z
	la_sec=(int)(d_time);					//�b�P�ʂő��
	la_min=(int)(la_sec/60);				//���P�ʂő��
	la_sec=la_sec%60;						//[�b�P�ʂ̎��ԁ�60]�̗]��=�b
	printf("Lapsed Time                         %3d:%2d\n",la_min,la_sec);
}

//////////////////////////////////////////////////////////////////////////////////
//				�z��̊m��
//					"makemat.h"�̒��ɂ���֐����Ăяo���z��̊m�ۂ��s���D
//////////////////////////////////////////////////////////////////////////////////
void MakeMatAll(void)						//�z��m��
{
	int NN;									//���ˊE�̌v�Z�_������p

	//--�`��z��
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
	LOADP = MakeMat1D_int( NLOAD );			//�e�C���s�[�_���X���׈ʒu
	LOADZ = MakeMat1D_complex( NLOAD );		//�e�C���s�[�_���X���גl

	//--�d���z��
	Zmn = MakeMat2D_complex( NSEG0 );		//complex Zmn[NSEG0][NSEG0]
	Im = MakeMat1D_complex( NSEG0 );		//complex Im[NSEG0]
	PhaseIm = MakeMat1D_double( NSEG0 );	//double PhaseIm[NSEG0]
	Zin = MakeMat1D_complex( NFED );		//�e���d�_�̓��̓C���s�[�_���X
	VSWR_50 = MakeMat1D_double( NFED );		//�e���d�_��50���ɑ΂���VSWR
	VSWR_75 = MakeMat1D_double( NFED );		//�e���d�_��75���ɑ΂���VSWR

	//--���ˊE�v�Z�p
	NN = int(DegWidth/DegDelta) + 1;		//���ˊE�̌v�Z�_��(+1�͗\��)
	TRAD = MakeMat1D_complex( NN );			//E�Ƃ̕��ˊE TRAD[  ]
	FRAD = MakeMat1D_complex( NN );			//E�ӂ̕��ˊE FRAD[  ]
	RRAD = MakeMat1D_complex( NN );			//ER �̕��ˊE RRAD[  ]
	LRAD = MakeMat1D_complex( NN );			//EL �̕��ˊE LRAD[  ]
	TPHASE = MakeMat1D_double( NN );		//E�Ƃ̈ʑ� TPHASE[  ]
	FPHASE = MakeMat1D_double( NN );		//E�ӂ̈ʑ� FPHASE[  ]
	RPHASE = MakeMat1D_double( NN );		//ER �̈ʑ� RPHASE[  ]
	LPHASE = MakeMat1D_double( NN );		//EL �̈ʑ� LPHASE[  ]
	TGAI = MakeMat1D_double( NN );			//E�Ƃ̗��� TGAI[  ]
	FGAI = MakeMat1D_double( NN );			//E�ӂ̗��� FGAI[  ]
	TFGAI = MakeMat1D_double( NN );			//E�Ƃ�E�ӂ̍��v�̗��� TFGAI[  ]
	RLGAI = MakeMat1D_double( NN );			//ER ��EL �̍��v�̗��� RLGAI[  ]
	RGAI = MakeMat1D_double( NN );			//ER �̗��� RGAI[  ]
	LGAI = MakeMat1D_double( NN );			//EL �̗��� LGAI[  ]
	TPAT = MakeMat1D_double( NN );			//E�Ƃ̃p�^�[�� TPAT[  ]
	FPAT = MakeMat1D_double( NN );			//E�ӂ̃p�^�[�� FPAT[  ]
	RPAT = MakeMat1D_double( NN );			//ER �̃p�^�[�� RPAT[  ]
	LPAT = MakeMat1D_double( NN );			//EL �̃p�^�[�� LPAT[  ]
	AR = MakeMat1D_double( NN );			//���� AR[ ]

	//--�K�E�X�ϕ��̃p�����[�^
	GaussWeightNor = MakeMat1D_double( GaussTenNor );	//�ʏ�p
	GaussBuntenNor = MakeMat1D_double( GaussTenNor );	//�ʏ�p
	GaussWeightSpe = MakeMat1D_double( GaussTenSpe );	//���ٓ_�p
	GaussBuntenSpe = MakeMat1D_double( GaussTenSpe );	//���ٓ_�p
}

//////////////////////////////////////////////////////////////////////////////////
//				�z��̊J��
//////////////////////////////////////////////////////////////////////////////////
void DelMatAll(void)			//�z��J��
{
	//--�`��z��
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
	DelMat1D_int( LOADP );			//�e�C���s�[�_���X���׈ʒu
	DelMat1D_complex( LOADZ );		//�e�C���s�[�_���X���גl

	//--�d���z��
	DelMat2D_complex(Zmn,NSEG0);	//complex Zmn[NSEG0][NSEG0]
	DelMat1D_complex(Im);			//complex Im[NSEG0]
	DelMat1D_double(PhaseIm);		//double PhaseIm[NSEG0]
	DelMat1D_complex(Zin);			//double Zin[NFED]
	DelMat1D_double(VSWR_50);		//double VSWR_50[NFED]
	DelMat1D_double(VSWR_75);		//double VSWR_75[NFED]

	//--���ˊE�v�Z�p
	DelMat1D_complex(TRAD);			//E�Ƃ̕��ˊE TRAD[  ]
	DelMat1D_complex(FRAD);			//E�ӂ̕��ˊE FRAD[  ]
	DelMat1D_complex(RRAD);			//ER �̕��ˊE RRAD[  ]
	DelMat1D_complex(LRAD);			//EL �̕��ˊE LRAD[  ]
	DelMat1D_double(TPHASE);		//E�Ƃ̈ʑ� TPHASE[  ]
	DelMat1D_double(FPHASE);		//E�ӂ̈ʑ� FPHASE[  ]
	DelMat1D_double(RPHASE);		//ER �̈ʑ� RPHASE[  ]
	DelMat1D_double(LPHASE);		//EL �̈ʑ� LPHASE[  ]
	DelMat1D_double(TGAI);			//E�Ƃ̗��� TGAI[  ]
	DelMat1D_double(FGAI);			//E�ӂ̗��� FGAI[  ]
	DelMat1D_double(TFGAI);			//E�Ƃ�E�ӂ̍��v�̗��� TFGAI[  ]
	DelMat1D_double(RLGAI);			//ER ��EL �̍��v�̗��� RLGAI[  ]
	DelMat1D_double(RGAI);			//ER �̗��� RGAI[  ]
	DelMat1D_double(LGAI);			//EL �̗��� LGAI[  ]
	DelMat1D_double(TPAT);			//E�Ƃ̃p�^�[�� TPAT[  ]
	DelMat1D_double(FPAT);			//E�ӂ̃p�^�[�� FPAT[  ]
	DelMat1D_double(RPAT);			//ER �̃p�^�[�� RPAT[  ]
	DelMat1D_double(LPAT);			//EL �̃p�^�[�� LPAT[  ]
	DelMat1D_double(AR);			//���� AR[ ]

	//--�K�E�X�ϕ��̃p�����[�^
	DelMat1D_double(GaussWeightNor);
	DelMat1D_double(GaussBuntenNor);
	DelMat1D_double(GaussWeightSpe);
	DelMat1D_double(GaussBuntenSpe);
}

//////////////////////////////////////////////////////////////////////////////////
//				�K�E�X�ϕ��̃p�����[�^�쐬
//////////////////////////////////////////////////////////////////////////////////
void MakeGaussPara(void)
{
	//----------------------------------------------------
	// �ʏ�p
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
	// ���ٓ_�p
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
//		����
//
//		ver 0_00	2D��3D�p�ɏ�������
//		ver 0_01	���ٓ_��40�_�Őϕ�
//		ver 0_02	�ς̌v�Z��������
//		ver 0_03	sin(dk0)�̌v�Z��2�d�ϕ��̒��Ɉړ�
//		ver 0_04	�����ݒ�֐���ǉ�
//		ver 0_05	�����p�^�[����ǉ�
//		ver 0_06	if��������v�Z���Ԃ�Z������
//					timer��ǉ�
//					pow�֐�����菜���v�Z���Ԃ�Z������
//		ver 0_07	�`��o�͂ɋ��d�_�ƕ��˃��[�h��ǉ�
//		ver 0_08	���f���v�Z��W�J���Ď��ԒZ�k
//		ver 0_09	�o�̓t�@�C������ς���
//					antenna_chara2.dat��ǉ�
//		ver 0_10	���g��������PARA1�ōs��
//		ver 0_11	�t�@�C������momawa�ɕύX
//		ver 0_12	hij()���̎���if���Ő؂�ւ��ɕύX
//		ver 0_13	sin(dk0)�̌v�Z���ϑ��_�Ɣg���̂��ꂼ��ɕ����Ďg�p
//		ver 0_14	�o�̓t�@�C��(char1,2)��ǉ�
//
//		ver 1_00	�ꕔ�錾���w�b�_�[��
//					�v�Z�����֐�Calculation( )��ǉ�
//		ver 1_01	�ϕ����ʌv�Z��菜��
//		ver 1_02	antenna.h�݂̂Ō`�����
//		ver 1_03	������臒l�v�Z������
//		ver 1_04	�`����͕��@��ύX
//		ver 1_05	�`����͂��~�X���Ă��`����o�͂�����
//		ver 1_06	�o�͏������v�Z�񂵂ɑΉ�
//		ver 1_07	���̓C���s�[�_���X�v�Z�̋��d�ʑ��l��
//		ver 1_08	���ˊE�v�Z�̃��C���[�[(�z��O�Q��)�C��
//					�`�󌟍��ɋ��d�֌W�ǉ�
//		ver 1_09	�C�ӏo�͂�Calculation( )�̒��Ɉړ�
//		ver 1_10	���̓C���s�[�_���X�v�Z�C��
//		ver 1_11	�C���s�[�_���X�s�񑀍�(��R���[�h)���C��
//		ver 1_12	�`��o�͂ɒ�R���[�h��ǉ�
//		ver 1_13	�K�E�X�̏����@ ���ʌv�Z��菜��
//		ver 1_14	���̓C���s�[�_���X�v�Z�C��(��Q���d�_�ȍ~�ɃC���[�W�@�K�p)
//					MakeGaussPara()�̌Ăяo���ʒu��ύX
//					Beta�̌v�Z�������l�v�Z�ɒǉ�
//					�����Δg�̈ʑ����̖��ʂȕ����폜
//		ver 1_15	[Division by Zero]�ł��`��o��
//					�R�����g�ǉ�
//					�����v�Z�ƃC���s�[�_���X�̕t���ɍŋߏ����^���������??
//
//============================================================================
//============================================================================





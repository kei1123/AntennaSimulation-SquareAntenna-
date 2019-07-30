//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//
//			配列動的確保ヘッダー　　makemat.h
//
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
//		配列確保									
//////////////////////////////////////////////////////////////////////////////////
double **MakeMat2D_double( int mat_max )		//2次元double型
{
	double **mat = (double **)malloc(mat_max * sizeof(double *));
	if(mat == NULL){
		printf("Make Matrix Error  2D_double  !!\n");
		getch();
		exit (1);
	}
	else{
		for(int i = 0; i < mat_max; i++){
			mat[i] = (double *)malloc(mat_max * sizeof(double));
			if(mat[i] == NULL){// mat[i] は行ベクトルを指すポインタ
				while(--i >= 0) free(mat[i]);
				free(mat);
				printf("Make Matrix Error  2D_double  !!\n");
				getch();
				exit (1);
			}
		}
	}
	return mat;
}


double  *MakeMat1D_double( int mat_max )		//1次元double型
{
	double *mat = (double *)malloc(mat_max * sizeof(double));
	if(mat == NULL){
		printf("Make Matrix Error  1D_double  !!\n");
		getch();
		exit (1);
	}
	return mat;
}

int **MakeMat2D_int( int mat_max  )				//2次元int型
{
	int **mat = (int **)malloc(mat_max * sizeof(int *));
	if(mat == NULL){
		printf("Make Matrix Error  2D_int  !!\n");
		getch();
		exit (1);
	}
    else{
		for(int i = 0; i < mat_max; i++){
			mat[i] = (int *)malloc(mat_max * sizeof(int));
			if(mat[i] == NULL){
				while(--i >= 0) free(mat[i]);
				free(mat);
				printf("Make Matrix Error  2D_int  !!\n");
				getch();
				exit (1);
			}
		}
	}
	return mat;
}

int  *MakeMat1D_int( int mat_max  )				//1次元int型
{
	int *mat = (int *)malloc(mat_max * sizeof(int));
	if(mat == NULL){
		printf("Make Matrix Error  1D_int  !!\n");
		getch();
		exit (1);
	}
	return mat;
}

complex **MakeMat2D_complex( int mat_max )		//2次元complex型
{
	complex **mat = (complex **)malloc(mat_max * sizeof(complex *));
	if(mat == NULL){
		printf("Make Matrix Error  2D_complex  !!\n");
		getch();
		exit (1);
	}
	else{
		for(int i = 0; i < mat_max; i++){
			mat[i] = (complex *)malloc(mat_max * sizeof(complex));
			if(mat[i] == NULL){
				while(--i >= 0) free(mat[i]);
				free(mat);
				printf("Make Matrix Error  2D_complex  !!\n");
				getch();
				exit (1);
			}
		}
	}
	return mat;
}

complex  *MakeMat1D_complex( int mat_max )		//1次元complex型
{
	complex *mat = (complex *)malloc(mat_max * sizeof(complex));
	if(mat == NULL){
		printf("Make Matrix Error  1D_complex!!\n");
		getch();
		exit (1);
	}
	return mat;
}

//////////////////////////////////////////////////////////////////////////////////
//		配列開放
//////////////////////////////////////////////////////////////////////////////////
void DelMat2D_double(double **mat, int mat_max)
{
	for(int i = 0; i < mat_max; i++){
		free(mat[i]);
	}
	free(mat);
}

void DelMat1D_double(double *mat )
{
	free(mat);
}

void DelMat2D_int(int **mat, int mat_max)
{
	for(int i = 0; i < mat_max; i++){
		free(mat[i]);
	}
	free(mat);
}

void DelMat1D_int(int *mat )
{
	free(mat);
}

void DelMat2D_complex(complex **mat, int mat_max)
{
	for(int i = 0; i < mat_max; i++){
		free(mat[i]);
	}
	free(mat);
}

void DelMat1D_complex(complex *mat )
{
	free(mat);
}
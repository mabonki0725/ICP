/*-------------------------------------------------------------------------
//
// File : ICPSample.c
//
// Discription : Sample code to estimate relative motion with 
//               Iterative Closest Point (ICP) algorithm.
//
// Environment : C
//
// Author : M.nakai (Atsushi Sakai)
//
// 
 -------------------------------------------------------------------------*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "mxlib.h"
#include "usrlib.h"

//ICP パラメータ

#define EPS     0.00001 //収束判定値
#define maxIter 100    //最大イタレーション数
//preError=0;//一つ前のイタレーションのerror値
//dError=1000;//エラー値の差分

double FindNearestPoint(double **,int,double **,int,int,int *,int *);
int    SVDMotionEstimation(double **,double **,int,int,int *,double **,double *);//特異値分解による移動量推定
int    svdcomp(double **a, int m, int n, double *w, double **v);
int    readFile(FILE *fp,double ***,int *,int *);
int    readFree(double **,int,int);
void   *mxAlloc(int n,int m);
/**********
  main
***********/
#if 0
#define NUM   100
#define SIZE   10
#define DIM    2
#define SIG    1
#define RO     10
#define SHIFT   3
int main(argc,argv)
int argc;
char *argv[];
{
    int i,j;
    FILE *fw1,*fw2;
    double rnd,ro;
#if 0    
    double data1[NUM][DIM];
    double data2[NUM][DIM];
    double R[DIM][DIM];
    double t[DIM];
#else
    double **data1;
    double **data2,**data2T;
    double **R;
    double t[DIM];
    double **wrk,**wrkT;
    double center[DIM];

    data1=(double **)comMxAlloc(NUM,DIM,sizeof(double));
    data2=(double **)comMxAlloc(NUM,DIM,sizeof(double));
    data2T=(double **)comMxAlloc(DIM,NUM,sizeof(double));

    R    =(double **)comMxAlloc(DIM,DIM,sizeof(double));
    wrk  =(double **)comMxAlloc(NUM,DIM,sizeof(double));
    wrkT =(double **)comMxAlloc(DIM,NUM,sizeof(double));

#endif

    if(argc < 3) {
      fprintf(stderr,"USAGE command outFile1 outFile2\n");
      return(-9);
    }
    if(!(fw1=fopen(argv[1],"w"))) {
      fprintf(stderr,"Cannot Write fileA=[%s]\n",argv[1]);
      return(-2);
    }
    if(!(fw2=fopen(argv[2],"w"))) {
      fprintf(stderr,"Cannot Write fileA=[%s]\n",argv[2]);
      return(-2);
    }

    for(j=0;j<DIM;j++) {
      for(i=0;i<DIM;i++) {
        R[j][i]=0.0;
      }
      t[j]=0.0;
    }

    for(i=0;i<NUM;i++) {
      for(j=0;j<DIM;j++) {
        data1[i][j] = rand()/(double)RAND_MAX * SIZE;
      }
    }
    /**
    for(j=0;j<DIM;j++) {
      R[j][j]=(5+(10.0*j))/180.0*PAI;
      R[j][j]=10.0/180.0*PAI;
    }
    **/
    ro = RO/180.0*PAI;
    R[0][0] =  cos(ro);
    R[0][1] =  sin(ro);
    R[1][0] = -sin(ro);
    R[1][1] =  cos(ro);
    //R[2][0] =  cos(ro);
    //R[2][1] =  sin(ro);
    
    /*
    for(j=0;j<DIM;j++) {
      t[j]= 2.0 + j*2.0;
    }
    */
    t[0] = 0.5*SHIFT;
    t[1] = -0.5*SHIFT;
    //t[2] = 0.0;

    for(j=0;j<DIM;j++) center[j]= 0.0;
    for(i=0;i<NUM;i++) {
      for(j=0;j<DIM;j++) {
        center[j] += data1[i][j];
      }
    }
    for(j=0;j<DIM;j++) center[j] /= NUM;

    for(i=0;i<NUM;i++) {
      for(j=0;j<DIM;j++) {
        wrk[i][j] = data1[i][j] - center[j];
      }
    }
#if 0
    mxTrns(wrk,NUM,DIM,wrkT);
    mxMult(R,DIM,DIM,NUM,wrkT,data2T);
    mxTrns(data2T,DIM,NUM,data2);

    for(i=0;i<NUM;i++) {
      for(j=0;j<DIM;j++) {
        data2[i][j] = data2[i][j] + center[j];
      }
    }
#else
    mxMult(data1,NUM,DIM,DIM,R,data2);
#endif

    for(i=0;i<NUM;i++) {
      for(j=0;j<DIM;j++) {
        data2[i][j] += t[j];
        rnd = rand()/(double)RAND_MAX;
        data2[i][j] += SIG*comPnorm(rnd);
      }
    }
    
    for(j=0;j<DIM;j++) {
      if(j == 0) fprintf(fw1,"$");
      fprintf(fw1,"dim%d",j);
      if(j+1 < DIM) fprintf(fw1,",");
      else          fprintf(fw1,"\n");

      if(j == 0) fprintf(fw2,"$");
      fprintf(fw2,"dim%d",j);
      if(j+1 < DIM) fprintf(fw2,",");
      else          fprintf(fw2,"\n");
    }

    for(i=0;i<NUM;i++) {
      for(j=0;j<DIM;j++) {
        fprintf(fw1,"%lf",data1[i][j]);
        if(j+1 < DIM) fprintf(fw1,",");
        else          fprintf(fw1,"\n");

        fprintf(fw2,"%lf",data2[i][j]);
        if(j+1 < DIM) fprintf(fw2,",");
        else          fprintf(fw2,"\n");
      }
    }
    fclose(fw1);
    fclose(fw2);

    comMxFree((void **)data1,NUM,DIM);
    comMxFree((void **)data2,NUM,DIM);
    comMxFree((void **)R,DIM,DIM);
    comMxFree((void **)wrk,NUM,DIM);

    return(0);
}
#else
int main(argc,argv)
int argc;
char *argv[];
{
    //int ix1,ix2;
    int i,j;
    int n,m;
    int n1,n2;
    double **R,**Ra;
    double *t,*ta,*vec;
    double error,preError;
    int count;

    double **data1,**data2,**data2T,**data2R;

    FILE *fp1,*fp2,*fw,*fa;
    int *index;

    if(argc < 5) {
      fprintf(stderr,"Usage command mainDataFile ShiftDataFile OutFile\n");
      exit(-9);
    }
    if(!(fp1=fopen(argv[1],"r"))) {
      fprintf(stderr,"Cannot read master Data File=[%s]\n",argv[1]);
      exit(-1);
    }
    if(!(fp2=fopen(argv[2],"r"))) {
      fprintf(stderr,"Cannot read Shift Data File=[%s]\n",argv[2]);
      exit(-1);
    }
    if(!(fw=fopen(argv[3],"w"))) {
      fprintf(stderr,"Cannot write output Data File=[%s]\n",argv[3]);
      exit(-2);
    }
    if(!(fa=fopen(argv[4],"w"))) {
      fprintf(stderr,"Cannot write param Data File=[%s]\n",argv[3]);
      exit(-2);
    }

    readFile(fp1,&data1,&n1,&m);
    fclose(fp1);
    readFile(fp2,&data2,&n2,&m);
    fclose(fp2);

    data2T = (double **)comMalloc(sizeof(double *)*m);
    data2R = (double **)comMalloc(sizeof(double *)*m);
    for(j=0;j<m;j++) {
      data2T[j]=(double *)comMalloc(sizeof(double)*n2);
      data2R[j]=(double *)comMalloc(sizeof(double)*n2);
    }

    R = (double **)comMalloc(sizeof(double *)*m);
    Ra = (double **)comMalloc(sizeof(double *)*m);
    for(i=0;i<m;i++) {
      R[i] = (double *)comMalloc(sizeof(double)*m);
      Ra[i] = (double *)comMalloc(sizeof(double)*m);
    }
    t = (double *)comMalloc(sizeof(double)*m);
    ta= (double *)comMalloc(sizeof(double)*m);
    vec=(double *)comMalloc(sizeof(double)*m);

    index=(int *)comMalloc(sizeof(int)*(n1 > n2 ? n1:n2));

    /* 蓄積値の初期値 */
    for(i=0;i<m;i++) {
      Ra[i][i]=1.0;
      ta[i] = 0.0;
    }

    preError=0;
    count = 0;
    //count=0;//ループカウンタ
    //R=eye(2);//回転行列
    //t=zeros(2,1);//並進ベクトル

    while (1) {
	  count=count+1;
    
      error=FindNearestPoint(data1,n1,data2,n2,m,&n,index); //最近傍点探索
      /* 終了条件 */
      if(fabs(preError-error) < EPS) break;
      //if(fabs(error) < EPS) break;
      preError = error;

      SVDMotionEstimation(data1,data2,n,m,index,R,t);//特異値分解による移動量推定

      //計算したRとtで点群とRとtの値を更新
      mxTrns(data2,n,m,data2T);
#if 1
      mxMult(R,m,m,n,data2T,data2R);
      //data2=R1*data2;
      for(i=0;i<n;i++) {
        for(j=0;j<m;j++) {
          data2R[j][i] += t[j];
        }
      }
#else
      //data2=R1*data2;
      for(i=0;i<n;i++) {
        for(j=0;j<m;j++) {
          data2T[j][i] -= t[j];
        }
      }
      mxMult(R,m,m,n,data2T,data2R);
#endif
      mxTrns(data2R,m,n,data2);
      //ta2=[data2(1,:)+t1(1) ; data2(2,:)+t1(2)];
      mxMult(Ra,m,m,m,R,Ra);

      mxMult(Ra,m,m,m,R,Ra);
      mxVecR(R,m,m,ta,vec);
      for(j=0;j<m;j++) {
        ta[j] = vec[j] + t[j];
      }
      //R = R1*R;
      //t = R1*t + t1; 
    
      if(count > maxIter) {//収束しなかった
         break;
      }
    }
    if(count < maxIter) fprintf(stderr,"Convergence:%d\n",count);
    else                fprintf(stderr,"Max Iteration");


    /***************
      結果の書出し
    ****************/
    for(j=0;j<m;j++) {
      if(j == 0) fprintf(fw,"$");
      fprintf(fw,"dim%d",j);
      if(j+1 < m) fprintf(fw,",");
      else         fprintf(fw,",");
    }
    fprintf(fw,"index\n");

      mxTrns(data1,n,m,data2T);
#if 1
      mxMult(Ra,m,m,n,data2T,data2R);
      //data2=R1*data2;
      for(i=0;i<n;i++) {
        for(j=0;j<m;j++) {
          data2R[j][i] += ta[j];
        }
      }
      mxTrns(data2R,m,n,data2);

#endif
    for(i=0;i<n;i++) {
      for(j=0;j<m;j++) {
        fprintf(fw,"%lf",data2[i][j]);
        if(j + 1 < m) fprintf(fw,",");
        else          fprintf(fw,",");
      }
      fprintf(fw,"%d\n",index[i]);
    }
    fclose(fw);

    /* write out parameter */
    for(j=0;j<m;j++) {
      if(!j) fprintf(fa,"$");
      fprintf(fa,"Rota%d,",j);
    }
    fprintf(fa,"Sift\n");


    for(i=0;i<m;i++) {
      for(j=0;j<m;j++) {
        fprintf(fa,"%lf,",Ra[i][j]);
      }
      fprintf(fa,"%lf\n",ta[i]);
    }
    fclose(fa);

    /* 領域の開放 **/
    for(i=0;i<m;i++) {
      free(R[i]);
      free(Ra[i]);
    }
    free(R);
    free(Ra);
    free(t);
    free(ta);
    free(vec);
    free(index);

    return(0);
}
#endif
/************

*************/
int readFile(fp,pData,pN,pM)
FILE *fp;
double ***pData;
int *pN;
int *pM;
{
    int i,j;
    int n,m;
    char record[2048],*pc;
    double **data;

    rewind(fp);
    n=0;
    m=0;
    while(fgets(record,2048,fp)) {
      record[strlen(record)-1]='\0';
      if(record[0] == '#') continue;
      if(record[0] == '$') continue;
      if(record[0] == '&') continue;

      pc=strtok(record,", ");
      j=0;
      while(pc) {
        j++;
        pc=strtok(NULL,", ");
      }
      if(m < j) m=j;
      n++;
    }

    /* 領域の作成 */
    data=(double **)comMalloc(sizeof(double *)*n);
    for(i=0;i<n;i++) {
      data[i]=(double *)comMalloc(sizeof(double)*m);
    }
    rewind(fp);

    i=0;
    j=0;
    while(fgets(record,2048,fp)) {
      record[strlen(record)-1]='\0';
      if(record[0] == '#') continue;
      if(record[0] == '$') continue;
      if(record[0] == '&') continue;

      pc=strtok(record,", ");
      j=0;
      while(pc) {
        data[i][j]=atof(pc);
        j++;
        pc=strtok(NULL,", ");
      }
      i++;
    }

    *pN = n;
    *pM = m;
    *pData = data;

    return(n);

}
/**********

***********/
int readFree(data,n,m)
double **data;
int n;
int m;
{
    int i;
    for(i=0;i<n;i++) {
      free(data[i]);
    }
    free(data);
    return(0);
}
/************
//function [index, error]=FindNearestPoint(data1, data2)
//data2に対するdata1の最近傍点のインデックスを計算する関数
*************/
double FindNearestPoint(data1,n1,data2,n2,m,pN,index)
double **data1; //n1*m
double **data2; //n2*m
int n1;
int n2;
int m;
int *pN;    //o num of data;
int index[];
{
    int i,j,k;
    int N;
    double *dist;
    double *dmin;
    //int *index;
    //double dm;
    //int is;
    double error;



    N = n1 < n2 ? n1 : n2;
    dist = (double *)comMalloc(sizeof(double)*N);
    dmin = (double *)comMalloc(sizeof(double)*N);
    //index= (int *)comMalloc(sizeof(int)*N);
    //index=[];
    //error=0;

    for(i=0;i<N;i++) {
      /* i番目との相手側全点での距離を算出 */
      for(k=0;k<N;k++) dist[k]=0.0;
      for(k=0;k<N;k++) {
        for(j=0;j<m;j++) {
          dist[k] += pow(data1[i][j]-data2[k][j],2.0);
        }
      }
      /* i番目と最小距離の算出 */
      dmin[i]=dist[0];
      for(k=0;k<N;k++) {
        if(dmin[i] >= dist[k]) {
          dmin[i] = dist[k];
          index[i] = k;
        }
      }
    }
    /* 相互の最小点を算出 
    dm = 0;
    for(i=0;i<N;i++) {
      if(dm < dmin[i]) {
        is=i;
        dm = dmin[i];
      }
    }

    error=0;
    for(i=0;i<N;i++) {
      for(k=0;k<N;k++) {
        for(j=0;j<m;j++) {
          error += pow(data1[i][j]-data2[k][j],2.0);
        }
      }
    }
    */


    error=0;
    for(i=0;i<N;i++) error += dmin[i];

    *pN = N;
    //dx=(data2-repmat(data1(:,i),1,m2));
    //dist=sqrt(dx(1,:).^2+dx(2,:).^2);
    //[dist, ii]=min(dist);
    //index=[index; ii];
    //error=error+dist;
    free(dmin);
    free(dist);
    return(error);
}

/***
 特異値分解法による並進ベクトルと、回転行列の計算
***/
int SVDMotionEstimation(data1,data2,n,m,index,R,t)
double **data1;  //n*m
double **data2;  //n*m
int n;    /* num of data1 and data2*/
int m;    /* dimentison */
int *index;
double **R;
double *t;
{
  //function [R, t]=SVDMotionEstimation(data1, data2, index)
 

  double *mm,*ms,*vec; //m
  int i,j;
  double **Ssifted,**Msifted; //m*n (Quation!!)
  double **W,**TR,**v,**s,**H;  //m*m
  double **Trn;
  double det;

  mm =  (double *)comMalloc(sizeof(double)*m);
  ms =  (double *)comMalloc(sizeof(double)*m);
  vec = (double *)comMalloc(sizeof(double)*m);

  Ssifted = (double **)comMalloc(sizeof(double *)*m);
  Msifted = (double **)comMalloc(sizeof(double *)*m);

  W       = (double **)comMalloc(sizeof(double *)*m);
  //R       = (double **)comMalloc(sizeof(double *)*m);
  TR      = (double **)comMalloc(sizeof(double *)*m);

  /* svd arg */
  v       = (double **)comMalloc(sizeof(double *)*m);
  s       = (double **)comMalloc(sizeof(double *)*m);
  H       = (double **)comMalloc(sizeof(double *)*m);
  Trn     = (double **)comMalloc(sizeof(double *)*n);

  for(i=0;i<m;i++) {
	  Ssifted[i] = (double *)comMalloc(sizeof(double)*n);
	  Msifted[i] = (double *)comMalloc(sizeof(double)*n);
  }

  for(i=0;i<m;i++) {
	W[i]       = (double *)comMalloc(sizeof(double)*m);
    TR[i]      = (double *)comMalloc(sizeof(double)*m);
    //R[i]       = (double *)comMalloc(sizeof(double)*m);
    v[i]       = (double *)comMalloc(sizeof(double)*m);
    s[i]       = (double *)comMalloc(sizeof(double)*m);  
    H[i]       = (double *)comMalloc(sizeof(double)*m); 
  }

    for(i=0;i<n;i++) {
    Trn[i]     = (double *)comMalloc(sizeof(double)*m);
  }

  //各点群の重心の計算
  for(i=0;i<n;i++) {
	for(j=0;j<m;j++) {
	  mm[j] += data1[i][j];
    }
  }
  for(j=0;j<m;j++) mm[j] /= n;


  for(i=0;i<n;i++) {
    for(j=0;j<m;j++) {
#if 0
      ms[j] += data2[index[i]][j];
#else
      ms[j] += data2[i][j];
#endif
    }
  }
  for(j=0;j<m;j++) ms[j] /= n;


  //M = data1; 
  //mm = mean(M,2);
  //S = data2(:,index);
  //ms = mean(S,2); 

  //各点群を重心中心の座標系に変換
  for(j=0;j<m;j++) {
    for(i=0;i<n;i++) {
      Ssifted[j][i] = data1[i][j] - mm[j];
      Msifted[j][i] = data2[i][j] - ms[j];
    }
  }	

  //Sshifted = [S(1,:)-ms(1); S(2,:)-ms(2);];
  //Mshifted = [M(1,:)-mm(1); M(2,:)-mm(2);];

  mxTrns(Msifted,m,n,Trn);
  //W = Sshifted*Mshifted';
  mxMult(Ssifted,m,n,m,Trn,W);

  /****************************************/
  /*[U,A,V] =*/ svdcomp(W,m,m,vec,v);//特異値分解
  /****************************************/

  mxMult(W,m,m,m,v,TR);
  mxTrns(TR,m,m,R);
  //R = (U*V')';//回転行列の計算
  det = fabs(mxDet(R,m));

  /* 回転行列の制約 */
  for(j=0;j<m;j++) {
    H[j][j]=1.0;
  }
  H[m-1][m-1]=det;
  mxMult(W,m,m,m,H,s);
  mxMult(s,m,m,m,v,TR);
  mxTrns(TR,m,m,R);



  mxVecR(R,m,m,ms,vec);
  for(j=0;j<m;j++) {
    t[j] = mm[j] - vec[j];
  }
  //t = mm - R*ms;//並進ベクトルの計算

  /* 領域の開放 */
  free(mm);
  free(ms);
  free(vec);

  for(j=0;j<m;j++) {
    free(W[j]);
    free(v[j]);
    free(s[j]);
    free(TR[j]);
    //free(R[j]);

    free(Ssifted[j]);
    free(Msifted[j]);
    free(H[j]);
  }
  free(W);
  free(v);
  free(s);
  free(TR);
  //free(R);
  free(Ssifted);
  free(Msifted);

  for(i=0;i<n;i++) {
    free(Trn[i]);
  }
  free(Trn);
  free(H);

  return(0);


}



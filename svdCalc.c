/************************
  SVD���f��
*************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "usrlib.h"

#define DLM "\t"
int fprint_d();

extern FILE *g_pFlog;  /* calc.h */
int readParmKey();
int scanParmKey();
int checkParmKey();
int freeParmKey();
//#define ONLY


#define IMIN(ix,iy) ix > iy ? iy : ix
#define FMAX(x, y ) x  >  y ?  x :  y
#define SQR(x) pow(x,2.0)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define LIMIT 30
#define EPS 0.000001
#define MAX_VAL 128

#ifndef ONLY
extern FILE *pFclc;
int is_blank();
int is_comment();
int nki_blank();
#endif

static double pythag();
int make_svd();
int svdcmp();
int mxRevGJ();

/*******************
  ��Main
********************/
#ifdef ONLY

int main(argc,argv)
int argc;
char *argv[];
{
   int i,j,ret;
   int n,m;
   char **names;

   FILE *fp,*fd,*fa;

   char record[2048],*pc;

   double **data;  /* ���̓f�[�^ */
   
   if(argc < 3) {
     fprintf(stderr,"Usage:command infile outfile\n");
     exit(-1);
   }

   if(!(fp=fopen(argv[1],"r"))) {
     fprintf(stderr,"Cannot read file=[%s]\n",argv[1]);
     exit(-2);
   }
   if(!(fd=fopen(argv[2],"w"))) {
     fprintf(stderr,"Cannot write file=[%s]\n",argv[1]);
     exit(-3);
   }

   strcpy(record,argv[2]);
   record[strlen(record)-4]='\0';
   strcat(record,".txt");
   if(!(fa=fopen(record,"w"))) {
     fprintf(stderr,"Cannot write file=[%s]\n",record);
     exit(-3);
   }

   m=0;
   n=0;
   while(fgets(record,2048,fp)) {
     record[strlen(record)-1]='\0';
     if(record[0] == '#') continue;
     if(record[0] == '&') continue;
     if(record[0] == '$') continue;
     
     pc=comStrtok(record,",");
     j=0;
     while(pc) {
       pc=comStrtok(NULL,",") ;
       j++;
     }
     if(m < j) m=j;
     n++;
  }
  

  /* �̈�̊l�� */
  data=(double **)malloc(sizeof(double *)*(n));

  for(i=0;i<n;i++) {
    data[i]=(double *)malloc(sizeof(double)*(m));
    memset(data[i],'\0',sizeof(double)*(m));
  }

  names=(char **)malloc(sizeof(char *)*(m+1));
   
  rewind(fp);
   i=0;
   while(fgets(record,2048,fp)) {
     record[strlen(record)-1]='\0';
     if(record[0] == '#') continue;
     if(record[0] == '&') continue;
     if(record[0] == '$') {
       pc=comStrtok(record,",");
       pc++;
       j=0;
       while(pc) {
         names[j]=(char *)malloc(sizeof(char)*(strlen(pc)+1));
         strcpy(names[j],pc);
         pc=comStrtok(NULL,",");
         j++;
       }
       continue;
     } 
     pc=comStrtok(record,",");
     j=0;
     while(pc) {
       data[i][j]=atof(pc);
       pc=comStrtok(NULL,",") ;
       j++;
     }
     i++;
  } 
  fclose(fp);
  
  ret=make_svd(data,n,m,names,0,m,0,fa,fd);
  
  fclose(fd);
  fclose(fa);

  for(i=0;i<n;i++) {
    free(data[i]);
  }
  free(data);
  for(j=0;j<m;j++) {
    free(names[j]);
  }
  free(names);
  
  return(0);

}
#endif
/******************
   �����t�s��
   Calcer�Ƃ�I/F
*******************/
int make_revmx(data,n,m,names,ntag,nexp,parm,fa,fw)
double **data; /* n�sm�� */
int n,m;
char *names[];
int ntag;  /* �ړI�ϐ��̐� */
int nexp;  /* �����ϐ��̐� */
int parm; 
FILE *fa;
FILE *fw;
{

   int i,j,ret;
   double **ans;

   int mR,nObj;
   double sum,*vec;
   double **wrk;

   /* ��͑Ώې� */
   if(!nexp) nObj=0;  /* ��͑ΏۂȂ� */
   else {
	 nObj = m - nexp;
   }
   if(nObj > 1) return(-1);  /* �ړI�ϐ����� */
   mR=m-nObj;

   if(n != mR) return(-1); /* �����s��̃`�F�b�N */
   if(!n)      return(-1);

   
   ans=(double **)malloc(sizeof(double)*n);
   for(i=0;i<n;i++) {
	 ans[i]=(double *)malloc(sizeof(double)*n);
   }
   vec=(double *)malloc(sizeof(double)*n);
   wrk=(double **)malloc(sizeof(double)*n);

   /* �K�E�X�E�W�����_���@ */
   for(i=0;i<n;i++) {
	 wrk[i]=&data[i][nObj];
   }
   ret=mxRevGJ(wrk,mR,ans);
#if 0
   /* �^�C�g�� */
   for(j=0;j<mR;j++) {
     if(!j) fprintf(fw,"$");
     fprintf(fw,"%s",names[j+nObj]);
     if(j < mR-1) fprintf(fw,",");
     else         fprintf(fw,"\n");
   }

   if(nObj) {
	 for(i=0;i<mR;i++) {
	   sum=0;
	   for(j=0;j<mR;j++) {
	     sum += data[j][0] * ans[i][j];
       }
	   vec[i]=sum;
     }
     for(j=0;j<mR;j++) {
	   fprintf(fw,"%lf",vec[j]);
	   if(j < mR-1) fprintf(fw,",");
	   else         fprintf(fw,"\n");
	 }
   }

   for(i=0;i<mR;i++) {
     for(j=0;j<mR;j++) {
       fprintf(fw,"%lf",ans[i][j]);
	   if(j < mR-1) fprintf(fw,",");
	   else         fprintf(fw,"\n");
	 }	 
   }
#else
   /* �^�C�g�� */
   for(j=0;j<mR;j++) {
     if(!j) fprintf(fw,"$");
     fprintf(fw,"%s",names[j+nObj]);
     if(j < mR-1) fprintf(fw,DLM);
     else         fprintf(fw,"\n");
   }

   if(nObj) {
	 for(i=0;i<mR;i++) {
	   sum=0;
	   for(j=0;j<mR;j++) {
	     sum += data[j][0] * ans[i][j];
       }
	   vec[i]=sum;
     }
     for(j=0;j<mR;j++) {
	   fprint_d(fw,vec[j]);
	   if(j < mR-1) fprintf(fw,DLM);
	   else         fprintf(fw,"\n");
	 }
   }

   for(i=0;i<mR;i++) {
     for(j=0;j<mR;j++) {
       fprint_d(fw,ans[i][j]);
	   if(j < mR-1) fprintf(fw,DLM);
	   else         fprintf(fw,"\n");
	 }	 
   }
#endif

   /* �̈�̊J�� */
   for(i=0;i<n;i++) {
     free(ans[i]);
   }
   free(ans);
   free(vec);
   free(wrk);

   return(n);
}
/******************
   SVD
   Calcer�Ƃ�I/F
*******************/
int make_svd(data,n,m,names,ntag,nexp,parm,fa,fw)
double **data; /* n�sm�� */
int n,m;
char *names[];
int ntag;  /* �ړI�ϐ��̐� */
int nexp;  /* �����ϐ��̐� */
int parm; 
FILE *fa;
FILE *fw;
{
   double **a;
   double *w;
   double **v;
   double **rev;
   double **wrk;
   int i,j,k,ret;
   double sum;


  char record[256],*pc;
  char *parms[MAX_VAL];    /* �I�v�V�������Q */
  char *parmv[MAX_VAL];    /* �I�v�V�����l�Q */
  int np,line;
  int outkind;
  char buff[256];

  int nout;
  int mR;
  int nObj;  /* ��͑Ώەϐ��� */

  outkind=0;
  np=0;


  /* �s���̃`�F�b�N */
  if(n < m) return(-3);

  /* ��͑Ώې� */
  if(!nexp) nObj=0;  /* ��͑ΏۂȂ� */
  else {
	nObj = m - nexp;
  }
  if(nObj > 1) return(-1);  /* �ړI�ϐ����� */
 
#ifndef ONLY
  if(parm) {
    line=0;
    while(1) {
      fprintf(stdout,"%2d>>",++line);
      pc=fgets(record,256,pFclc);
      if(!pc) {
        /* �p�����[�^���r���łȂ� */
        return(-1);
      }
      if(is_blank(pc) || is_comment(pc)) continue;

      if(pc == NULL || record[0] == ';') break;
      pc=strtok(record,":=");
      if(pc) {
        strcpy(buff,pc);
        nki_blank(buff);
        parms[np]=(char *)malloc(strlen(buff)+1);
        strcpy(parms[np],buff);
        if((pc=strtok(NULL,";\n"))) {
          strcpy(buff,pc);
          nki_blank(buff);
          parmv[np]=(char *)malloc(strlen(buff)+1);
          strcpy(parmv[np],buff);
        }
        else parmv[np]=NULL;
        np++;
      }
    }
  }

  readParmKey("svd");
  for(i=0;i<np;i++) {
    if(scanParmKey(g_pFlog,"svd",parms[i]) < 0) continue;
    if(checkParmKey(g_pFlog,"svd",parms[i],parmv[i]) < 0) continue;

	if(!strcmp(parms[i],"output")) {
	  if(!strcmp(parmv[i],"rev")) outkind=0;
	  if(!strcmp(parmv[i],"uwv")) outkind=1;
	}
  }
  freeParmKey("svd");
   
#endif
   /* �̈�̊m�� */
   a=(double **)malloc(sizeof(double *)*(n+1));
   for(i=0;i<=n;i++) {
     a[i]=(double *)malloc(sizeof(double)*(m+1));
     memset(a[i],'\0',sizeof(double)*(m+1));
   }
   w=(double *)malloc(sizeof(double)*(m+1));
   v=(double **)malloc(sizeof(double *)*(m+1));
   for(i=0;i<=m;i++) {
     v[i]=(double *)malloc(sizeof(double)*(m+1));
   }
   /* a��FORTRAN�̔z�� */
   for(i=0;i<n;i++) {
     for(k=0,j=nObj;j<m;j++) {
	   a[i+1][k+1]=data[i][j];
	   k++;
     }
   }
   mR=k;
   
   /****************
    ���ْl���� SVD
   *****************/
   ret=svdcmp(a,n,mR,w,v);
   
   /* ���ʂ̕\���@*/
   fprintf(fa,"** Singular Value Decompsition **\n");
   fprintf(fa,"name        Singular Value \n");
   for(j=0;j<mR;j++) {
	 fprintf(fa,"%-12s %lf\n",names[j+nObj],w[j+1]);
   }

   rev=(double **)malloc(sizeof(double *)*m);
   wrk=(double **)malloc(sizeof(double *)*m);
   for(j=0;j<m;j++) {
	 rev[j]=(double *)malloc(sizeof(double)*n);
	 wrk[j]=(double *)malloc(sizeof(double)*n);
   }

   /* SVD�̏��s���\������ */
   if(outkind == 1) {
     for(j=0;j<mR;j++) {
       if(!j) fprintf(fw,"$type\t");
	   fprintf(fw,"%s",names[nObj+j]);
	   if(j < mR-1) fprintf(fw,DLM);
	   else         fprintf(fw,"\n");
	 }
     for(j=0;j<mR;j++) {
	   if(!j) fprintf(fw,"w\t");
	   fprint_d(fw,w[1+j]);
	   if(j < mR-1) fprintf(fw,DLM);
	   else         fprintf(fw,"\n");
	 }
     for(i=0;i<mR;i++) {
       for(j=0;j<mR;j++) {
	     if(!j) fprintf(fw,"v\t");
         fprint_d(fw,v[i+1][j+1]);
	     if(j < mR-1) fprintf(fw,DLM);
	     else         fprintf(fw,"\n");
	   }	 
	 }
     for(i=0;i<n;i++) {
	   for(j=0;j<mR;j++) {
	     if(!j) fprintf(fw,"u\t");
         fprint_d(fw,a[i+1][j+1]);
	     if(j < mR-1) fprintf(fw,DLM);
	     else         fprintf(fw,"\n");
	   }	 
	 }
	 nout=1+mR+n;
   }
   else 
   if(outkind == 0) {
     /* ���^�s��̋t�s����Z�o����@*/
     for(i=0;i<mR;i++) {
       for(j=0;j<n;j++) {
         if(fabs(w[i+1]) > EPS) wrk[i][j] = a[j+1][i+1] / w[i+1];  /* �]�uU�s��/���ْl�Ίp�s�� */
	     else                   wrk[i][j] = 0.0;
	   }
	 }

     for(i=0;i<mR;i++) {
	   for(j=0;j<n;j++) {
         sum=0;
         for(k=0;k<mR;k++) {
		   sum += v[i+1][k+1]*wrk[k][j];
		 }
	     rev[i][j]=sum;
	   }
	 }
	 /* �ړI�ϐ����� ���̎Z�o rev * tag */
	 if(nObj) {
	   for(i=0;i<mR;i++) {
		 sum=0;
		 for(j=0;j<n;j++) {
		   sum += data[j][0] * rev[i][j];
         }
		 wrk[i][0]=sum;
       }
       for(j=0;j<mR;j++) {
         if(!j) fprintf(fw,"$");
	     fprintf(fw,"%s",names[nObj+j]);
	     if(j < mR-1) fprintf(fw,DLM);
	     else         fprintf(fw,"\n");
	   }
       for(j=0;j<mR;j++) {
	     fprint_d(fw,wrk[j][0]);
	     if(j < mR-1) fprintf(fw,DLM);
	     else         fprintf(fw,"\n");
	   }
	   return(2);
     }
     else {
       for(j=0;j<n;j++) {
	     if(!j) fprintf(fw,"$");
	     fprintf(fw,"col%d",j+1);
	     if(j < n -1) fprintf(fw,DLM);
	     else         fprintf(fw,"\n");
	   }
       for(i=0;i<mR;i++) {
	     for(j=0;j<n;j++) {
           fprint_d(fw,rev[i][j]);
	       if(j < n-1) fprintf(fw,DLM);
	       else        fprintf(fw,"\n");
		 }	 
	   }
	   nout=mR;
     }
   }
   
   /* �̈�̊J�� */
   for(i=0;i<n+1;i++) {
     free(a[i]);
   }
   free(a);
   for(i=0;i<m+1;i++) {
     free(v[i]);
   }
   free(v);
   free(w);
   for(i=0;i<m;i++) {
	 free(rev[i]);
	 free(wrk[i]);
   }
   free(rev);
   free(wrk);

   for(i=0;i<np;i++) {
	 free(parms[i]);
	 free(parmv[i]);
   }
   return(nout);

}
/******************
  SVD����
  ���j
  �f�[�^��m�sn��i���ʂƋt)
  a,w,v��FORTRAN�̔z��
  ���p��)
  Numerical Recipe
*******************/
int svdcmp(a,m,n,w,v)
double **a; /* ���͑Ώۍs�� */
int n;  /* �� */
int m;  /* �s�� */
double *w; /* ���ْl�Ίp���� m�� */
double **v; /* m�sm�� */
{
   int flag,i,its,j,jj,k,l,nm;
   double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

   rv1=(double *)malloc(sizeof(double)*(n+1));
   memset(rv1,'\0',sizeof(double)*(n+1));

   g=scale=anorm=0.0;
   /* HouseHolder�@�̂Q�d�Ίp�̌`�ɒ��� */
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				if (i != n) {
					for (j=l;j<=n;j++) {
						for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
						f=s/h;
						for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
					}
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				if (i != m) {
					for (j=l;j<=m;j++) {
						for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
						for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
					}
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=n;i>=1;i--) {
		l=i+1;
		g=w[i];
		if (i < n)
			for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			if (i != n) {
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
					f=(s/a[i][i])*g;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else {
			for (j=i;j<=m;j++) a[j][i]=0.0;
		}
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if (fabs(rv1[l])+anorm == anorm) {
					flag=0;
					break;
				}
				if (fabs(w[nm])+anorm == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					if (fabs(f)+anorm != anorm) {
						g=w[i];
						h=pythag(f,g);
						w[i]=h;
						h=1.0/h;
						c=g*h;
						s=(-f*h);
						for (j=1;j<=m;j++) {
							y=a[j][nm];
							z=a[j][i];
							a[j][nm]=y*c+z*s;
							a[j][i]=z*c-y*s;
						}
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k]=(-v[j][k]);
				}
				break;
			}
			if (its == 30) {
			  fprintf(stderr,"No convergence in 30 SVDCMP iterations\n");
			  return(-its);
			}
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y=y*c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=(c*g)+(s*y);
				x=(c*y)-(s*g);
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	return(m);
}

/* (a**2 + b**2)**0.5�̌v�Z�@�I�[�o�t���[���N�����ɂ��� */
static double pythag(a,b)
double a;
double b;
{
   double absa,absb;
   absa = fabs(a);
   absb = fabs(b);
   if(absa > absb) return(absa*sqrt(1.0+SQR(absb/absa)));

   return(absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

#undef IMIN
#undef FMAX
#undef SIGN
#undef SQR
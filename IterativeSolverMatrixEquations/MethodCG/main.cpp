#include "mtxreader.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <iostream>
using namespace std;

bool isSrandCalled = false;
const double EPSILON = 0.0000001;
const int MAX_VAL = 1000;

struct crsMatrix{
  int N;// размер матрицы
  int NZ;// число ненулевых элементов
  double *Value;// массив значений
  int *Col;//массив номеров столбцов
  int *RowIndex;// массив указателей позиций, с которых начинается описание очередной строки в массивах Value и Col
};

void InitMatrix(int N, int NZ, crsMatrix &mtx){
  mtx.N = N;
  mtx.NZ = NZ;
  mtx.Value = new double[NZ];
  mtx.Col = new int[NZ];
  mtx.RowIndex = new int[N + 1];
}

void FreeMatrix(crsMatrix &mtx){
  delete[] mtx.Value;
  delete[] mtx.Col;
  delete[] mtx.RowIndex;
}
int WriteMatrix(crsMatrix mtx, char *fileName){
  FILE *f = fopen(fileName, "w+");
  fprintf(f, "%i\n", mtx.NZ);
  fprintf(f, "%i\n", mtx.N);
  for (int i = 0; i < mtx.NZ; i++)
    fprintf(f, "%lf;%i\n", mtx.Value[i], mtx.Col[i]);
  for (int i = 0; i < mtx.N + 1; i++)
    fprintf(f, "%i\n", mtx.RowIndex[i]);
  fclose(f);
  return 0;
}

int ReadMatrix(crsMatrix &mtx, char *fileName){
  int N, NZ;
  if (fileName == NULL) return -1;
  FILE *f = fopen(fileName, "r");
  if (f == NULL) return -1;
  fscanf(f, "%d", &NZ);
  fscanf(f, "%d", &N);
  InitMatrix(N, NZ, mtx);
  for (int i = 0; i < NZ; i++)
    fscanf(f, "%lf;%i", &(mtx.Value[i]), &(mtx.Col[i]));
  for (int i = 0; i < N + 1; i++)
    fscanf(f, "%d", &(mtx.RowIndex[i]));
  fclose(f);
  return 0;
}
int ReadFromMTX(crsMatrix &mtx, char *fileName){
  int N, NZ;
  if (fileName == NULL) return -1;
  FILE *f = fopen(fileName, "r");
  if (f == NULL) return -1;
  char str[100];
  //int c;
  while(getc(f)=='%')
    fgets(str,100,f);
  fscanf(f, "%d", &N);
  fscanf(f, "%d", &N);
  fscanf(f, "%d", &NZ);
  int* row = new int[NZ];
  int* col = new int[NZ];
  double* val = new double[NZ];
  int nz_diag=0;
  for (int i = 0; i < NZ; i++) {
    fscanf(f, "%i %i %lf", &(col[i]),&(row[i]), &(val[i]));
    if (col[i]==row[i])
      nz_diag++;
  }
  InitMatrix(N,NZ*2-nz_diag,mtx);
  int ind=0;
  int iind=0;
  int k=0;
  mtx.RowIndex[0]=0;
  for (int i=0; i<NZ; i++){
    if (row[i]==1){
      mtx.Col[k]=col[i]-1;
      mtx.Value[k]=val[i];
      ind++;
      k++;
    }
    else{
      if (row[i]!=row[i-1]){
        iind++;
        mtx.RowIndex[iind]=mtx.RowIndex[iind-1]+ind;
        ind=0;
        for (int t=0;t<i;t++){
          if (col[t]==row[i]){
            mtx.Col[k]= row[t]-1;
            mtx.Value[k]=val[t];
            k++;
            ind++;
          }
        }
        mtx.Col[k]=col[i]-1;
        mtx.Value[k]= val[i];
        ind++;
        k++;
      }
      else{
        mtx.Col[k]=col[i]-1;
        mtx.Value[k]= val[i];
        ind++;
        k++;
      }
    }
  }
  mtx.RowIndex[iind+1]=mtx.RowIndex[iind]+ind;
  delete[] row;
  delete[] col;
  delete[] val;
  return 0;
}

int ReadEps(double &eps, char *fileName){
  if (fileName == NULL) return -1;
  FILE *f = fopen(fileName, "r");
  if (f == NULL) return -1;
  fscanf(f, "%lf", &eps);
  fclose(f);
  return 0;
}

void MultMV(crsMatrix M,const double *v, double *res){
  for (int i = 0; i < M.N; i++){
    res[i] = 0.0;
    for (int j=M.RowIndex[i]; j<M.RowIndex[i + 1]; j++)
      res[i] += M.Value[j] * v[M.Col[j]];
  }
}
double ScalarProduct(const double *v1,const double *v2, int size){
  double result = 0;
  for (int i = 0; i < size; i++)
    result += v1[i]*v2[i];
  return result;
}

double Dest(double *x1, double *x2, int size){
  double norm = 0;
  for(int i=0; i<size; i++)
    norm += (x1[i]-x2[i]) * (x1[i]-x2[i]);
  return sqrt(norm);
}
void PrintVec(double* v, int size){
  for(int i=0; i<size;i++)
    cout<<v[i]<<" ";
  cout<<endl;
}
void Method_CG_1(crsMatrix A, double *b, double eps, double *result, int &iter){
  int size = A.N;
  double * currX = new double[size];
  double * prevX = new double[size];
  double * currZ = new double[size];
  double * prevZ = new double[size];
  double * currR = new double[size];
  double * prevR = new double[size];
  double * Az = new double[size];
  double alpha, betta;
  iter = 1;
  int maxIter = size * size;
  for (int i=0; i<size; i++) {
    prevX[i] = 0;
    prevZ[i] = b[i];
    prevR[i] = b[i];
  }
  do{
    if (iter > 1){//текущий делаем предидущим
      memcpy(prevX,currX,size * sizeof(double));
      memcpy(prevZ,currZ,size * sizeof(double));
      memcpy(prevR,currR,size * sizeof(double));
    }
    MultMV(A,prevZ,Az);
    alpha = ScalarProduct(prevR,prevR,size)/ScalarProduct(Az,prevZ,size);
    for (int i=0; i<size; i++){
      currX[i] = prevX[i] + alpha * prevZ[i];
      currR[i] = prevR[i] - alpha * Az[i];
    }
    betta = ScalarProduct(currR,currR,size)/ScalarProduct(prevR,prevR,size);
    for (int i=0; i<size; i++)
      currZ[i] = currR[i] + betta * prevZ[i];
    iter++;
  }while
    ((Dest(prevX, currX, size) > eps) && (iter < maxIter));
  for (int i=0; i<size; i++)
    result[i] = currX[i];
  delete[] currX;
  delete[] prevX;
  delete[] currZ;
  delete[] prevZ;
  delete[] currR;
  delete[] prevR;
  delete[] Az;
}
void Method_CG(crsMatrix A, double *b, double eps, double *result, int &iter){
  int size = A.N;
  double * currX = new double[size];
  double * prevX = new double[size];
  double * currGrad = new double[size];
  double * prevGrad = new double[size];
  double * currDir = new double[size];
  double * prevDir = new double[size];
  double * temp = new double[size];
  double step;
  iter = 1;
  int maxIter = size *size;
  for (int i=0; i<size; i++) {
    prevX[i] = 0;
    prevDir[i] = 0;
    prevGrad[i] = -b[i];
  }
  do{
    if (iter > 1){//текущий делаем предидущим
      memcpy(prevX,currX,size * sizeof(double));
      memcpy(prevGrad,currGrad,size * sizeof(double));
      memcpy(prevDir,currDir,size * sizeof(double));
    }

    //вычисляем градиент
    MultMV(A, prevX, currGrad);
    for(int i=0; i<size; i++)
      currGrad[i] -= b[i];

    //вычисляем направление
    double coeff = ScalarProduct(currGrad, currGrad, size)/ScalarProduct(prevGrad, prevGrad, size);
    for (int i=0; i<size; i++)
      currDir[i] = -currGrad[i] + prevDir[i] * coeff;

    //вычисление величины смещения по выбранному направлению
    MultMV(A, currDir, temp);
    step = ScalarProduct(currGrad, currGrad, size)/ScalarProduct(temp, currDir, size);

    //вычисление нового приближения решения
    for (int i=0; i<size; i++)
      currX[i] = prevX[i] + step * currDir[i];

    iter++;
  }while
    ((Dest(prevX, currX, size) > eps) && (iter < maxIter));
  for (int i=0; i<size; i++)
    result[i] = currX[i];
  delete[] currX;
  delete[] prevX;
  delete[] currGrad;
  delete[] prevGrad;
  delete[] currDir;
  delete[] prevDir;
  delete[] temp;
}
double Next(){
  return ((double)rand() / (double)RAND_MAX);
}
void GenerateVector(int seed, int N, double *vec){
  srand(seed);
  for (int i = 0; i < N; i++)
    vec[i] = Next() * MAX_VAL;
}
void KeyboardMatrixInit(crsMatrix &mtx){
  int n,nz;
  cout<< "N= ";
  cin>> n;
  cout<< "NZ= ";
  cin>> nz;
  InitMatrix(n,nz,mtx);
  for(int i=0; i<nz; i++){
    cout<<"col["<<i+1<<"]= ";
    cin>>mtx.Col[i];
  }
  for(int i=0; i<nz; i++){
    cout<<"val["<<i+1<<"]= ";
    cin>>mtx.Value[i];
  }
  for(int i=0; i<n+1; i++){
    cout<<"rowindex["<<i+1<<"]= ";
    cin>>mtx.RowIndex[i];
  }
}
void main(int argc, char* argv[]){
  crsMatrix A;
  double eps = EPSILON;
  char * outputfile = NULL;
  char * timefile = NULL;
  if (argc == 5){
    ReadMatrix(A, argv[1]);
    ReadEps(eps,argv[2]);
    outputfile = argv[3];
    timefile = argv[4];
  }
  //KeyboardMatrixInit(A);
  //ReadFromMTX(A,"tmt_sym.mtx");
  ReadMatrixFromFile("tmt_sym.mtx", &A.N, &A.Col, &A.RowIndex, &A.Value);
  cout<<"reding is OK"<<endl;
  //InitMatrix(n,nz,A);
  A.NZ = A.RowIndex[A.N];
  cout<<"NZ = "<<A.RowIndex[A.N]<<endl;
  cout<<"N = "<<A.Col[A.NZ-1]+1<<endl;
  double * x = new double[A.N]; //эталонное, точное решение
  double * b = new double[A.N]; //b=Ax
  double * result = new double[A.N];// результат
  int iter;

  for(int i=0; i<A.N;i++)
    x[i]=1;
  MultMV(A,x,b);
  Method_CG(A, b, eps, result, iter);
  cout<<"iteration = "<<iter<<endl;
  //for(int i=0; i<A.NZ;i++)
  //cout<<A.Col[i]<<" ";
  //cout<<endl;
  //for(int i=0; i<A.NZ;i++)
  //cout<<A.Value[i]<<" ";
  //cout<<endl;
  //for(int i=0; i<=A.N;i++)
  //cout<<A.RowIndex[i]<<" ";
  for(int i=0; i<A.N;i=i+A.N/10)
    cout<<result[i]<<" ";
  FreeMatrix(A);
  delete[] x;
  delete[] b;
  delete[] result;
}

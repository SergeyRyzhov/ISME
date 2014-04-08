#ifndef __MTX_READER__
#define __MTX_READER__

#define MTX_TYPE double
#define ISME_OK 0
#define ISME_ERROR 1

int ReadMatrixFromFile(char* fileName, int* n, int** column, int** row, MTX_TYPE** val);

#endif
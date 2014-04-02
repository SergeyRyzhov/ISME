#include <stdio.h> 
#include <stdlib.h>  
#include <string.h>

#include "mtxreader.h"


typedef struct
{
  char * fileName;

}Parameters;
Parameters input;

#define empty 0

void ParseArguments(int argc, char** argv);
/*void Clean(Parameters self);
void Clean(Parameters self)
{
if (empty != self.fileName)
{
free(self.fileName);
self.fileName = empty;
}
}*/


void ParseArguments(int argc, char** argv)
{
  const int fileIndx = 1;

  if (argc > fileIndx)
  {
    input.fileName = (char*)malloc(strlen(argv[fileIndx]));
    strcpy(input.fileName, argv[fileIndx]);
  }

}

int main(int argc, char* argv[])
{
  int n;  
  int* col;
  int* row;
  double* val;

  ParseArguments(argc, argv);

  printf("%d\n", argc);
  printf("%s\n", input.fileName);
  ReadMatrixFromFile(input.fileName, &n, &col, &row,&val );
  printf("%d\n", n);

  //Clean(options);
  return 0;
}

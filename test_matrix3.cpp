#include <cstdio>

#include <core/matrix3.h>

int
main()
{
  matrix3<double> A;

  int i, j;
  // fill matrix
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      A(i, j) = i + j;
    }
  }

  printf("Matrix:\n");
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      printf("%f ", A(i, j));
    }
    printf("\n");
  }

  printf("Program ended normally\n");
  return 0;
}

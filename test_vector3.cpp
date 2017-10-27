#include <iostream>
#include <cstdio>
#include "jdftx/core/vector3.h"

using namespace std;

int main()
{
  // Using explicit constructor
  auto a = vector3<>(1.0, 2.0, 3.0); // can use default template parameter
  auto b = vector3<>(1.2, 2.1, 2.1);

  // access elements
  printf("a = [%f %f %f]\n", a[0], a[1], a[2]); // no operator() ??
  printf("b = [%f %f %f]\n", b[0], b[1], b[2]);

  // test operators
  auto c = a + b;
  printf("c = [%f %f %f]\n", c[0], c[1], c[2]);

  //
  printf("length a = %f\n", a.length());

  printf("Program ended normally.\n");
  return 0;
}

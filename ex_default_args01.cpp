#include <iostream>

using namespace std;

void my_func(int i, const char* prefix="Halo") {
  cout << prefix << ", this is number " << i << endl;
}

int main() {
  my_func(1, "Hello");
  my_func(3);
}
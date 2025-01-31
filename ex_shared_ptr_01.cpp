#include <iostream>
#include <memory>

using namespace std;

struct MyStruct {
  int x;
};

void func01(shared_ptr<MyStruct> p) {
  // change the underlaying object
  p->x = 20;
}

int main() {
  // initialize
  auto one = shared_ptr<MyStruct>(new MyStruct());
  // another shared_ptr pointing nowhere
  shared_ptr<MyStruct> two;

  // modify
  one->x = 10;

  // read
  cout << "x: " << one->x << endl;

  return 0;
}
#include <cstdio>

int main()
{
  int * buf = new int[10000];
  printf("%s", "allocated");
  delete buf;
}


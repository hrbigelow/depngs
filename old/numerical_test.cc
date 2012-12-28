#include <cmath>
#include <cstdio>
#include <limits>

int main(int argc, char ** argv)
{

    double u;
    sscanf(argv[1], "%lf", &u);

    double e = exp2l(u);
    
    printf("%10.10lf, exp2f: %lg\n", u, e);

    if (e == std::numeric_limits<double>::infinity())
    {
        printf("e is infinite");
    }
}

    

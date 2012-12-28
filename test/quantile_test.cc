#include <algorithm>
#include <vector>

#include "henry/sampling.h"

int main()
{
    size_t N = 1000;
    std::vector<int> ints(N);
    double nums[N];
    std::vector<double *> pnums(N);

    for (size_t i = 0; i != N; ++i)
    {
        nums[i] = static_cast<double>(i);
        pnums[i] = nums + i;
    }
    std::random_shuffle(pnums.begin(), pnums.end());

    for (size_t i = 0; i != N; ++i)
    {
        printf("%g ", *pnums[i]);
    }

    double quantiles[] = { 0.005, 0.05, 0.5, 0.95, 0.995 };
    double quantile_values[5];

    find_integral_bounds(& pnums, 0, quantiles, 5, quantile_values);

    return 0;
}

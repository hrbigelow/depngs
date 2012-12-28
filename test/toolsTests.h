#ifndef __TOOLSTEST_H
#define __TOOLSTEST_H

#include <cfloat>
#include <cstdio>

#include <cxxtest/TestSuite.h>
#include "tools.h"

class toolsTests : public CxxTest::TestSuite
{
public:

    // Test fixture setup code.
    void setUp()
    {
        const int size = 20;
        this->numbers.resize(size);
        this->running_total.resize(size);

        this->sum = 0.0;
        for (int i = 0; i != size; ++i)
        {
            float value = static_cast<float>(i);
            numbers[i] = log2f(value);
            sum += value;
            running_total[i] = log2f(sum);
        }
    }

    // Test fixture cleanup code.
    void tearDown()
    {

    }

    //test quality to quality code translation
    void test_QualityTranslation()
    {
        for (int i = 0; i != 50; ++i)
        {
            TS_ASSERT_EQUALS(i, QualityCodeToQuality(QualityToQualityCode(i)));
        }

        float const delta = 0.00000001;
        for (int i = 1; i != 7; ++i)
        {
            float error_prob = exp10(-static_cast<float>(i));
            int quality = i * 10;
            float delta = error_prob / 10000.0f;
            
            TS_ASSERT_DELTA(error_prob, QualityToErrorProb(quality), delta);
        }
    }
        
    //test SumLog2ArrayTruncated, Log2Accumulate, and SumLog2NumbersTruncated
    void test_SumLog2()
    {

        //TS_ASSERT_EQUALS(-INFINITY, -NAN);
        //test SumLog2Numbers extreme behavior
        TS_ASSERT_EQUALS(log2f(0.0), -INFINITY);
        //float values[] = { 0.0, FLT_MIN, 5.2, 52983984.23984f, FLT_MAX };
        float values[] = { FLT_MIN, 5.2, 52983984.23984f };
        int size = sizeof(values) / sizeof(float);
        for (int i = 0; i != size; ++i)
        {
            for (int j = i; j != size; ++j)
            {
                printf("Testing %f added to %f\n", log2f(values[i]), log2f(values[j]));
                TS_ASSERT_EQUALS(log2f(values[i] + values[j]),
                                 SumLog2NumbersTruncated(log2f(values[i]), 
                                                         log2f(values[j])));
            }
        }

    }



//    void test_SumLog2ArrayTruncated()
//    {
//        TS_ASSERT_EQUALS(log2f(sum), SumLog2Truncated(numbers));
//    }


//    void testLog2Accumulat()
//    {
//        std::vector<float> test_running_total =
//            Log2Accumulate(numbers.begin(), numbers.end());

        //!!! This is necessary to avoid the TS_ASSERT_EQUALS(inf, *) hanging bug
        //for (int i = 0; i != test_running_total.size(); ++i)
//        for (int i = 1; i != test_running_total.size(); ++i)
//        {
//            TS_ASSERT_DELTA(test_running_total[i],
//                            this->running_total[i],
//                            test_running_total[i] * 0.00001);
//        }
//    }

    class SquareFunction {
    public:
        float operator()(float x)
        {
            return x * x;
        }
    };

    class InverseSquareFunction {
    public:
        float operator()(float x)
        {
            return -(x * x);
        }
    };

    void testBinarySearchFunctional()
    {
        float x = 5.23984938;
        SquareFunction square;
        float y = square(x);
        float delta = 0.0000001;
        float x_bar = BinarySearchFunctional(0.0, 20.0, y, delta, delta, square);
        TS_ASSERT_DELTA(x_bar, x, delta);
        printf("xbar = %10.10f, x = %10.10f, x_bar - x = %10.10f, delta=%10.10f\n", x_bar, x, x_bar - x, delta);

        float x2 = x;
        InverseSquareFunction inv_square;
        float y2 = inv_square(x2);
        float x2_bar = BinarySearchFunctional(20.0, 0.0, y2, delta, delta, inv_square);
        TS_ASSERT_DELTA(x2_bar, x2, delta);
        printf("x2_bar = %10.10f, x2 = %10.10f, x2_bar - x2 = %10.10f, delta=%f10.10\n", x2_bar, x2, x2_bar - x2, delta);

    };


    //test with known analytical functions
    void testLowerUpperBoundNormLog2()
    {
        SAMPLE log2_ecurve;
        const size_t limit = 10000;
        for (size_t i = 0; i <= limit; ++i)
        {
            float x = static_cast<float>(i) / static_cast<float>(limit);
            log2_ecurve.insert(std::make_pair(x, log2f(expf(x))));
        }

        float delta = 0.0001;

        SAMPLE log2_ecurve_integral_lower = 
            LowerBoundNormLog2(log2_ecurve.begin(), log2_ecurve.end());

        TS_ASSERT_DELTA(log2f(expf(1.0) - 1.0), log2_ecurve_integral_lower[1.0], delta);

        SAMPLE log2_ecurve_integral_upper = 
            UpperBoundNormLog2(log2_ecurve.begin(), log2_ecurve.end());

        TS_ASSERT_DELTA(log2f(expf(1.0) - 1.0), log2_ecurve_integral_upper[1.0], delta);

    }


 private:
    std::vector<float> numbers, running_total;
    float sum;
    
};

#endif // __TOOLSTEST_H

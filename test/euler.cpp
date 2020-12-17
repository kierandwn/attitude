#include <gtest/gtest.h>

#include "include/helpers.h"
#include "euler.h"

using namespace attitude;

TEST(EulerInstantiation, InstantiateEuler)
{
    // Ensures ability to create concrete class.
    euler_<int> theta;
}

TEST(EulerInstantiation, DefaultIs0Rotation)
{
    // Zero method returns instatiates and returns identity.
    euler_<int> theta;
    ASSERT_TRUE(theta[0] == 0 && theta[1] == 0 && theta[2] == 0);
}

TEST(EulerDebugging, STDOUTDisplay)
{
    // Displays rotation matrix to console without error.
    euler_<int> theta;
    display(theta);
}

TEST(EulerMath, ConversionFromDCM)
{
    dcm_<double> R = R1(random_angle<double>()) + 
                     R2(random_angle<double>()) +
                     R3(random_angle<double>());
    euler_<double> theta(R, randomise_order());
    ASSERT_TRUE(theta == R);
}

TEST(EulerMath, ComparisonWithDCM)
{
    // Adding rotations.
    dcm_<double> R = R1(random_angle<double>());
    euler_<double> theta(R, randomise_order());
    ASSERT_TRUE(theta == R);
}

TEST(EulerMath, AddRotations)
{
    double alpha = random_angle<double>();
    double beta  = random_angle<double>();
    
    // Adding rotations.
    euler_<double> T1(R2(alpha), randomise_order());
    euler_<double> T2(R2(beta),  randomise_order());
    euler_<double> T3 = T1 + T2;
    ASSERT_TRUE(T3 == R2(alpha + beta));
}

TEST(EulerMath, SubtractRotations)
{
    double alpha = random_angle<double>();
    double beta = random_angle<double>();

    // Subtracting rotations.
    dcm_<double> T1(R3(alpha));
    dcm_<double> T2(R3(beta));

    dcm_<double> T3 = T1 - T2;
    ASSERT_TRUE(T3 == R3(alpha - beta));
}

TEST(EulerUtility, Reverse)
{
    double alpha = random_angle<double>();

    // Reversing DCM should be the same creating one in the negative direction.
    euler_<double> T1_pos(R1(alpha), randomise_order());
    euler_<double> T2_neg(R1(-1 * alpha), randomise_order());
    ASSERT_TRUE(T2_neg == T1_pos.reverse());
}


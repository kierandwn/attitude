#include <gtest/gtest.h>

#include "include/helpers.h"
#include "quaternion.h"

using namespace attitude;

TEST(QuaternionInstantiation, InstantiateQuaternion)
{
    // Ensures ability to create concrete class.
    quaternion_<int> q;
}

TEST(QuaternionInstantiation, DefaultIs0Rotation)
{
    // Zero method returns instatiates and returns identity.
    quaternion_<int> q;
    ASSERT_TRUE(q == ZERO<int>());
}

TEST(QuaternionDebugging, STDOUTDisplay)
{
    // Displays rotation matrix to console without error.
    quaternion_<int> q;
    display(q);
}

TEST(QuaternionMath, ConversionFromDCM)
{
    dcm_<double> R = R1(random_angle<double>()) +
        R2(random_angle<double>()) +
        R3(random_angle<double>());
    quaternion_<double> q(R);
    ASSERT_TRUE(q == R);
}

TEST(QuaternionMath, AddRotations)
{
    // Adding rotations.
    double alpha = random_angle<double>();
    double beta  = random_angle<double>();

    quaternion_<double> T1(R2<double>(alpha));
    quaternion_<double> T2(R2<double>(beta));
    quaternion_<double> T3 = T1 + T2;
    ASSERT_TRUE(T3 == R2(alpha + beta));
}

TEST(QuaternionMath, SubtractRotations)
{
    // Adding rotations.
    double alpha = random_angle<double>();
    double beta = random_angle<double>();

    // Subtracting rotations.
    quaternion_<double> T1(R2(alpha));
    quaternion_<double> T2(R2(beta));
    quaternion_<double> T3 = T1 - T2;
    ASSERT_TRUE(T3 == R2<double>(alpha - beta));
}

TEST(QuaternionUtility, Reverse)
{
    double alpha = random_angle<double>();

    // Reversing DCM should be the same creating one in the negative direction.
    quaternion_<double> T1_pos(R1(alpha));
    quaternion_<double> T2_neg(R1(-1 * alpha));
    ASSERT_TRUE(T2_neg == T1_pos.reverse());
}


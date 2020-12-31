#include <gtest/gtest.h>

#include "include/helper.h"
#include "crp.h"

using namespace attitude;

TEST(CRPInstantiation, InstantiateCRP)
{
  // Ensures ability to create concrete class.
  crp<double> q;
}

TEST(CRPInstantiation, DefaultIs0Rotation)
{
  // Zero method returns instatiates and returns identity.
  crp<double> q(0, 0, 0);
  ASSERT_TRUE(q == dcm::ZERO<double>());
}

TEST(CRPDebugging, STDOUTDisplay)
{
    // Displays rotation matrix to console without error.
    crp<double> q;
    // display(q); TODO
}

TEST(CRPMath, ConversionFromDCM)
{
  matrix<double, 3, 3> R = random_rotation<double>();
  crp<double> q(R);
  ASSERT_TRUE(q == R);
}

TEST(CRPMath, AddRotations)
{
    // Adding rotations.
    double alpha = random_angle<double>();
    double beta  = random_angle<double>();

    crp<double> T1(dcm::R2<double>(alpha));
    crp<double> T2(dcm::R2<double>(beta));
    crp<double> T3 = T1 + T2;
    ASSERT_TRUE(T3 == dcm::R2(alpha + beta));
}

TEST(CRPMath, SubtractRotations)
{
    // Adding rotations.
    double alpha = random_angle<double>();
    double beta = random_angle<double>();

    // Subtracting rotations.
    crp<double> T1(dcm::R2(alpha));
    crp<double> T2(dcm::R2(beta));
    crp<double> T3 = T1 - T2;
    ASSERT_TRUE(T3 == dcm::R2<double>(alpha - beta));
}

TEST(CRPUtility, Reverse)
{
    double alpha = random_angle<double>();

    // Reversing DCM should be the same creating one in the negative direction.
    crp<double> T1_pos(dcm::R1(alpha));
    crp<double> T2_neg(dcm::R1(-1 * alpha));
    ASSERT_TRUE(T2_neg == T1_pos.reverse());
}


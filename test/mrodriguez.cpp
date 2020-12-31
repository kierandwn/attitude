#include <gtest/gtest.h>

#include "include/helper.h"
#include "mrp.h"

using namespace attitude;

TEST(MRPInstantiation, InstantiateMRP)
{
  // Ensures ability to create concrete class.
  mrp<double> sigma;
}

TEST(MRPInstantiation, DefaultIs0Rotation)
{
  // Zero method returns instatiates and returns identity.
  mrp<double> q(0, 0, 0);
  ASSERT_TRUE(q == dcm::ZERO<double>());

}

TEST(MRPDebugging, STDOUTDisplay)
{
    // Displays rotation matrix to console without error.
    mrp<int> q;
    // display(q);
}

TEST(MRPMath, ConversionFromDCM)
{
  matrix<double, 3, 3> R = random_rotation<double>();
  mrp<double> q(R);
  ASSERT_TRUE(q == R);
}

TEST(MRPMath, AddRotations)
{
    // Adding rotations.
    double alpha = random_angle<double>();
    double beta  = random_angle<double>();

    mrp<double> T1(dcm::R2<double>(alpha));
    mrp<double> T2(dcm::R2<double>(beta));
    mrp<double> T3 = T1 + T2;
    ASSERT_TRUE(T3 == dcm::R2(alpha + beta));
}

TEST(MRPMath, SubtractRotations)
{
    // Adding rotations.
    double alpha = random_angle<double>();
    double beta = random_angle<double>();

    // Subtracting rotations.
    mrp<double> T1(dcm::R2(alpha));
    mrp<double> T2(dcm::R2(beta));
    mrp<double> T3 = T1 - T2;
    ASSERT_TRUE(T3 == dcm::R2<double>(alpha - beta));
}

TEST(MRPUtility, Reverse)
{
    double alpha = random_angle<double>();

    // Reversing DCM should be the same creating one in the negative direction.
    mrp<double> T1_pos(dcm::R1(alpha));
    mrp<double> T2_neg(dcm::R1(-1 * alpha));
    ASSERT_TRUE(T2_neg == T1_pos.reverse());
}


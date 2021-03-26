#include <gtest/gtest.h>

#include "include/helper.h"
#include "attitude/quaternion.h"

using namespace attitude;

TEST(QuaternionInstantiation, InstantiateQuaternion)
{
  // Ensures ability to create concrete class.
  quaternion<double> q;
}

TEST(QuaternionInstantiation, DefaultIs0Rotation)
{
  // Zero method and returns identity.
  quaternion<double> q(1, 0, 0, 0);
  ASSERT_TRUE(q == dcm::ZERO<int>());
}

TEST(QuaternionDebugging, STDOUTDisplay)
{
    // Displays rotation matrix to console without error.
    quaternion<int> q;
    // display(q);
}

TEST(QuaternionMath, ConversionFromDCM)
{
  matrix<double, 3, 3> R = random_rotation<double>();
  quaternion<double> q(R);
  ASSERT_TRUE(q == R);
}

TEST(QuaternionMath, AddRotations)
{
  // Adding rotations.
  double alpha = random_angle<double>();
  double beta  = random_angle<double>();

  quaternion<double> T1(dcm::R2<double>(alpha));
  quaternion<double> T2(dcm::R2<double>(beta));
  quaternion<double> T3 = T1 + T2;
  ASSERT_EQUAL_WITHIN_NUMERICAL_PRECISION(T3.matrix(), dcm::R2(alpha + beta));
}

TEST(QuaternionMath, SubtractRotations)
{
  // Adding rotations.
  double alpha = random_angle<double>();
  double beta = random_angle<double>();

  // Subtracting rotations.
  quaternion<double> T1(dcm::R2(alpha));
  quaternion<double> T2(dcm::R2(beta));
  quaternion<double> T3 = T1 - T2;
  ASSERT_EQUAL_WITHIN_NUMERICAL_PRECISION(T3.matrix(), dcm::R2(alpha - beta));
}

TEST(QuaternionUtility, Reverse)
{
  double alpha = random_angle<double>();

  // Reversing DCM should be the same creating one in the negative direction.
  quaternion<double> T1_pos(dcm::R1(alpha));
  quaternion<double> T2_neg(dcm::R1(-1 * alpha));
  ASSERT_TRUE(T2_neg == T1_pos.reverse());
}


#include <gtest/gtest.h>

 #include "attitude/dcm.h"
 #include "include/helper.h"

 using namespace attitude;

TEST(DCMProperties, CreateRandomRotation)
{
  // Can we create a random rotation, expressed as rotation matrix?
  matrix<double, 3, 3> T = random_rotation<double>();
}

TEST(DCMProperties, ZeroRotationIsIdentity) {
    // Zero method returns instatiates and returns identity.
    matrix<double, 3, 3> T = dcm::R1(0.) * dcm::R2(0.) * dcm::R3(0.);
    ASSERT_TRUE(dcm::ZERO<double>() == T);
}

TEST(DCMProperties, MatrixDeterminant)
{
  matrix<double, 3, 3> T = random_rotation<double>();
  ASSERT_TRUE(determinant(T) == 1.);
}

TEST(DCMProperties, MatrixOrthogonal)
{
  matrix<double, 3, 3> T = random_rotation<double>();
  ASSERT_EQUAL_WITHIN_NUMERICAL_PRECISION(T.transpose(), inverse(T));
}

TEST(DCMMath, Add180DegRotations)
{
  matrix<double, 3, 3> T1_180 = dcm::R1(deg2rad(180.));
  matrix<double, 3, 3> T2_180 = dcm::R1(deg2rad(180.));
  matrix<double, 3, 3> T3 = T1_180 * T2_180;

  ASSERT_EQUAL_WITHIN_NUMERICAL_PRECISION(T3, dcm::ZERO<double>());
}

TEST(DCMMath, AddOppositeRotations)
{
    double alpha = random_angle<double>();

    // Adding equal and opposite rotations cancel to zero.
    matrix<double, 3, 3> T1_pos = dcm::R1(alpha);
    matrix<double, 3, 3> T2_neg = dcm::R1(-1 * alpha);
    matrix<double, 3, 3> T3 = T1_pos * T2_neg;

    ASSERT_EQUAL_WITHIN_NUMERICAL_PRECISION(T3, dcm::ZERO<double>());
}

TEST(DCMUtility, Reverse)
{
    double alpha = random_angle<double>();

    // Tranposing DCM is the same creating one in the negative direction?
    matrix<double, 3, 3> T_positive = dcm::R3(alpha);
    matrix<double, 3, 3> T_negative = dcm::R3(-1 * alpha);

    ASSERT_EQUAL_WITHIN_NUMERICAL_PRECISION(T_positive.transpose(), T_negative);
}


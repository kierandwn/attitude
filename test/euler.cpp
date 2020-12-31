#include <gtest/gtest.h>

#include "include/helper.h"
#include "euler.h"

using namespace attitude;

TEST(EulerInstantiation, InstantiateEuler)
{
    // Ensures ability to create concrete class.
    euler<int> theta;
}

TEST(EulerInstantiation, DefaultIs0Rotation)
{
    // Zero method returns instatiates and returns identity.
    euler<int> theta;
    ASSERT_TRUE(theta[0] == 0 && theta[1] == 0 && theta[2] == 0);
}

TEST(EulerDebugging, STDOUTDisplay)
{
    // Displays rotation matrix to console without error.
    euler<int> theta;
    // display(theta);
}

TEST(EulerMath, ConversionFromDCM)
{
  matrix<double, 3, 3> R = dcm::R1(random_angle<double>()) * 
                           dcm::R2(random_angle<double>()) *
                           dcm::R3(random_angle<double>());

  euler<double> theta(R, randomise_order());
  ASSERT_TRUE(theta == R);
}

TEST(EulerMath, ComparisonWithDCM)
{
  // Adding rotations.
  matrix<double, 3, 3> R = dcm::R1(random_angle<double>());
  euler<double> theta(R, randomise_order());
  ASSERT_TRUE(theta == R);
}

TEST(EulerMath, AddRotations)
{
    double alpha = random_angle<double>();
    double beta  = random_angle<double>();
    
    // Adding rotations.
    euler<double> T1(dcm::R2(alpha), randomise_order());
    euler<double> T2(dcm::R2(beta),  randomise_order());
    euler<double> T3 = T1 + T2;

    ASSERT_EQUAL_WITHIN_NUMERICAL_PRECISION(T3.matrix(), dcm::R2(alpha + beta));
}

TEST(EulerMath, SubtractRotations)
{
    double alpha = random_angle<double>();
    double beta  = random_angle<double>();

    // Subtracting rotations.
    euler<double> T1(dcm::R3(alpha), randomise_order());
    euler<double> T2(dcm::R3(beta), randomise_order());

    euler<double> T3 = T1 - T2;
    ASSERT_EQUAL_WITHIN_NUMERICAL_PRECISION(T3.matrix(), dcm::R3(alpha - beta));
}

TEST(EulerUtility, Reverse)
{
    double alpha = random_angle<double>();

    // Reversing DCM should be the same creating one in the negative direction.
    euler<double> T1_pos(dcm::R1(alpha), randomise_order());
    euler<double> T2_neg(dcm::R1(-1 * alpha), randomise_order());
    ASSERT_TRUE(T2_neg == T1_pos.reverse());
}

TEST(EulerMath, DifferentialKinematicRelation) {
    // Zero rotation results in 1:1 mapping between ang. vel and euler rates
    euler<double> T(dcm::ZERO<double>(), 123);
    // display(T.dke());
    ASSERT_EQUAL_WITHIN_NUMERICAL_PRECISION(T.dke(), dcm::ZERO<double>());
}


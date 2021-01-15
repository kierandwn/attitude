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

TEST(MRPConsistency, DCMConversion) 
{
  double alpha = random_angle<double>();
  double beta = random_angle<double>();

  matrix<double, 3, 3> T1 = dcm::AXIS(2, alpha);
  matrix<double, 3, 3> T2 = dcm::AXIS(2, beta);
  matrix<double, 3, 3> T3 = T1 * T2;

  mrp<double> sigma(T3);
  sigma[0] = sigma[0]; // marks parameters as modified to indiicated DCM should be updated.

  ASSERT_EQUAL_WITHIN_NUMERICAL_PRECISION(sigma.matrix(), T3);
}

TEST(MRPMath, ConversionFromDCM)
{
  matrix<double, 3, 3> R = random_rotation<double>();
  mrp<double> q(R);
  ASSERT_TRUE(q == R);
}

TEST(MRPMath, AddRotations) {
  // Adding rotations.
  double alpha = random_angle<double>();
  double beta = random_angle<double>();

  mrp<double> T1(dcm::R3<double>(alpha));
  mrp<double> T2(dcm::R3<double>(beta));
  mrp<double> T3 = T1 + T2;
  ASSERT_EQUAL_WITHIN_NUMERICAL_PRECISION(T3.matrix(), dcm::R3(alpha + beta));
}

TEST(MRPMath, SubtractRotations) {
  // Adding rotations.
  double alpha = random_angle<double>();
  double beta = random_angle<double>();

  mrp<double> T1(dcm::R1<double>(alpha));
  mrp<double> T2(dcm::R1<double>(beta));
  mrp<double> T3 = T1 - T2;

  ASSERT_EQUAL_WITHIN_NUMERICAL_PRECISION(T3.matrix(), dcm::R1(alpha - beta));
}

TEST(MRPMath, DifferentialKinematicRelation) {
  // Zero rotation results in 1:1 mapping between ang. vel and MRP rates
  mrp<double> sigma(dcm::ZERO<double>());
  ASSERT_EQUAL_WITHIN_NUMERICAL_PRECISION(sigma.dke(), eye<double, 3>() * 0.25);
}

TEST(MRPUtility, Reverse)
{
  //double alpha = random_angle<double>();
  matrix<double, 3, 3> M = random_rotation<double>();

  // Reversing DCM should be the same creating one in the negative direction.
  mrp<double> sigma(M);
  ASSERT_EQUAL_WITHIN_NUMERICAL_PRECISION(sigma.reverse().matrix(), M.transpose());
}


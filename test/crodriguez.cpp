#include <gtest/gtest.h>

#include "include/helper.h"
#include "attitude/crp.h"

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
    ASSERT_EQUAL_WITHIN_NUMERICAL_PRECISION(T3.matrix(), dcm::R2(alpha + beta));
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
    ASSERT_EQUAL_WITHIN_NUMERICAL_PRECISION(T3.matrix(), dcm::R2<double>(alpha - beta));
}

TEST(CRPUtility, Reverse)
{
  double alpha = random_angle<double>();
  matrix<double, 3, 3> M = random_rotation<double>();

  // Reversing DCM should be the same creating one in the negative direction.
  crp<double> q(M);

  /*display(q.reverse().matrix());
  display(M.transpose());*/

  ASSERT_EQUAL_WITHIN_NUMERICAL_PRECISION(q.reverse().matrix(), M.transpose());
}

TEST(CRPMath, DifferentialKinematicRelation) {
  // Zero rotation results in 1:1 mapping between ang. vel and CRP rates
  crp<double> q(dcm::ZERO<double>());
  ASSERT_EQUAL_WITHIN_NUMERICAL_PRECISION(q.dke(), eye<double, 3>() * 0.5);
}


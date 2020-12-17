#include <gtest/gtest.h>

#include <cstdlib>
#include <iostream>
#include <ctime>

#include "dcm.h"
#include "matrix.h"
#include "include/helpers.h"

using namespace attitude;


TEST(DCMInstantiation, InstantiateDCM)
{
    // Ensures ability to create concrete class.
    dcm_<int> T = ZERO<int>();
}

TEST(DCMInstantiation, DefaultIsIdentity)
{
    // Zero method returns instatiates and returns identity.
    dcm_<double> T = ZERO<double>();
    linalg::square<double> I( 3, 1., 0., 0., 0., 1., 0., 0., 0., 1.);
    ASSERT_TRUE(I == T.matrix()); // copy const. used? why?
}

TEST(DCMDebugging, IndexOperator)
{
    // Row (list of type T) is accessed by first index.
    dcm_<int> T = ZERO<int>();
    ASSERT_TRUE(T[0][0] == 1 && T[0][1] == 0);
}

TEST(DCMDebugging, STDOUTDisplay)
{
    // Displays rotation matrix to console without error.
    dcm_<int> T = ZERO<int>();
    display(T);
}

TEST(DCMProperties, MatrixDeterminant)
{
    dcm_<double> T = random_rotation<double>();
    linalg::square<double> M = T.matrix();
    ASSERT_TRUE(std::abs(M.determinant() - 1.) < 1E-06);
}

TEST(DCMProperties, MatrixOrthogonal)
{
    dcm_<double> T = random_rotation<double>();
    linalg::square<double> M = T.matrix();
    ASSERT_EQUAL_WITHIN_NUMERICAL_PRECISION(M.transpose(), M.inverse());
}

TEST(DCMMath, AddRotations)
{
    double alpha = random_angle<double>();
    double beta  = random_angle<double>();

    dcm_<double> T1 = R1<double>(alpha);
    dcm_<double> T2 = R1<double>(beta);
    
    // Adding rotations.
    dcm_<double> T3 = T1 + T2;
    ASSERT_TRUE(T3 == R1<double>(alpha + beta));
}

TEST(DCMMath, Add180DegRotations)
{
    dcm_<double> T1_180 = R1(deg2rad(180.));
    dcm_<double> T2_180 = R1(deg2rad(180.));
    dcm_<double> T3 = T1_180 + T2_180;
    ASSERT_TRUE(T3 == ZERO<double>());
}

TEST(DCMMath, AddOppositeRotations)
{
    double alpha = random_angle<double>();

    // Adding equal and opposite rotations cancel to zero.
    dcm_<double> T1_pos = R1(alpha);
    dcm_<double> T2_neg = R1(-1 * alpha);
    dcm_<double> T3 = T1_pos + T2_neg;
    ASSERT_TRUE(T3 == ZERO<double>());
}

TEST(Math, SubtractRotations)
{
    double alpha = random_angle<double>();
    double beta  = random_angle<double>();

    dcm_<double> T1 = R2<double>(alpha);
    dcm_<double> T2 = R2<double>(beta);
    
    // Subtract rotation.
    dcm_<double> T3 = T1 - T2;
    ASSERT_TRUE(T3 == R2<double>(alpha - beta));
}

TEST(DCMUtility, Reverse)
{
    double alpha = random_angle<double>();

    // Reversing DCM is the same creating one in the negative direction.
    dcm_<double> T_pos = R3(alpha);
    dcm_<double> T_neg = R3(-1 * alpha);
    ASSERT_TRUE(T_neg == T_pos.reverse());
}


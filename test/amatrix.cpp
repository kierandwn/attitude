#include <gtest/gtest.h>

#include "matrix.h"

using namespace attitude::linalg;

TEST(MatrixInstantiation, InstantiateMatrix)
{
    // Ensures ability to create concrete class.
    matrix<int> M(3, 3);
    square<double> I(3, 1., 0., 0., 0., 1., 0., 0., 0., 1.);
    
    ASSERT_TRUE(I[0][0] == 1.);
    ASSERT_TRUE(I[1][1] == 1.);
}

TEST(MatrixInstantiation, TypeConversion)
{
    // Zero method returns instantiates and returns identity.
    square<double> M(3, 1., 0., 0., 0., 1, 0., 0., 0., 1.);
    ASSERT_TRUE(M[1][1] == 1.);
}

TEST(MatrixOperations, CompareMatrixWithSquare)
{
    // Ensures ability to compare base and derived classes.
    square<int> M(3, 1, 0, 0, 0, 1, 0, 0, 0, 1);
    matrix<int> I(3, 3, 1, 0, 0, 0, 1, 0, 0, 0, 1);

    ASSERT_TRUE(M == I);
    ASSERT_TRUE(I == M);
}

TEST(MatrixOperations, DeterminantMatrix)
{
    // Ensures ability to create concrete class.
    square<int> M(3, 1, 0, 0, 0, 1, 0, 0, 0, 1);
    ASSERT_TRUE(M.determinant() == 1);
}

TEST(MatrixOperations, InverseMatrix)
{
    // Ensures ability to create concrete class.
    square<int> M(3, 1, 0, 0, 0, 1, 0, 0, 0, 1);
    ASSERT_TRUE(M.inverse() == M);
}

#include <gtest/gtest.h>

#include "attitude/matrix.h"

using namespace attitude;

TEST(MatrixInstantiation, DefautConstructor)
{
  // Defaut constructor is valid?
  matrix<double, 3, 3> a();
}

TEST(MatrixInstantiation, SingleValInitList)
{
  matrix<int, 2, 2> a{1};
}

TEST(MatrixInstantiation, FullMatrixInitList) 
{
  matrix<double, 3, 3> a{1., 0., 0.,
                         0., 1,  0.,
                         0., 0., 1.};
}

TEST(MatrixInstantiation, ConstructFromArray) 
{
  int arr[9]{1, 0, 0, 0, 1, 0, 0, 0, 1};

  matrix<int, 3, 3> a(arr);
}

TEST(MatrixInstantiation, CopyConstruction) {
  int arr[9]{1, 0, 0, 0, 1, 0, 0, 0, 1};
  matrix<int, 3, 3> a(arr);

  matrix<int, 3, 3> b(a);
}

//TEST(MatrixConversion, TypeConversionOnArrayConstruction) {
//  int arr[9]{1, 0, 0, 0, 1, 0, 0, 0, 1};
//
//  matrix<double, 3, 3> a(arr);
//}

TEST(MatrixConversion, TypeConversionOnCopy) {
  int arr[9]{1, 0, 0, 0, 1, 0, 0, 0, 1};
  matrix<int, 3, 3> a(arr);

  matrix<double, 3, 3> b(a);
}

TEST(MatrixConversion, DimensionSwap) {
  // Can another matrix of different dimensions be created from the first.
  // NxM must match.
  matrix<int, 3, 4> a{1, 0, 0, 0,
                      0, 1, 0, 0,
                      0, 0, 1, 0};

  matrix<int, 4, 3> b(a);
}

TEST(IdentityProperties, Orthogonal) {
  // Can another matrix of different dimensions be created from the first.
  // NxM must match.
  matrix<double, 3, 3> I{1, 0, 0, 0, 1, 0, 0, 0, 1};

  ASSERT_TRUE(I == I.transpose());
}

TEST(IdentityProperties, Determinant) {
  // Is determinant of identity == 1?
  matrix<double, 3, 3> I{1., 0., 0., 
                         0., 1., 0., 
                         0., 0., 1.};

  matrix<double, 2, 2> reduced = remove_rowi_colj(0, 0, I);
  ASSERT_TRUE(determinant(I) == 1.);
}

TEST(IdentityProperties, Inverse) {
  // Inverse of identity is identity?
  matrix<double, 3, 3> I{1., 0., 0., 
                         0., 1., 0., 
                         0., 0., 1.};

  ASSERT_TRUE(inverse(I) == I);
}



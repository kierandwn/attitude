#include <gtest/gtest.h>


// #include "euler.h"
#include "matrix.h"

using namespace attitude;

TEST(VectorInstantiation, DefautConstructor) {
  // Defaut constructor is valid?
  vector<double, 3> a;
}

TEST(VectorInstantiation, InitListVector) {
  vector<double, 3> a{0, 0, 0};
}

TEST(VectorMath, MultVector) {
  vector<double, 3> a{1, 2, 3};
  vector<double, 3> b = eye<double, 3>() * a;
}

//TEST(Test, DescriptionTest) {
//  euler<double> a;
//}


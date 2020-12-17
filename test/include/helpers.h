#ifndef ATT_TEST_HELPERS_H
#define ATT_TEST_HELPERS_H

#include <gtest/gtest.h>

#include "dcm.h"

const double PI = 3.14159;

int random_not(int);
int randomise_order();

template<typename T>
T random_angle() { return T(4. * PI * std::rand() / RAND_MAX - 2. * PI); };

template<typename T>
attitude::dcm_<T> random_rotation() {
	T alpha = random_angle<T>();
	T beta  = random_angle<T>();
	T gamma = random_angle<T>();
	return attitude::R1(alpha) + attitude::R2(beta) + attitude::R3(gamma);
};

template <typename T>
T deg2rad(T deg) { return deg * PI / 180; };

template <typename T>
T rad2deg(T rad) { return rad * 180 / PI; };

template <typename T>
void ASSERT_EQUAL_WITHIN_NUMERICAL_PRECISION(attitude::linalg::matrix<T> M1, attitude::linalg::matrix<T> M2)
{
	attitude::linalg::matrix<T> d = M2 - M1;
	ASSERT_TRUE(d == 0);
}

#endif // ATT_TEST_HELPERS_H

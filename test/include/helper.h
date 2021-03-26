#ifndef ATT_TEST_HELPERS_H
#define ATT_TEST_HELPERS_H

#include <gtest/gtest.h>

#include "attitude/dcm.h"

using namespace attitude;

const double kPi = 3.14159;

int random_not(int);
int randomise_order();

template<typename Tp>
Tp random_angle() { return Tp(4. * kPi * std::rand() / RAND_MAX - 2. * kPi); };

template<typename Tp, size_t rows, size_t cols>
matrix<Tp, rows, cols> abs(matrix<Tp, rows, cols> M) {
  for (int i = 0; i < rows; ++i)
    for (int j = 0; j < cols; ++j) M[i][j] = std::abs(M[i][j]);

	return M;
}

template<typename Tp>
matrix<Tp, 3, 3> random_rotation() {
	Tp alpha = random_angle<Tp>();
	Tp beta  = random_angle<Tp>();
	Tp gamma = random_angle<Tp>();
	return dcm::R1(alpha) * dcm::R2(beta) * dcm::R3(gamma);
};

template <typename Tp>
Tp deg2rad(Tp deg) { return deg * kPi / 180.0l; };

template <typename Tp>
Tp rad2deg(Tp rad) { return rad * 180 / kPi; };

template <typename Tp, size_t rows, size_t cols>
void ASSERT_EQUAL_WITHIN_NUMERICAL_PRECISION(matrix<Tp, rows, cols> M1,
                                             matrix<Tp, rows, cols> M2)
{
	matrix<Tp, rows, cols> d = M2 - M1;
	ASSERT_TRUE(d <= 1E-05); // TODO: improve precision: accuracy 1e-05 with double?
}

#endif // ATT_TEST_HELPERS_H

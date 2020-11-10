#ifndef ATT_COMMON_H
#define ATT_COMMON_H

#include "matrix9.h"

namespace attitude
{
	namespace math {
		const double PI = 3.14159;

		template <typename T>
		T deg2rad(T deg)
		{
			return deg * PI / 180;
		};

		template <typename T>
		T rad2deg(T rad)
		{
			return rad * 180 / PI;
		};
	}
}

#endif // ATT_COMMON_H
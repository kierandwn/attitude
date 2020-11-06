#ifndef ATT_COMMON_H
#define ATT_COMMON_H

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

		template <typename T>
		std::array<T, 9> matrix_multiply(std::array<T, 9> M1, std::array<T, 9> M2)
		{
			return {
				// Row 1
				M1[0] * M2[0] + M1[1] * M2[3] + M1[2] * M2[6], // M3[0,0]
				M1[0] * M2[1] + M1[1] * M2[4] + M1[2] * M2[7], // M3[0,1]
				M1[0] * M2[2] + M1[1] * M2[5] + M1[2] * M2[8], // M3[0,2]
				// Row 2
				M1[3] * M2[0] + M1[4] * M2[3] + M1[5] * M2[6], // M3[1,0]
				M1[3] * M2[1] + M1[4] * M2[4] + M1[5] * M2[7], // M3[1,1]
				M1[3] * M2[2] + M1[4] * M2[5] + M1[5] * M2[8], // M3[1,2]
				// Row 3
				M1[6] * M2[0] + M1[7] * M2[3] + M1[8] * M2[6], // M3[2,0] 
				M1[6] * M2[1] + M1[7] * M2[4] + M1[8] * M2[7], // M3[2,1]
				M1[6] * M2[2] + M1[7] * M2[5] + M1[8] * M2[8]  // M3[2,2]
			};
		}

		template <typename T>
		std::array<T, 9> matrix_transpose(std::array<T, 9> M)
		{
			std::array<T, 9> transposed;

			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					transposed[i * 3 + j] = M[j * 3 + i];
				}
			}
			return transposed;
		}

		template <typename T>
		bool matrix_is_value(std::array<T, 9> M, T val)
		{
			bool result = true;

			for (int i = 0; i < M.size(); i++)
			{
				if (std::abs(M[i] - val) > 1e-05)
					result = false; break;
			}

			return result;
		}

		template <typename T>
		std::array<T, 9> array_diff(std::array<T, 9> A1, std::array<T, 9> A2)
		{
			std::array<T, 9> result;

			for (int i = 0; i < 9; i++)
			{
				result[i] = A1[i] - A2[i];
			}
			return result;
		}
	}
}



#endif // ATT_COMMON_H
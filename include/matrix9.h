#ifndef ATT_MATRIX_H
#define ATT_MATRIX_H

#include <array>

namespace attitude
{
	namespace math
	{
		template <typename T>
		class matrix9
		{
		private:
			T _items[9];

		public:
			matrix9() : _items{ 0 }
			{}

			matrix9(T _1, T _2, T _3,
					T _4, T _5, T _6,
					T _7, T _8, T _9 ) : _items{ _1, _2, _3, _4, _5, _6, _7, _8, _9 }
			{}

			int sizeofthis() const { return 9; }

			std::array<T, 3> operator[] (int i) 
			{ 
				return { _items[i * 3], _items[i * 3 + 1], _items[i * 3 + 2] };
			}

			const std::array<T, 3> operator[] (int i) const 
			{ 
				return { _items[i * 3], _items[i * 3 + 1], _items[i * 3 + 2] };
			}

			matrix9<T> operator- (matrix9<T> M)
			{
				matrix9<T> result;

				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						result[i][j] = operator[](i)[j] - M[i][j];
					}
				}
				return result;
			}

			matrix9 operator* (matrix9<T> M)
			{
				return matrix9({
					// Row 1
					operator[](0)[0] * M[0][0] + operator[](0)[1] * M[1][0] + operator[](0)[2] * M[2][0], // M3[0,0]
					operator[](0)[0] * M[0][1] + operator[](0)[1] * M[1][1] + operator[](0)[2] * M[2][1], // M3[0,1]
					operator[](0)[0] * M[0][2] + operator[](0)[1] * M[1][2] + operator[](0)[2] * M[2][2], // M3[0,2]
					// Row 2
					operator[](1)[0] * M[0][0] + operator[](1)[1] * M[1][0] + operator[](1)[2] * M[2][0], // M3[1,0]
					operator[](1)[0] * M[0][1] + operator[](1)[1] * M[1][1] + operator[](1)[2] * M[2][1], // M3[1,1]
					operator[](1)[0] * M[0][2] + operator[](1)[1] * M[1][2] + operator[](1)[2] * M[2][2], // M3[1,2]
					// Row 3
					operator[](2)[0] * M[0][0] + operator[](2)[1] * M[1][0] + operator[](2)[2] * M[2][0], // M3[2,0] 
					operator[](2)[0] * M[0][1] + operator[](2)[1] * M[1][1] + operator[](2)[2] * M[2][1], // M3[2,1]
					operator[](2)[0] * M[0][2] + operator[](2)[1] * M[1][2] + operator[](2)[2] * M[2][2]  // M3[2,2]
				});
			}

			bool operator==(matrix9<T> M)
			{
				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						if (operator[](i)[j] != M[i][j])
							return false;
					}
				}

				return true;
			}

			bool operator==(T val)
			{
				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						if (std::abs(operator[](i)[j] - val) > 1e-05)
							return false;
					}
				}

				return true;
			}

			matrix9<T> transpose()
			{
				matrix9<T> transposed;

				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						transposed[i * 3 + j] = operator[](j * 3 + i);
					}
				}
				return transposed;
			}
		};

	}
}
#endif // ATT_MATRIX_H
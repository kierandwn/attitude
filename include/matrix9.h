#ifndef ATT_MATRIX_H
#define ATT_MATRIX_H

namespace attitude
{
	namespace math
	{
		template <typename T>
		class matrix9
		{
		private:
			T _items[9];
			T& _get(int i, int j) { return operator[](i)[j]; }
			
			T _cofactor(int i, int j)
			{
				int* horz_dims = new int[2];
				int* vert_dims = new int[2];

				horz_dims = remove_dim(i, horz_dims);
				vert_dims = remove_dim(j, vert_dims);

				T a = _get(horz_dims[0], vert_dims[0]);
				T b = _get(horz_dims[1], vert_dims[0]);
				T c = _get(horz_dims[0], vert_dims[1]);
				T d = _get(horz_dims[1], vert_dims[1]);

				// determinant of cofactor matrix
				return a * d - b * c;
			}

		public:
			matrix9() : _items{ 0 }{}

			matrix9(T _1, T _2, T _3,
					T _4, T _5, T _6,
					T _7, T _8, T _9 ) : _items{ _1, _2, _3, _4, _5, _6, _7, _8, _9 }{}

			matrix9<T> transpose()
			{
				matrix9<T> transposed;

				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						transposed[i][j] = operator[](j)[i];
					}
				}
				return transposed;
			}

			matrix9<T> trace() 
			{ 
				return operator[](0)[0] **2 + operator[](1)[1] **2 + operator[](2)[2] **2; 
			}

			T det() 
			{
				return operator[](0)[0] * (operator[](1)[1] * operator[](2)[2] - operator[](1)[2] * operator[](2)[1]) -
					   operator[](0)[1] * (operator[](1)[0] * operator[](2)[2] - operator[](1)[2] * operator[](2)[0]) +
					   operator[](0)[2] * (operator[](1)[0] * operator[](2)[1] - operator[](1)[1] * operator[](2)[0]);
			}

			matrix9<T> inverse() {
				// Form matrix of cofactors
				matrix9<T> cofactors = *this;
				int scalar;

				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++) {
						scalar = ((i * 3 + j) % 2 == 0) ? 1 : -1;
						cofactors[i][j] = scalar * _cofactor(i, j);
					}
				}
				
				// Compute inverse
				return cofactors.transpose() * (T(1.) / det());
			}

			// Elementwise operations.
			matrix9<T> operator+ (matrix9<T> M)
			{
				matrix9<T> result;

				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						result[i][j] = operator[](i)[j] + M[i][j];
					}
				}
				return result;
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

			matrix9<T> operator+ (T val)
			{
				matrix9<T> result;

				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						result[i][j] = operator[](i)[j] + val;
					}
				}
				return result;
			}

			matrix9<T> operator- (T val)
			{
				matrix9<T> result;

				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						result[i][j] = operator[](i)[j] - val;
					}
				}
				return result;
			}

			matrix9<T> operator* (T val)
			{
				matrix9<T> result;

				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						result[i][j] = operator[](i)[j] * val;
					}
				}
				return result;
			}

			matrix9<T> operator/ (T val)
			{
				matrix9<T> result;

				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						result[i][j] = operator[](i)[j] / val;
					}
				}
				return result;
			}

			matrix9<T>* operator+= (matrix9<T> M)
			{
				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						operator[](i)[j] += M[i][j];
					}
				}
				return this;
			}

			matrix9<T>* operator-= (matrix9<T> M)
			{
				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						operator[](i)[j] -= M[i][j];
					}
				}
				return this;
			}

			matrix9<T>* operator+= (T val)
			{
				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						operator[](i)[j] += val;
					}
				}
				return this;
			}

			matrix9<T>* operator*= (T val)
			{
				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						operator[](i)[j] *= val;
					}
				}
				return this;
			}

			matrix9<T>* operator/= (T val)
			{
				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						operator[](i)[j] /= val;
					}
				}
				return this;
			}

			// Matrix multiplication
			matrix9<T> operator* (matrix9<T> M)
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

			matrix9<T>* operator*= (matrix9<T> M)
			{
				// Row 1
				operator[](0)[0] = operator[](0)[0] * M[0][0] + operator[](0)[1] * M[1][0] + operator[](0)[2] * M[2][0]; // M3[0,0]
				operator[](0)[1] = operator[](0)[0] * M[0][1] + operator[](0)[1] * M[1][1] + operator[](0)[2] * M[2][1]; // M3[0,1]
				operator[](0)[2] = operator[](0)[0] * M[0][2] + operator[](0)[1] * M[1][2] + operator[](0)[2] * M[2][2]; // M3[0,2]
				// Row 2
				operator[](1)[0] = operator[](1)[0] * M[0][0] + operator[](1)[1] * M[1][0] + operator[](1)[2] * M[2][0]; // M3[1,0]
				operator[](1)[1] = operator[](1)[0] * M[0][1] + operator[](1)[1] * M[1][1] + operator[](1)[2] * M[2][1]; // M3[1,1]
				operator[](1)[2] = operator[](1)[0] * M[0][2] + operator[](1)[1] * M[1][2] + operator[](1)[2] * M[2][2]; // M3[1,2]
				// Row 3
				operator[](2)[0] = operator[](2)[0] * M[0][0] + operator[](2)[1] * M[1][0] + operator[](2)[2] * M[2][0]; // M3[2,0] 
				operator[](2)[1] = operator[](2)[0] * M[0][1] + operator[](2)[1] * M[1][1] + operator[](2)[2] * M[2][1]; // M3[2,1]
				operator[](2)[2] = operator[](2)[0] * M[0][2] + operator[](2)[1] * M[1][2] + operator[](2)[2] * M[2][2]; // M3[2,2]
				return this;
			}

			// Division is inversion followed by multiplication 
			matrix9<T> operator/ (matrix9<T> M)
			{
				return M.inverse() * this;
			}

			matrix9<T>* operator/= (matrix9<T> M)
			{
				return M.inverse() *= this;
			}

			// Comparison
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

			// Indexing
			T* operator[] (int i) { return &_items[i * 3]; }
		};

		inline int* remove_dim(int k, int* dims)
		{
			if (k == 0)
			{
				dims[0] = 1;
				dims[1] = 2;
			}
			else if (k == 1)
			{
				dims[0] = 0;
				dims[1] = 2;
			}
			else if (k == 2)
			{
				dims[0] = 0;
				dims[1] = 1;
			}
			return dims;
		}

		// 3x3 Identity matrix
		template <typename T>
		matrix9<T> EYE3{ T(1.), T(0.), T(0.), T(0.), T(1.), T(0.), T(0.), T(0.), T(1.) };
	}
}

template<typename T>
void display(attitude::math::matrix9<T> M)
{
	for (int i = 0; i < 3; i++)
	{
		T* row = M[i];

		std::cout
			<< "[ "
			<< row[0] << ","
			<< row[1] << ","
			<< row[2]
			<< " ]"
			<< std::endl;
	}
}
#endif // ATT_MATRIX_H
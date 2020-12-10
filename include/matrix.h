#ifndef ATT_MATRIX_H
#define ATT_MATRIX_H

#include <cstdarg>

namespace attitude
{
	namespace math
	{
		inline void remove_dim(int k, const int length, int* dims)
		{
			int j = 0;

			for (int i = 0; i < length; i++)
				if (k != i) { dims[j] = i; j += 1; }
		};

		template<typename T, int N>
		class matrix // square matrix
		{
		protected:
			const int _size = N;
			T _items[N * N] = { 0 };

			T _get(int i, int j) { return operator[](i)[j]; }
			void _set(int i, int j, T val) { operator[](i)[j] = val; }
			
			void _set(T arr[N * N]) { for (int i = 0; i < N * N; i++) { _items[i] = arr[i]; } }

			matrix<T, N> _cofactor()
			{
				if (_size > 2) {
					T* cof = new T[_size * _size];
					int scalar;

					for (int i = 0; i < _size; i++) {
						int* rows = new int[_size - 1];
						remove_dim(i, 3, rows);

						for (int j = 0; j < _size; j++) {
							int* cols = new int[_size - 1];
							remove_dim(j, 3, cols);

							scalar = ((i * 3 + j) % 2 == 0) ? 1 : -1;
							cof[i * _size + j] = scalar * determinant(_size - 1, rows, cols);
							delete[] cols;
						}
						delete[] rows;
					}
					matrix<T, N> _cof(cof);
					delete[] cof;
					return _cof;
				}
			}

		public:
			matrix<T, N> cofactors() { return _cofactor(); }

			matrix() : _items{ 0 } {}

			matrix(int sz, ...) : _size(sz)
			{
				va_list args;
				va_start(args, sz);

				for (int k = 0; k < _size * _size; k++)
					_items[k] = va_arg(args, T);

				va_end(args);
			}

			matrix(T* arr) {
				for (int i = 0; i < _size; i++)
					for (int j = 0; j < _size; j++)
						_set(i, j, arr[i * _size + j]);
			}

			int size() { return _size; }

			T trace() { 
				T trace = 0;
				for (int k = 0; k < _size; k++)
					trace += this[k][k] ** 2;

				return T;
			}

			matrix<T, N> transpose()
			{
				matrix<T, N> transposed(_size);

				for (int i = 0; i < _size; i++)
					for (int j = 0; j < _size; j++)
						transposed[i][j] = _get(j, i);

				return transposed;
			}

			T determinant(int size, int* rows, int* cols)
			{
				if (size == 2) {
					T a = _get(rows[0], cols[0]);
					T b = _get(rows[0], cols[1]);
					T c = _get(rows[1], cols[0]);
					T d = _get(rows[1], cols[1]);
					return  a * d - b * c;
				}
				else
				{
					int* rows = new int[size - 1];
					remove_dim(0, size, rows); // TODO also need to keep track of previously removed dims in the case of N > 3.

					T det = 0.;
					int scalar;

					for (int i = 0; i < size; i++)
					{	
						int* cols = new int[size - 1];
						remove_dim(i, size, cols);

						scalar = ((i * 3) % 2 == 0) ? 1 : -1;
						det += scalar * _get(0, i) * determinant(size - 1, rows, cols);
						delete cols;
					}
					delete rows;
					return det;
				}
			}

			T determinant()
			{
				int* rows = new int[_size];
				int* cols = new int[_size];

				for (int i = 0; i < _size; i++) { rows[i] = i; cols[i] = i; }
				
				T det = determinant(_size, rows, cols);
				delete rows, cols;
				return det;
			}

			matrix<T, N> inverse() {
				return _cofactor().transpose() * (T(1.) / determinant());
			}

			// Elementwise operations.
			matrix<T, N> operator+ (matrix<T, N> M)
			{
				matrix<T, N> result;

				for (int i = 0; i < _size; i++)
					for (int j = 0; j < _size; j++)
						result[i][j] = _get(i, j) + M[i][j];

				return result;
			}

			matrix<T, N> operator- (matrix<T, N> M)
			{
				matrix<T, N> result;

				for (int i = 0; i < _size; i++)
					for (int j = 0; j < _size; j++)
						result[i][j] = _get(i, j) - M[i][j];

				return result;
			}

			matrix<T, N> operator+ (T val)
			{
				matrix<T, N> result(_size);

				for (int i = 0; i < _size; i++)
					for (int j = 0; j < _size; j++)
						result[i][j] = operator[](i)[j] + val;

				return result;
			}

			matrix<T, N> operator- (T val)
			{
				matrix<T, N> result(_size);

				for (int i = 0; i < _size; i++)
					for (int j = 0; j < _size; j++)
						result[i][j] = operator[](i)[j] - val;
					
				return result;
			}

			// Elementwise multiplication
			matrix<T, N> operator* (T val)
			{
				matrix<T, N> result(_size);

				for (int i = 0; i < _size; i++)
					for (int j = 0; j < _size; j++)
						result[i][j] = operator[](i)[j] * val;
					
				return result;
			}

			matrix<T, N> operator/ (T val)
			{
				matrix<T, N> result(_size);

				for (int i = 0; i < _size; i++)
					for (int j = 0; j < _size; j++)
						result[i][j] = operator[](i)[j] / val;
					
				return result;
			}

			// Matrix multiplication
			matrix<T, N> operator* (matrix<T, N> M)
			{
				matrix<T, N> result(_size);

				for (int i = 0; i < _size; i++) {
					for (int j = 0; j < _size; j++) {
						result[i][j] = T(0.);
						for (int k = 0; k < _size; k++)
							result[i][j] += _get(i, k) * M[k][j];
					}
				}
				return result;
			}

			matrix<T, N>* operator*= (matrix<T, N> M)
			{
				T* items = new T[N * N]{ T(0.) };

				for (int i = 0; i < _size; i++) {
					for (int j = 0; j < _size; j++) {
						for (int k = 0; k < _size; k++)
							items[i * _size + j] += _get(i, k) * M[k][j];
					}
				}
				
				_set(items);
				delete items;
				
				return this;
			}

		//	// Division is inversion followed by multiplication 
		//	matrix9<T> operator/ (matrix9<T> M)
		//	{
		//		return M.inverse() * this;
		//	}

		//	matrix9<T>* operator/= (matrix9<T> M)
		//	{
		//		return M.inverse() *= this;
		//	}

			// Comparison
			bool operator==(matrix<T, N> M)
			{
				for (int i = 0; i < _size; i++)
					for (int j = 0; j < _size; j++)
						if (_get(i, j) != M[i][j]) { return false; }

				return true;
			}

			bool operator==(T val)
			{
				for (int i = 0; i < _size; i++)
					for (int j = 0; j < _size; j++)
						if (std::abs(_get(i, j) - val) > 1e-05) { return false; }

				return true;
			}

			// Indexing
			T* operator[] (int i) { return &_items[i * _size]; }
		};
	}
}

template<typename T, int n>
void display(attitude::math::matrix<T, n> M)
{
	for (int i = 0; i < n; i++)
	{
		std::cout << "[ ";
		for (int j = 0; j < n; j++)
			std::cout << M[i][j] << ",";

		std::cout << " ]" << std::endl;
	}
}
#endif // ATT_MATRIX_H
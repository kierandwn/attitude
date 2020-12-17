#ifndef ATT_MATRIX_H
#define ATT_MATRIX_H

#include <cstdarg>

namespace attitude
{
	namespace linalg
	{
		inline void remove_dim(int k, const int length, int* dims, int* reduced_dims)
		{
			int j = 0;

			for (int i = 0; i < length; i++)
				if (k != i) { reduced_dims[j] = dims[i]; j += 1; }
		};

		template<typename T>
		class matrix
		{
		protected:
			const int _shape[2];
			const int _size;

			const int _rows() { return _shape[0]; }
			const int _cols() { return _shape[1]; }

			T * _items; // enforce _size limit in memory?
			bool _allocated = false;

			T _get(int i, int j) { return _items[i * _rows() + j]; }
			void _set(int i, int j, T val) { operator[](i)[j] = val; }
			
			void _set(T* arr) { for (int i = 0; i <  _size; i++) { _items[i] = arr[i]; } }

			void _allocate() { if ( !_allocated ) { _items = new T[_size]; _allocated = true; } }
			void _deallocate() { if ( _allocated ) { delete[] _items;  _allocated = false; } }

		public:
			//template<typename... elements>
			//matrix(const int rows, const int cols, elements... items) : _shape{ rows, cols },
			//															_size(rows* cols)
			//{
			//	const int n = sizeof...(items);
			//	// TODO : Raise exception if rows * cols != n

			//		//T* items = new T[shape[0] * shape[1]]{ T(0.) }
			//	_items = { items... };
			//	//elements = [ items... ];
			//}

			matrix(const int rows, const int cols, ...) : _shape{ rows, cols },
														  _size(rows * cols)
			{
				_allocate();

				va_list args;
				va_start(args, cols);

				for (int k = 0; k < _size; k++)
					_items[k] = va_arg(args, T); // doesn't seem to correctly convert to type

				va_end(args);
			}

			matrix(const int rows, const int cols, T val) : _shape{ rows, cols },
															_size(rows * cols)
			{
				_allocate();

				for (int i = 0; i < rows; i++)
					for (int j = 0; j < cols; j++)
						_items[i * rows + j] = val;
			}

			matrix(const int rows, const int cols, T* arr) : _shape{ rows, cols },
															 _size(rows * cols) 
			{
				_allocate();

				// TODO : Raise exception if rows * cols != len(arr)
				for (int i = 0; i < rows; i++)
					for (int j = 0; j < cols; j++)
						_set(i, j, arr[i * rows + j]);
			}

			// Copy constructor
			matrix(matrix<T>& M) : _shape{ M.shape()[0], M.shape()[1] },
								   _size(M.size())
			{
				_allocate();

				// TODO : Raise exception if rows * cols != len(arr)
				for (int i = 0; i < _rows(); i++)
					for (int j = 0; j < _cols(); j++)
						_set(i, j, M[i][j]);
			}

			~matrix() { _deallocate(); }

			int size() { return _size; }
			const int * shape() { return _shape; }

			matrix<T> transpose()
			{
				matrix<T> transposed(_cols(), _rows());

				for (int i = 0; i < _rows(); i++)
					for (int j = 0; j < _cols(); j++)
						transposed[i][j] = _get(j, i);

				return transposed;
			}

			// Elementwise operations.
			matrix<T> operator+ (matrix<T> M)
			{
				matrix<T> result(_rows(), _cols());

				for (int i = 0; i < _rows(); i++)
					for (int j = 0; j < _cols(); j++)
						result[i][j] = _get(i, j) + M[i][j];

				return result;
			}

			matrix<T> operator- (matrix<T> M)
			{
				matrix<T> result(_rows(), _cols());

				for (int i = 0; i < _rows(); i++)
					for (int j = 0; j < _cols(); j++)
						result[i][j] = _get(i, j) - M[i][j];

				return result;
			}

			matrix<T> operator+ (T val)
			{
				matrix<T> result(_rows(), _cols());

				for (int i = 0; i < _rows(); i++)
					for (int j = 0; j < _cols(); j++)
						result[i][j] = operator[](i)[j] + val;

				return result;
			}

			matrix<T> operator- (T val)
			{
				matrix<T> result(_rows(), _cols());

				for (int i = 0; i < _rows(); i++)
					for (int j = 0; j < _cols(); j++)
						result[i][j] = operator[](i)[j] - val;
					
				return result;
			}

			// Elementwise multiplication
			matrix<T> operator* (T val)
			{
				matrix<T> result(_rows(), _cols());

				for (int i = 0; i < _rows(); i++)
					for (int j = 0; j < _cols(); j++)
						result[i][j] = operator[](i)[j] * val;
					
				return result;
			}

			matrix<T> operator/ (T val)
			{
				matrix<T> result(_rows(), _cols());

				for (int i = 0; i < _rows(); i++)
					for (int j = 0; j < _cols(); j++)
						result[i][j] = operator[](i)[j] / val;
					
				return result;
			}

			// Matrix multiplication
			matrix<T> operator* (matrix<T> M)
			{
				// TODO: raise exception if inner dimensions don't match
				const int shape[2]{ _rows(), M.shape()[1] };

				matrix<T> result(shape[0], shape[1]);

				for (int i = 0; i < shape[0]; i++) {
					for (int j = 0; j < shape[1]; j++) {
						result[i][j] = T(0.);
						for (int k = 0; k < _cols(); k++)
							result[i][j] += _get(i, k) * M[k][j];
					}
				}
				return result;
			}

			matrix<T>* operator*= (matrix<T> M)
			{
				// TODO: raise exception if inner dimensions don't match
				const int shape[2]{ _rows(), M.shape()[1] };

				// create new items array while
				T* items = new T[shape[0] * shape[1]]{ T(0.) };

				for (int i = 0; i < shape[0]; i++) {
					for (int j = 0; j < shape[1]; j++) {
						for (int k = 0; k < _cols(); k++)
							items[i * shape[0] + j] += _get(i, k) * M[k][j];
					}
				}
				
				_set(items);
				delete items;
				
				return this;
			}

			// Comparison
			bool operator==(matrix<T> M)
			{
				// TODO: raise exception if dimensions don't match (or simply return false?)
				for (int i = 0; i < _rows(); i++)
					for (int j = 0; j < _cols(); j++)
						if (_get(i, j) != M[i][j]) { return false; }

				return true;
			}

			bool operator==(T val)
			{
				for (int i = 0; i < _rows(); i++)
					for (int j = 0; j < _cols(); j++)
						if (std::abs(_get(i, j) - val) > 1e-05) { return false; }

				return true;
			}

			// Indexing
			T* operator[] (int i) { return &_items[i * _rows()]; }
		};

		 template<typename T>
		 class square : virtual public matrix<T> {
		 public:
			 // template<typename... elements>
			 square(const int n, ...) : matrix(n, n)
			 {
				 va_list args;
				 va_start(args, n);

				 for (int k = 0; k < _size; k++)
					 _items[k] = va_arg(args, T);

				 va_end(args);
			 }

			 square(const int n, T* arr) : matrix(n, n, arr) {}

			 // convert base-class matrix into square
			 square(matrix<T> M) : matrix( M.shape()[0], M.shape()[1] )
			 {
				 _allocate();

				 // TODO : Raise exception if rows * cols != len(arr)
				 for (int i = 0; i < _rows(); i++)
					 for (int j = 0; j < _cols(); j++)
						 _set(i, j, M[i][j]);
			 }

			 T trace() {
				 T trace = 0;
				 for (int k = 0; k < _rows(); k++)
					 trace += this[k][k] ** 2;

				 return T;
			 }

			 T determinant(int n, int* rows, int* cols)
			 {
				 if (n == 2) {
					 T a = _get(rows[0], cols[0]);
					 T b = _get(rows[0], cols[1]);
					 T c = _get(rows[1], cols[0]);
					 T d = _get(rows[1], cols[1]);
					 return  a * d - b * c;
				 }
				 else
				 {
					 int* reduced_rows = new int[n - 1];
					 remove_dim(0, n, rows, reduced_rows);

					 T det = 0.;
					 int scalar;

					 for (int i = 0; i < n; i++)
					 {
						 int* reduced_cols = new int[n - 1];
						 remove_dim(i, n, cols, reduced_cols);

						 scalar = ((i * 3) % 2 == 0) ? 1 : -1;
						 det += scalar * _get(0, i) * determinant(n - 1, reduced_rows, reduced_cols);
						 delete reduced_cols;
					 }
					 delete reduced_rows;
					 return det;
				 }
			 }

			 T determinant()
			 {
				 int * rows = new int[_rows()];
				 int * cols = new int[_cols()];

				 for (int i = 0; i < _rows(); i++) { rows[i] = i; }
				 for (int i = 0; i < _cols(); i++) { cols[i] = i; }

				 T det = determinant(_rows(), rows, cols);
				 delete rows, cols;
				 return det;
			 }

			 square<T> cofactors()
			 {
				 if (_rows() > 2) {
					 T* cof = new T[_size];
					 int scalar;

					 int * full_rows = new int[_rows()];
					 int * full_cols = new int[_cols()];

					 for (int i = 0; i < _rows(); i++) { full_rows[i] = i; }
					 for (int i = 0; i < _cols(); i++) { full_cols[i] = i; }

					 for (int i = 0; i < _rows(); i++) {
						 int * reduced_rows = new int[_rows() - 1];
						 remove_dim(i, 3, full_rows, reduced_rows);

						 for (int j = 0; j < _cols(); j++) {
							 int * reduced_cols = new int[_cols() - 1];
							 remove_dim(j, 3, full_cols, reduced_cols);

							 scalar = ((i * 3 + j) % 2 == 0) ? 1 : -1;
							 cof[i * _rows() + j] = scalar * determinant(_rows() - 1, reduced_rows, reduced_cols);
							 delete[] reduced_cols;
						 }
						 delete[] reduced_rows;
					 }
					 square<T> _cof(_rows(), cof);

					 delete[] cof, full_rows, full_cols;
					 return _cof;
				 }
				 else
				 {
					 // raise exception: not possible to obtain cofactor matrix of a 2x2
				 }
			 }

			 square<T> inverse() { return cofactors().transpose() * (T(1.) / determinant()); }
		 }; // square matrix

		 template<typename T>
		 class vector : virtual public matrix<T>
		 {
		 private:
			 T _get(int i) { return _items[i]; }
			 void _set(int i, T val) { operator[](i) = val; }

		 public:
			 vector(const int len, ...) : _shape{ 1, len },
										  _size( len )
			 {
				 _allocate();

				 va_list args;
				 va_start(args, len);

				 for (int i = 0; i < len; i++)
					 _set(k, va_arg(args, T)); // doesn't seem to correctly convert to type

				 va_end(args);
			 }

			 vector(const int len, T val) : _shape{ 1, len },
											_size( len )
			 {
				 _allocate();

				 for (int i = 0; i < len; i++)
					 _set(i, val);
			 }

			 vector(const int len, T* arr) : _shape{ 1, len },
											 _size( len )
			 {
				 _allocate();

				 // TODO : Raise exception if rows * cols != len(arr)
				 for (int i = 0; i < len; i++)
					_set(1, arr[i]);
			 }

			 vector(matrix<T> M) : _shape{ 1, M.size() },
								   _size( M.size() )
			 {
				 _allocate();

				 // TODO : Raise exception if rows * cols != len(arr)
				 for (int i = 0; i < M.shape()[0]; i++)
					 for (int j = 0; j < M.shape()[1]; j++)
						 _set(i * M.shape()[0] + j, M[i][j]);
			 }

			 int length() { return _size; }
			
			 // Vector operations
			 T inner(vector<T> V)
			 {
				 // TODO: Check length of vectors match, throw if not

				 T result(0.);
				 for (int i = 0; i < length(); i++) { result += _get(i) * V[i]; }

				 return result;
			 }

			 matrix<T> outer(vector<T> V)
			 {
				 square<T> result(length(), T(0.));

				 for (int i = 0; i < length(); i++)
					 for (int j = 0; j < length(); j++)
						 result[i][j] = _get(i) * V[j];

				 return result;
			 }
		 };

		 template<typename T> 
		 vector<T> cross(vector<T> V1, vector<T> V2)
		 {
			 if (V1.length() == 3 && V2.length == 3)
			 {
				 return vector<T>(3,
					 V1[1] * V2[2] - V1[2] * V2[1],
					 V1[0] * V2[2] - V1[2] * V2[0],
					 V1[0] * V2[1] - V1[1] * V2[0]
					 );
			 
			 } else {
				 // raise exception
			 }
		 }

		 template<typename T>
		 square<T> skew(vector<T> V)
		 {
			 if (V.length() == 3)
			 {
				 return square<T>(3,
					 T(0.), -V[2],  V[1],
					  V[2], T(0.), -V[0],
					 -V[1], V[0],   T(0.)
					 );
			 }
			 else {
				 // raise exception
			 }
		 }
	}
}

template<typename T>
void display(attitude::linalg::matrix<T> M)
{
	for (int i = 0; i < M.shape()[0]; i++)
	{
		std::cout << "[ ";
		for (int j = 0; j < M.shape()[1]; j++)
			std::cout << M[i][j] << ",";

		std::cout << " ]" << std::endl;
	}
}
#endif // ATT_MATRIX_H
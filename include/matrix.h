// Copyright (c) 2020 Kieran Downie. All rights reserved.
//  
// This file is part of attitude.
//  
// attitude is free software : you can redistribute it and /
// or modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation,
// either version 3 of the License,
// or (at your option) any later version.
// 
// attitude is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public License
// along with attitude.  If not, see <https://www.gnu.org/licenses/>.
// 
// - 
// matrix.h
// Contains implementations for linear algebra formulations matrix and vector, 
// used in the attitude description parameter sets. Header-only.
// 
// File namespace attitude::linalg::
//
#ifndef ATT_MATRIX_H_
#define ATT_MATRIX_H_

#include <cstdarg>

namespace attitude {
namespace linalg {

// generate_submatrix_indices_along_dimension (function)
// Returns an array of indices that describe a submatrix formed by ommitting the
// perpindicular dimension at the speciifed index.
// Function inlined to [???]
// 
// Usage:
//	> int * full_row_indices = [0, 1, 2] // 3x3 matrix
//	> int * reduced_row_indices = generate_submatrix_indices_along_dimension(
//				2, // index of column to be removed from matrix description.
//				2, // size of row dimension.
//				rows // array with each of the full matrix row indices present.
//		); // returns [0, 3]
// 
// TODO: is inlining necessary? Complete doc comment.
// TODO: this is clumsy, need a rethink?
inline int * generate_submatrix_indices_along_dimension (
	int index_to_remove,
	const int dimension_size,
	int * current_indices
) {
  int j = 0;
	int * new_indices = new int[dimension_size - 1];

  for (int i = 0; i < dimension_size; i++)
    if (index_to_remove != current_indices[i]) {
      new_indices[j] = current_indices[i];
      j += 1;
    }
  return new_indices;
}

// matrix (class)
// Implementation of MxN matrix for attitude library. Matrix multiplication 
// & 'division' operations encoded. Templated so that can be defined for an
// element of any (primitive) type, precision can be application specific.
//
// Contruction based on the dimension sizes, followed by either a single 
// parameter, or MxN list (funciton parameters or c-style list, row-major):
//  - matrix<Tp> T(M, N, <MxN further args> );
//  - matrix<Tp> T(M, N, Tp(val) );
//  - matrix<Tp> T(M, N, carray[MxN] );
// 
// Usage:
//	> matrix<int> M(2, 2, 0, 1, -1, 0);
//  > M == M.transpose(); // returns true
//
template<typename Tp>
class matrix
{
 protected:
	const int shape_[2];
	const int size_;

	const int rows_ () { return shape_[0]; }
	const int cols_ () { return shape_[1]; }

	Tp * items_; // enforce size_ limit in memory?

	Tp get_ (int i, int j) { return items_[i * rows_() + j]; }

	void set_ (int i, int j, Tp val) { operator[](i)[j] = val; }
	void set_ (Tp * arr) { for (int i = 0; i <  size_; i++) { items_[i] = arr[i]; } }

	bool allocated_ = false; // Indicates whether memory has been allocated. TODO: check whether required now?

	void allocate_ () { if ( !allocated_ ) { items_ = new Tp[size_]; allocated_ = true; } }
	void deallocate_ () { if ( allocated_ ) { delete[] items_;  allocated_ = false; } }

 public:
	matrix(const int rows, const int cols, ...)
			: shape_{ rows, cols },
				size_(rows * cols)
	{
		allocate_();

		va_list args;
		va_start(args, cols);

		for (int k = 0; k < size_; k++)
			items_[k] = va_arg(args, Tp); // doesn't seem to correctly convert to type

		va_end(args);
	}

	matrix(const int rows, const int cols, Tp val)
			: shape_{ rows, cols },
				size_(rows * cols)
	{
		allocate_();

		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				items_[i * rows + j] = val;
	}

	matrix(const int rows, const int cols, Tp * arr)
			: shape_{ rows, cols },
				size_(rows * cols) 
	{
		allocate_();

		// TODO : Raise exception if rows * cols != len(arr)
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				set_(i, j, arr[i * rows + j]);
	}

	// Copy constructor
	matrix(matrix<Tp>& M)
			: shape_{ M.shape()[0], M.shape()[1] },
				size_(M.size())
	{
		allocate_();

		// TODO : Raise exception if rows * cols != len(arr)
		for (int i = 0; i < rows_(); i++)
			for (int j = 0; j < cols_(); j++)
				set_(i, j, M[i][j]);
	}

	~matrix() { deallocate_(); }

	int size() { return size_; }
	const int * shape() { return shape_; }

	// transpose (function)
	// Returns a new MxN transpose matrix (type Tp) of current NxM matrix instance.
	// If A = B^T, A_{i, j} = B_{j, i}
	// Accepts no arguments.
	// 
	// Usage:
	//	- matrix<Tp> M2 = M1.transpose();
	//
	matrix<Tp> transpose() 
	{
		const int new_rows = cols_();
    const int new_cols = rows_();

		// Must create new object to store new values to avoid overwriting old
    // values while needed.
    matrix<Tp> transposed(new_rows, new_cols);

		for (int i = 0; i < new_rows; i++)
			for (int j = 0; j < new_cols; j++)
				transposed[i][j] = get_(j, i);

		return transposed;
	}

	// -------------------- Elementwise operations. --------------------

  // + (function, operator)
  // Addition between two matrices computed on an elementwise basis. New matrix
	// created object and returned.
	// 
	matrix<Tp> operator+ (matrix<Tp> M)
	{
		matrix<Tp> result(rows_(), cols_());

		for (int i = 0; i < rows_(); i++)
			for (int j = 0; j < cols_(); j++)
				result[i][j] = get_(i, j) + M[i][j];

		return result;
	}

  // - (function, operator)
  // Subtraction between two matrices computed on an elementwise basis. New
  // matrix object created and returned.
	//
	matrix<Tp> operator- (matrix<Tp> M)
	{
		matrix<Tp> result(rows_(), cols_());

		for (int i = 0; i < rows_(); i++)
			for (int j = 0; j < cols_(); j++)
				result[i][j] = get_(i, j) - M[i][j];

		return result;
	}

  // + (function, operator)
	// Single value added to each element of the matrix instance. New matrix object
	// created and returned.
  //
	matrix<Tp> operator+ (Tp val)
	{
		matrix<Tp> result(rows_(), cols_());

		for (int i = 0; i < rows_(); i++)
			for (int j = 0; j < cols_(); j++)
				result[i][j] = operator[](i)[j] + val;

		return result;
	}

  // - (function, operator)
	// Single value subtracted from each element of the matrix instance. New matrix
	// object created and returned.
  //
	matrix<Tp> operator- (Tp val)
	{
		matrix<Tp> result(rows_(), cols_());

		for (int i = 0; i < rows_(); i++)
			for (int j = 0; j < cols_(); j++)
				result[i][j] = operator[](i)[j] - val;
					
		return result;
	}

  // * (function, operator)
	// Single factor applied to each matrix element. New matrix object created and
	// returned.
  //
	matrix<Tp> operator* (Tp val)
	{
		matrix<Tp> result(rows_(), cols_());

		for (int i = 0; i < rows_(); i++)
			for (int j = 0; j < cols_(); j++)
				result[i][j] = operator[](i)[j] * val;
					
		return result;
	}

  // / (function, operator)
	// Single divisor applied to each matrix element. New matrix object created and
	// returned.
  //
	matrix<Tp> operator/ (Tp val)
	{
		matrix<Tp> result(rows_(), cols_());

		for (int i = 0; i < rows_(); i++)
			for (int j = 0; j < cols_(); j++)
				result[i][j] = operator[](i)[j] / val;
					
		return result;
	}

	// += (function, operator)
  // Addition between two matrices computed on an elementwise basis. Current 
	// instance modified, pointer returned.
	// 
	matrix<Tp> * operator+= (matrix<Tp> M)
	{
		for (int i = 0; i < rows_(); i++)
			for (int j = 0; j < cols_(); j++)
				set_(i, j, get_(i, j) + M[i][j]);

		return this;
	}

  // -= (function, operator)
  // Subtraction between two matrices computed on an elementwise basis. Current 
	// instance modified, pointer returned.
	//
	matrix<Tp> * operator-= (matrix<Tp> M)
	{
		for (int i = 0; i < rows_(); i++)
			for (int j = 0; j < cols_(); j++)
				set_(i, j, get_(i, j) - M[i][j]);

		return this;
	}

  // += (function, operator)
	// Single value added to each element of the matrix instance. Current instance
	// modified, pointer returned.
  //
	matrix<Tp> * operator+= (Tp val)
	{
		for (int i = 0; i < rows_(); i++)
			for (int j = 0; j < cols_(); j++)
				set_(i, j, operator[](i)[j] + val);

		return this;
	}

  // -= (function, operator)
	// Single value subtracted from each element of the matrix instance. Current 
	// instance modified, pointer returned.
  //
	matrix<Tp> * operator-= (Tp val)
	{
		for (int i = 0; i < rows_(); i++)
			for (int j = 0; j < cols_(); j++)
				set_(i, j, operator[](i)[j] - val);
					
		return this;
	}

  // *= (function, operator)
	// Single factor applied to each matrix element. Current instance modified,
	// pointer returned.
  //
	matrix<Tp> * operator*= (Tp val)
	{
		for (int i = 0; i < rows_(); i++)
			for (int j = 0; j < cols_(); j++)
				set_(i, j, operator[](i)[j] * val);
					
		return this;
	}

  // /= (function, operator)
	// Single divisor applied to each matrix element. Current instance modiified,
	// pointer returned.
  //
	matrix<Tp> * operator/= (Tp val)
	{
		for (int i = 0; i < rows_(); i++)
			for (int j = 0; j < cols_(); j++)
				set_(i, j, operator[](i)[j] / val);
					
		return this;
	}

	// -------------------- Matrix operations. --------------------

	// * (function, operator)
  // Matrix multiplication. New instance created and returned.
  //
	matrix<Tp> operator* (matrix<Tp> M)
	{
		// TODO: raise exception if inner dimensions don't match
		const int shape[2]{ rows_(), M.shape()[1] };

		matrix<Tp> result(shape[0], shape[1]);

		for (int i = 0; i < shape[0]; i++) {
			for (int j = 0; j < shape[1]; j++) {
				result[i][j] = Tp(0.);
				for (int k = 0; k < cols_(); k++)
					result[i][j] += get_(i, k) * M[k][j];
			}
		}
		return result;
	}

  // *= (function, operator)
	// Matrix multiplication. Current instance modiified, pointer returned.
  //
	matrix<Tp> * operator*= (matrix<Tp> M)
	{
		// TODO: raise exception if inner dimensions don't match
		const int shape[2]{ rows_(), M.shape()[1] };

		// Must create new object to store new values to avoid overwriting old
		// values while needed.
		Tp * items = new T[shape[0] * shape[1]]{ T(0.) };

		for (int i = 0; i < shape[0]; i++) {
			for (int j = 0; j < shape[1]; j++) {
				for (int k = 0; k < cols_(); k++)
					items[i * shape[0] + j] += get_(i, k) * M[k][j];
			}
		}
				
		set_(items);
		delete items;
				
		return this;
	}

	// / (function, operator)
  // Matrix division computed as divisor's inverse applied to current instance.
	// i.e. [C] = [A]/[B] = [B]^{-1} * [A]
	// New instance created and returned.
  //
	matrix<Tp> operator/ (matrix<Tp> M) { return M.inverse() * *(this); }

	// /= (function, operator)
  // Matrix division computed as divisor's inverse applied to current instance.
	// i.e. [C] = [A]/[B] = [B]^{-1} * [A]
	// New instance created and pointer returned, current instance destructed[?].
  //
	matrix<Tp> * operator/= (matrix<Tp> M) { 
		matrix<Tp> inverted = M.inverse * *(this);
		return &inverted; // TODO: check that current instance (this) is deleted now.
	}

	// -------------------- Comparison operations. --------------------

  // == (function, operator)
  // Current instance is compared with right hand matrix. Returns true if all 
	// elements are the same.
	// i.e. true if all A_{i,j} == B{i, j}, otherwise false
	// 
	bool operator== (matrix<Tp> M)
	{
		// TODO: raise exception if dimensions don't match (or simply return false?)
		for (int i = 0; i < rows_(); i++)
			for (int j = 0; j < cols_(); j++)
				if (!(get_(i, j) == M[i][j])) { return false; }

		return true;
	}

  // == (function, operator)
  // Current instance is compared with right hand value. Returns true if all 
	// elements match the value.
	// i.e. true if all A_{i,j} == val, otherwise false
	// 
	bool operator== (Tp val)
	{
		for (int i = 0; i < rows_(); i++)
			for (int j = 0; j < cols_(); j++)
				if (!(get_(i, j) == val)) { return false; }

		return true;
	}

	// > (function, operator)
  // Current instance is compared with right hand value. Returns true if all 
	// elements are greater than the value.
	// i.e. true if all A_{i,j} > val, otherwise false
	// 
	bool operator> (Tp val)
	{
		for (int i = 0; i < rows_(); i++)
			for (int j = 0; j < cols_(); j++)
				if (!(get_(i, j) > val)) { return false; }

		return true;
	}

	// >= (function, operator)
  // Current instance is compared with right hand value. Returns true if all 
	// elements are greater than or equal to the value.
	// i.e. true if all A_{i,j} >= val, otherwise false
	// 
	bool operator>= (Tp val)
	{
		for (int i = 0; i < rows_(); i++)
			for (int j = 0; j < cols_(); j++)
				if (!(get_(i, j) >= val)) { return false; }

		return true;
	}

	// < (function, operator)
  // Current instance is compared with right hand value. Returns true if all 
	// elements are less than the value.
	// i.e. true if all A_{i,j} < val, otherwise false
	// 
	bool operator< (Tp val)
	{
		for (int i = 0; i < rows_(); i++)
			for (int j = 0; j < cols_(); j++)
				if (!(get_(i, j) < val)) { return false; }

		return true;
	}

	// <= (function, operator)
  // Current instance is compared with right hand value. Returns true if all 
	// elements are less than or equal to the value.
	// i.e. true if all A_{i,j} <= val, otherwise false
	// 
	bool operator<= (Tp val)
	{
		for (int i = 0; i < rows_(); i++)
			for (int j = 0; j < cols_(); j++)
				if (!(get_(i, j) <= val)) { return false; }

		return true;
	}

	// -------------------- Indexing. --------------------
	Tp * operator[] (int i) { return &items_[i * rows_()]; }

}; // matrix

// square (class)
// Implementation of NxN matrix for attitude library. Specific square matrix
// operations encoded. Templated so that can be defined for an element of any
// (primitive) type, precision can be application specific. Dervied from more
// general MxN matrix class.
// 
// Contruction based on the dimension size, followed by either a single
// parameter, or MxN list (funciton parameters or c-style list, row-major):
//  - square<Tp> T(N, <NxN further args> );
//  - square<Tp> T(N, Tp(val) );
//  - square<Tp> T(N, carray[NxN] );
//
// Usage:
//	> square<int> M(2, 2, 1, 0, 0, 1);
//  > 1 == M.determinant(); // returns true
//
template<typename Tp>
class square : virtual public matrix<Tp> {
	public:
	square(const int n, ...) 
			: matrix(n, n)
	{
		va_list args;
		va_start(args, n);

		for (int k = 0; k < size_; k++)
			items_[k] = va_arg(args, Tp);

		va_end(args);
	}
	square(const int n, Tp * arr) : matrix(n, n, arr) {}

	// convert base-class matrix into square
	square(matrix<Tp> M)
			: matrix(M.shape()[0], M.shape()[1])
	{
		allocate_();

		// TODO : Raise exception if rows * cols != len(arr)
		for (int i = 0; i < rows_(); i++)
			for (int j = 0; j < cols_(); j++)
				set_(i, j, M[i][j]);
	}

	// -------------------- Square Matrix Properties. --------------------

	// trace (function)
	// Returns trace (type Tp) of matrix instance.
	// Accepts no arguments.
	// 
	// Usage:
	//	- Tp tr = M.trace();
	//
	Tp trace() {
		// Trace is the sum of the square of each term on the matrix diagonal.
		Tp trace = 0;
		for (int k = 0; k < rows_(); k++)
			trace += this[k][k] ** 2;

		return trace;
	}

	// determinant (function)
	// Returns scalar-valued determinant (element type Tp) of the sub-matrix described
	// by arrays of row- & column- indices of the global matrix instance.
	//
	// Usage:
	//	- Tp det = M.determinant(N - 1, indices_of_rows, indices_of_cols);
	//
	private:
	Tp determinant(int n, int * rows, int * cols)
	{
		// Determinant of a 2x2 matrix is determinant as ad - bc. The determinant of 
		// larger matrices is found by iterating across the first row of the matrix, 
		// creating a sub-matrix description based on the elements remaining after 
		// removing the corresponding row & column for each element. Function is recursive.
		// 
		if (n == 2) {
			Tp a = get_(rows[0], cols[0]);
			Tp b = get_(rows[0], cols[1]);
			Tp c = get_(rows[1], cols[0]);
			Tp d = get_(rows[1], cols[1]);
			return  a * d - b * c;
		}
		else
		{
      int * reduced_rows = generate_submatrix_indices_along_dimension(0, n, rows);

			Tp det = 0.;
			int scalar;

			for (int i = 0; i < n; i++)
			{
				int * reduced_cols = generate_submatrix_indices_along_dimension(i, n, cols);

				scalar = ((i * 3) % 2 == 0) ? 1 : -1;
				det += scalar * get_(0, i) * determinant(n - 1, reduced_rows, reduced_cols);
				delete reduced_cols;
			}
			delete reduced_rows;
			return det;
		}
	}

	// determinant (function)
	// Returns scalar-valued determinant (element type Tp) of the current matrix instance.
	// Accepts no arguments.
	//
	// Usage:
	//	- Tp det = M.determinant();
	//
	public:
	Tp determinant()
	{
		int * rows = new int[rows_()];
		int * cols = new int[cols_()];

		for (int i = 0; i < rows_(); i++) { rows[i] = i; }
		for (int i = 0; i < cols_(); i++) { cols[i] = i; }

		// Entry point for lower-level (recursive) function defined above. Full matrix
		// description passed into next level.
		Tp det = determinant(rows_(), rows, cols);
		delete rows, cols;
		return det;
	}

	// cofactors (function)
	// Returns (square) matrix (element type Tp) of co-factors. Alternating factors
	// in [-1, 1] is applied to matrix of minors.
	// Accepts no arguments.
	// 
	// Usage:
	//	- square<Tp> cof = M.cofactors();
	// 
	square<Tp> cofactors()
	{
		if (rows_() > 2) {
			Tp * cof = new Tp[size_];
			int scalar;

			int * full_rows = new int[rows_()];
			int * full_cols = new int[cols_()];

			for (int i = 0; i < rows_(); i++) { full_rows[i] = i; }
			for (int i = 0; i < cols_(); i++) { full_cols[i] = i; }

			for (int i = 0; i < rows_(); i++) {
				int * reduced_rows = generate_submatrix_indices_along_dimension(i, rows_(), full_rows);

				for (int j = 0; j < cols_(); j++) {
					int * reduced_cols = generate_submatrix_indices_along_dimension(j, cols_(), full_cols);

					scalar = ((i * 3 + j) % 2 == 0) ? 1 : -1;
					cof[i * rows_() + j] = scalar * determinant(rows_() - 1, reduced_rows, reduced_cols);
					delete[] reduced_cols;
				}
				delete[] reduced_rows;
			}
			square<Tp> _cof(rows_(), cof);

			delete[] cof, full_rows, full_cols;
			return _cof;
		}
		else
		{
			// raise exception: not possible to obtain cofactor matrix of a 2x2
		}
	}

	// inverse (function)
	// Returns the inverse (element type Tp) of current matrix instance.
	// Accepts no arguments.
	// 
	// Usage:
	//	- square<Tp> inv = M.inverse();
	// 
	// Matrix inverse is (1 / det) * [C] where [C] is cofactor matrix.
	square<Tp> inverse() { return cofactors().transpose() * (Tp(1.) / determinant()); }

}; // square matrix


// vector (class)
// Implementation of a N-dimensional vector as an 1xN dimensional matrix representation
// for attitude library. Veector operations encoded as matrix operations. Templated so 
// it can be defined for an element of any (primitive) type, precision can be application 
// specific. Derived from more general MxN matrix class.
//
// Contruction based on the length, followed by either a single parameter, or MxN list
// (function parameters or c-style list, row-major):
//  - vector<Tp> T(N, <NxN further args> );
//  - vector<Tp> T(N, Tp(val) );
//  - vector<Tp> T(N, carray[NxN] );
//
// Usage:
//	> square<int> M(2, 2, 1, 0, 0, 1);
//  > 1 == M.determinant(); // returns true
//
template<typename Tp>
class vector : virtual public matrix<Tp>
{
	private:
	Tp get_(int i) { return items_[i]; }

	void set_(int i, Tp val) { operator[](i) = val; }

	public:
	vector(const int len, ...)
			: shape_{ 1, len },
				size_(len)
	{
		allocate_();

		va_list args;
		va_start(args, len);

		for (int i = 0; i < len; i++)
			set_(k, va_arg(args, Tp)); // TODO: doesn't seem to correctly enforce type

		va_end(args);
	}

	vector(const int len, Tp val)
			: shape_{ 1, len },
				size_(len)
	{
		allocate_();

		for (int i = 0; i < len; i++)
			set_(i, val);
	}

	vector(const int len, Tp * arr)
			: shape_{ 1, len },
				size_(len)
	{
		allocate_();

		// TODO : Raise exception if rows * cols != len(arr)
		for (int i = 0; i < len; i++)
		set_(1, arr[i]);
	}

	vector(matrix<Tp> M)
			: shape_{ 1, M.size() },
				size_( M.size() )
	{
		allocate_();

		// TODO : Raise exception if rows * cols != len(arr)
		for (int i = 0; i < M.shape()[0]; i++)
			for (int j = 0; j < M.shape()[1]; j++)
				set_(i * M.shape()[0] + j, M[i][j]);
	}

	int length() { return size_; }
			
	// -------------------- Vector Operations. --------------------

	// inner (function)
	// Returns inner (or dot) product of current instance and argument vectors.
	// Scalar result. Equivalent of matrix multiplication between vector's trans-
	// pose and itself.
	// i.e. [Minner] = v^T * v
	// 
	// Usage:
	//	- Tp vi = v.inner(i);
	//
	Tp inner(vector<Tp> V)
	{
		// TODO: Check length of vectors match

		Tp result(0.);
		for (int i = 0; i < length(); i++) { result += get_(i) * V[i]; }

		return result;
	}

	// outer (function)
	// Returns outer product of current instance and argument vector. Result is
	// symmetric matrix. Equivalent of the matrix multiplication of vector and it's
	// transpose.
	// i.e. [Mouter] = v * v^T
	// 
	// Usage:
	//	- square<Tp> M = i.outer(j);
	//
	matrix<Tp> outer(vector<Tp> V)
	{
		square<Tp> result(length(), T(0.));

		for (int i = 0; i < length(); i++)
			for (int j = 0; j < length(); j++)
				result[i][j] = get_(i) * V[j];

		return result;
	}
};

// cross (function)
// Returns cross product of two inputted vectors. Result is a new vector
//
// Usage:
//	- vector<Tp> k = cross(i, j);
//
template<typename Tp> 
vector<Tp> cross(vector<Tp> V1, vector<Tp> V2)
{
	if (V1.length() == 3 && V2.length == 3)
	{
		return vector<Tp>(3,
			V1[1] * V2[2] - V1[2] * V2[1],
			V1[0] * V2[2] - V1[2] * V2[0],
			V1[0] * V2[1] - V1[1] * V2[0]
			);
			 
	} else {
		// TODO: raise exception
	}
}

// tilde (function)
// Returns tilde matrix based inputted vector. This is the matrix
// equivalent pre-multiplier to the vector's cross product operation.
// i.e. cross(v1, v2) == tilde(v1) * v2
// 
// Usage:
//	- matrix<Tp> Mtilde = tilde(v1)
//
template<typename Tp>
square<Tp> tilde(vector<Tp> V)
{
	if (V.length() == 3) {
		return square<Tp>(3,
			Tp(0.), -V[2],  V[1],
			V[2], Tp(0.), -V[0],
			-V[1], V[0],   Tp(0.)
			);
	}
	else {
		// raise exception
	}
}

} // linalg
} // attitude

// display (function)
// Displays inputted matrix/square/vector in stdout. Printed line-by-line.
//
// Usage:
//	> square<int> I(3, 1, 0, 0, 0, 1, 0, 0, 0, 1);
//	> display(I); // displays identity matrix line-by-line in cout.
//
template<typename Tp>
void display(attitude::linalg::matrix<Tp> M)
{
	for (int i = 0; i < M.shape()[0]; i++)
	{
		std::cout << "[ ";
		for (int j = 0; j < M.shape()[1]; j++)
			std::cout << M[i][j] << ",";

		std::cout << " ]" << std::endl;
	}
}

#endif // ATT_MATRIX_H_

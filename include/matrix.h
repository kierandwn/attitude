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
// File namespace attitude::
//
#ifndef ATT_MATRIX_H_
#define ATT_MATRIX_H_

#include <initializer_list>
#include <algorithm>

namespace attitude {


// matrix (class)
// MxN matrix implementation. Matrix multiplication & 'division' operations 
// defined. Templated so that can be defined for an element of any (primitive)
// type, precision can be application specific.
template <typename Tp, size_t rows, size_t cols>
class matrix {
  // ==================== matrix ====================
 private:
  Tp get_(int i) {
    if (i < size_) {
      return items_[i];
    } /* else { throw InvalidIndex; } */
  }
  Tp get_(int i, int j) {
    if (i < shape_[0] && j < shape_[1]) {
      return items_[i * shape_[1] + j];
    } /* else { throw InvalidIndex; } */
  }

  void set_(int i, int j, Tp val) {
    if (i < shape_[0] && j < shape_[1]) {
      items_[i * shape_[1] + j] = val;
    } /* else { throw InvalidIndex; } */
  }

  void set_(int i, Tp val) {
    if (i < size_) {
      items_[i] = val;
    } /* else { throw InvalidIndex; } */
  }

  void set_(Tp* arr) {
    for (int i = 0; i < size_; i++) {
      items_[i] = arr[i];
    }
  }

 protected:
  Tp items_[rows * cols];

  const size_t size_ = rows * cols;
  const size_t shape_[2]{rows, cols};

 public:
  matrix() {}
  matrix(Tp arr[rows * cols]) { set_(arr); }

  matrix(std::initializer_list<Tp> il) {
    if (il.size() == 1) {
      for (int i = 0; i < size_; ++i) {
        items_[i] = *(il.begin());
      }

    } else if (il.size() == size_) {
      std::copy(il.begin(), il.end(), items_);

    } else {
      throw("List initialiser length does not match matrix dimensions.");
    }
  }

  template <typename Tp2>
  matrix(std::initializer_list<Tp2> il) {
    if (il.size() == 1) {
      for (int i = 0; i < size_; ++i) {
        items_[i] = *(il.begin());
      }

    } else if (il.size() == size_) {
      std::copy(il.begin(), il.end(), items_);

    } else {
      throw("List initialiser length does not match matrix dimensions.");
    }
  }

  template <typename Tp2, size_t rows2>
  matrix(matrix<Tp2, rows2, (rows * cols) / rows2>& M) {
    for (int i = 0; i < size_; i++) {
      set_(i, M(i));
    }
  }

  template <typename Tp2, size_t rows2>
  matrix(matrix<Tp2, rows2, (rows * cols) / rows2>* M) {
    for (int i = 0; i < size_; i++) {
      set_(i, M->operator()(i));
    }
  }

  const size_t size() { return size_; }
  const size_t* shape() { return &shape_[0]; }

  // transpose (function)
  // Returns a new MxN transpose matrix (type Tp) of current NxM matrix
  // instance. If A = B^T, all A_{i, j} = B_{j, i}.
  matrix<Tp, cols, rows> transpose() {
    matrix<Tp, cols, rows> transposed;

    // Must create new object to store new values to avoid overwriting old
    // values while needed.
    for (int i = 0; i < shape_[0]; i++)
      for (int j = 0; j < shape_[1]; j++)
        transposed[j][i] = items_[i * cols + j];

    return transposed;
  }

  // -------------------- Elementwise operations. --------------------

  // + (function, operator)
  // Elementwise matrix addition. New object created & returned.
  matrix<Tp, rows, cols> operator+(matrix<Tp, rows, cols> M) {
    matrix<Tp, rows, cols> result(this);

    for (int i = 0; i < size_; ++i) result(i) += M(i);
    return result;
  }

  // - (function, operator)
  // Elementwise matrix subtraction. New object created & returned.
  matrix<Tp, rows, cols> operator-(matrix<Tp, rows, cols> M) {
    matrix<Tp, rows, cols> result(this);

    for (int i = 0; i < size_; ++i) result(i) -= M(i);
    return result;
  }

  // + (function, operator)
  // Single value added to each matrix element. New object created & returned.
  matrix<Tp, rows, cols> operator+(Tp val) {
    matrix<Tp, rows, cols> result(this);

    for (int i = 0; i < size_; ++i) result(i) += val;
    return result;
  }

  // - (function, operator)
  // Single value subtracted from each matrix element. New object created &
  // returned.
  matrix<Tp, rows, cols> operator-(Tp val) {
    matrix<Tp, rows, cols> result(this);

    for (int i = 0; i < size_; ++i) result(i) -= val;
    return result;
  }

  // * (function, operator)
  // Single factor applied to each matrix element. New object created &
  // returned.
  matrix<Tp, rows, cols> operator*(Tp val) {
    matrix<Tp, rows, cols> result(this);

    for (int i = 0; i < size_; ++i) result(i) *= val;
    return result;
  }

  // / (function, operator)
  // Single divisor applied to each matrix element. New object created &
  // returned.
  matrix<Tp, rows, cols> operator/(Tp val) {
    matrix<Tp, rows, cols> result(this);

    for (int i = 0; i < size_; ++i) result(i) /= val;
    return result;
  }

  // += (function, operator)
  // Elementwise matrix addition. Current instance modified, pointer returned.
  matrix<Tp, rows, cols>* operator+=(matrix<Tp, rows, cols> M) {
    for (int i = 0; i < size_; ++i) set_(i, get_(i) + M(i));
    return this;
  }

  // -= (function, operator)
  // Elementwise matrix subtraction. Current instance modified, pointer
  // returned.
  matrix<Tp, rows, cols>* operator-=(matrix<Tp, rows, cols> M) {
    for (int i = 0; i < size_; ++i) set_(i, get_(i) - M(i));
    return this;
  }

  // += (function, operator)
  // Single value added to each matrix element. Current instance modified,
  // pointer returned.
  matrix<Tp, rows, cols>* operator+=(Tp val) {
    for (int i = 0; i < size_; ++i) set_(i, get_(i) + val);
    return this;
  }

  // -= (function, operator)
  // Single value subtracted from each matrix element. Current instance
  // modified, pointer returned.
  matrix<Tp, rows, cols>* operator-=(Tp val) {
    for (int i = 0; i < size_; ++i) set_(i, get_(i) - val);
    return this;
  }

  // *= (function, operator)
  // Single factor applied to each matrix element. Current instance modified,
  // pointer returned.
  matrix<Tp, rows, cols>* operator*=(Tp val) {
    for (int i = 0; i < size_; ++i) set_(i, get_(i) * val);
    return this;
  }

  // /= (function, operator)
  // Single divisor applied to each matrix element. Current instance modiified,
  // pointer returned.
  //
  matrix<Tp, rows, cols>* operator/=(Tp val) {
    for (int i = 0; i < size_; ++i) set_(i, get_(i) / val);
    return this;
  }

  // power (function)
  // Elementise power function. Raises each element to specified exponent.
  matrix<Tp, rows, cols> power(int exp) {
    matrix<Tp, rows, cols> result(*(this));

    for (int i = 0; i < size_; ++i) { result(i) = pow(get_(i), exp); }
    return result;
  }
  // -------------------- Matrix operations. --------------------

  // * (function, operator)
  // Matrix multiplication. New instance created and returned.
  template <size_t cols2>
  matrix<Tp, rows, cols2> operator* (matrix<Tp, cols, cols2> M) {
    const int shape[2]{shape_[0], M.shape()[1]};
    matrix<Tp, rows, cols2> result{0.};

    for (int i = 0; i < shape[0]; ++i) {
      for (int j = 0; j < shape[1]; ++j) {
        for (int k = 0; k < shape_[1]; ++k)
          result[i][j] += get_(i, k) * M[k][j];
      }
    }
    return result;
  }

  // *= (function, operator)
  // Matrix multiplication. Current instance modiified, pointer returned.
  template <size_t cols2>
  matrix<Tp, rows, cols2>* operator*=(matrix<Tp, cols, cols2> M) {
    const int shape[2]{shape_[0], M.shape()[1]};

    // Must create new object to store new values to avoid overwriting old
    // values while needed.
    Tp items[rows * cols2]{0.};

    for (int i = 0; i < shape[0]; ++i) {
      for (int j = 0; j < shape[1]; ++j) {
        for (int k = 0; k < shape_[1]; ++k)
          items[i * shape[1] + j] += get_(i, k) * M[k][j];
      }
    }

    set_(items);
    return this;
  }

  // / (function, operator)
  // Matrix division computed as divisor's inverse applied to current instance.
  // i.e. [C] = [A]/[B] = [B]^{-1} * [A].
  // Only possible with square right-hand matrix.
  // New instance created and returned.
  matrix<Tp, rows, cols> operator/(matrix<Tp, rows, rows> M) {
    return inverse(M) * *(this);
  }

  // /= (function, operator)
  // Matrix division computed as divisor's inverse applied to current instance.
  // i.e. [C] = [A]/[B] = [B]^{-1} * [A]
  // Only possible with square right-hand matrix.
  // New instance created and pointer returned, current instance destructed.
  matrix<Tp, rows, cols>* operator/=(matrix<Tp, rows, rows> M) {
    matrix<Tp, rows, cols> result = inverse(M) * *(this);
    *(this) = result;
    return this;
  }

  // -------------------- Comparison operations. --------------------

  // == (function, operator)
  // Returns true if all elements match the value.
  bool operator==(Tp val) {
    for (int i = 0; i < size_; ++i)
      if (!(get_(i) == val)) { return false; } 

    return true;
  }

  // == (function, operator)
  // Returns true if all elements are the same.
  bool operator==(matrix<Tp, cols, rows> M) {
    for (int i = 0; i < size_; ++i)
      if (!(get_(i) == M(i))) { return false; }
     
    return true; 
  }

  // > (function, operator)
  // Returns true if all elements are greater than the value.
  bool operator>(Tp val) {
    for (int i = 0; i < size_; ++i)
      if (!(get_(i) > val)) { return false; }

    return true;
  }

  // > (function, operator)
  // Returns true if all elements are greater than all corresponding right-hand
  // elements.
  bool operator>(matrix<Tp, cols, rows> M) {
    for (int i = 0; i < size_; ++i)
      if (!(get_(i) > M(i))) { return false; }

    return true;
  }

  // >= (function, operator)
  // Returns true if all elements are greater than or equal to the value.
  bool operator>=(Tp val) {
    for (int i = 0; i < size_; ++i)
      if (!(get_(i) >= val)) { return false; }

    return true;
  }

  // >= (function, operator)
  // Returns true if all elements are greater than or equal to all
  // corresponding right-hand elements.
  bool operator>=(matrix<Tp, cols, rows> M) {
    for (int i = 0; i < size_; ++i)
      if (!(get_(i) >= M(i))) { return false; }

    return true;
  }

  // < (function, operator)
  // Returns true if all elements are less than the value.
  bool operator<(Tp val) {
    for (int i = 0; i < size_; ++i)
      if (!(get_(i) < val)) { return false; }
     
    return true;
  }

  // < (function, operator)
  // Returns true if all elements are less than all corresponding right-hand
  // elements.
  bool operator<(matrix<Tp, cols, rows> M) {
    for (int i = 0; i < size_; ++i)
      if (!(get_(i) < M(i))) { return false; }

    return true;
  }

  // <= (function, operator)
  // Returns true if all elements are less than or equal to the value.
  bool operator<=(Tp val) {
    for (int i = 0; i < size_; ++i)
      if (!(get_(i) <= val)) { return false; }

    return true;
  }

  // <= (function, operator)
  // Returns true if all elements are less than or equal to all corresponding
  // right-hand elements.
  bool operator<=(matrix<Tp, cols, rows> M) {
    for (int i = 0; i < size_; ++i)
      if (!(get_(i) <= M(i))) { return false; }

    return true;
  }

  // -------------------- Indexing. --------------------
  Tp & operator()(int i) { return items_[i]; }  // row major linear index. enables A(i) = val;
  Tp * operator[](int i) { return &items_[i * shape_[1]]; }  // enables A[i][j] = val;

  // -------------------- Assignment. --------------------
  matrix<Tp, rows, cols>& operator= (matrix<Tp, rows, cols>& M) {
    for (int i = 0; i < size_; ++i) { set_(i, M(i)); }
    return *(this);
  }

};  // matrix

// remove_rowi_colj (function)
// Returns scalar-valued matrix determinant (type Tp)
template <typename Tp, size_t n>
matrix<Tp, n - 1, n - 1> remove_rowi_colj(size_t i, size_t j,
                                          matrix<Tp, n, n> square) {
  Tp reduced[(n - 1) * (n - 1)];
  int kk = 0;

  for (int k = 0; k < square.size(); ++k) {
    if (!((k / n == i) || (k % n == j))) {
      reduced[kk] = square(k);
      kk += 1;
    }
  }
  return matrix<Tp, n - 1, n - 1>(reduced);
}

// ==================== Matrix Properties. ====================
// ================== (non-member functions) ==================

// trace (function)
// Returns trace (type Tp) of matrix instance.
template <typename Tp, size_t n>
Tp trace(matrix<Tp, n, n> square) {
  Tp trace = 0.;
  for (int i = 0; i < n; ++i) trace += square[i][i];
  return trace;
}

// determinant (function)
// Returns scalar-valued matrix determinant (type Tp). Recursive.
// Determinant is sum along one dimension of matrix of minors.
template <typename Tp, size_t n>
Tp determinant(matrix<Tp, n, n> square) {
  Tp det = 0.;
  int scalar;

  for (int k = 0; k < n; ++k) {
    scalar = (k % 2 == 0) ? 1 : -1;
    det += scalar * square[0][k] * determinant(remove_rowi_colj(0, k, square));
  }
  return det;
}

// Base case for matrix determinant function: 2x2.
// Determinant is ad - bc.
template <typename Tp>
Tp determinant(matrix<Tp, 2, 2> square) {
  return square[0][0] * square[1][1] - square[0][1] * square[1][0];
}

// cofactors (function)
// Returns matrix of co-factors. Alternating factors
// in [-1, 1] is applied to matrix of minors.
template <typename Tp, size_t n>  // TODO: n > 2 enforced?
matrix<Tp, n, n> cofactors(matrix<Tp, n, n> square) {
  matrix<Tp, n, n> cofactor_matrix;
  int scalar;

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      scalar = ((i * 3 + j) % 2 == 0) ? 1 : -1;

      cofactor_matrix[i][j] =
          scalar * determinant(remove_rowi_colj(i, j, square));
    }
  }
  return cofactor_matrix;
}

// inverse (function)
// Returns the inverse (element type Tp) of current matrix instance.
// Matrix inverse is (1 / det) * [C] where [C] is cofactor matrix.
template <typename Tp, size_t n>
matrix<Tp, n, n> inverse(matrix<Tp, n, n> square) {
  return cofactors(square).transpose() * (Tp(1.) / determinant(square));
}



// vector (class)
// Nx1 matrix representation of vector. Vector operations encoded.
// Derived from matrix class.
template <typename Tp, size_t len>
class vector : virtual public matrix<Tp, len, 1> {
  // ==================== vector ====================
 private:
  Tp get_(int i) {
    if (i < size_) {
      return items_[i];
    } /* else { throw InvalidIndex; } */
  }

  void set_(int i, Tp val) {
    if (i < size_) {
      items_[i] = val;
    } /* else { throw InvalidIndex; } */
  }

  void set_(Tp arr[len]) {
    for (int i = 0; i < size_; ++i) {
      items_[i] = arr[i];
    }
  }

 public:
  vector() {}
  vector(Tp arr[len]) { set_(arr); }

  vector(std::initializer_list<Tp> il) {
    if (il.size() == 1) {
      for (int i = 0; i < size_; ++i) {
        items_[i] = *(il.begin());
      }

    } else if (il.size() == size_) {
      std::copy(il.begin(), il.end(), items_);

    } else {
      throw("List initialiser length does not match matrix dimensions.");
    }
  }

  vector(matrix<Tp, len, 1> M) {
    for (int i = 0; i < M.size(); ++i) set_(i, M(i));
  }

  template <typename Tp2, size_t cols>
  vector(matrix<Tp2, len / cols, cols> M) {
    for (int i = 0; i < M.size(); ++i) set_(i, M(i));
  }

  const size_t length() { return size_; }

  Tp norm(){
    Tp norm_sq = Tp(0.);
    for (int i = 0; i < len; ++i) { norm_sq += pow(get_(i), 2); }
    return sqrt(norm_sq);
  }

  // -------------------- Vector Operations. --------------------

  // inner (function)
  // Returns inner (or dot) product of left & right-hand vectors. Scalar result.
  Tp inner(vector<Tp, len> V) {
    Tp result(0.);
    for (int i = 0; i < length(); ++i) result += get_(i) * V[i];
    return result;
  }

  // outer (function)
  // Returns outer product of left & right hand vector. Symmetric matrix result.
  matrix<Tp, len, len> outer(vector<Tp, len> V) {
    matrix<Tp, len, len> result{0.};

    for (int i = 0; i < length(); ++i)
      for (int j = 0; j < length(); ++j) result[i][j] = get_(i) * V[j];

    return result;
  }

  // -------------------- Indexing. --------------------
  Tp & operator[](int i) { return items_[i]; }  // enables A[i][j] = val;


};  // vector


// cross (function)
// Returns cross product of two inputted vectors. Result is a new vector
//
template <typename Tp>
vector<Tp, 3> cross(vector<Tp, 3> lhs, vector<Tp, 3> rhs) {
  return vector<Tp, 3>{
    lhs[1] * rhs[2] - lhs[2] * rhs[1],
    -1 * (lhs[0] * rhs[2] - lhs[2] * rhs[0]),
    lhs[0] * rhs[1] - lhs[1] * rhs[0]
  };
}

// tilde (function)
// Returns tilde matrix based inputted vector. This is the matrix
// equivalent pre-multiplier to the vector's cross product operation.
// i.e. cross(v1, v2) == tilde(v1) * v2
//
template <typename Tp>
matrix<Tp, 3, 3> tilde(vector<Tp, 3> vec) {
  return matrix<Tp, 3, 3>{
    Tp(0.), -vec[2],  vec[1], 
    vec[2],  Tp(0.), -vec[0],
   -vec[1],  vec[0], Tp(0.)
  };
}

// eye (function)
// Returns nxn identity matrix. 
//
template <typename Tp, size_t n>
matrix<Tp, n, n> eye() {
  matrix<Tp, n, n> identity{Tp(0.)};
  for (int i = 0; i < n; ++i) {
    identity[i][i] = 1.;
  }
  return identity;
}


}  // namespace attitude

template <typename Tp, size_t r, size_t c>
void display(attitude::matrix<Tp, r, c> M) {
  printf("[ ");
  for (int i = 0; i < M.size(); ++i) {
    printf("%.8f, ", M(i));
    if (((i + 1) % M.shape()[1]) == 0) {
      printf("] \n");
      if ((i + 1) != M.size()) {
        printf("[ ");
      }
    }
  }
}

#endif ATT_MATRIX_H_

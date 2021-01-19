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
// base.h
// Contains abstract base definition for attitude description parameter sets.
//
// File namespace attitude::
//
#ifndef ATT_BASE_H_
#define ATT_BASE_H_

#include <initializer_list>
#include <iostream>

#include "matrix.h"

namespace attitude {


// description_set (class, abstract)
// Desribes commonality between all attitude parameterisations, and indicates 
// specific functionalities to implemented in each subclass.
template<typename Tp, size_t n_items>
class description_set {
 protected:
  vector<Tp, n_items> items_;
  
  Tp get_(int i) { return items_[i]; } 
  void set_(int i, Tp val) { items_[i] = val; } 
  void set_(Tp vals[n_items]) { items_ = vals; } 
  void set_(vector<Tp, n_items> vals) { items_ = vals; } 

  bool update_dcm_ = false;
  virtual void dcm_from_parameters_() = 0;
  virtual void parameters_from_dcm_() = 0;

  matrix<Tp, 3, 3> matrix_;

  description_set() {}
  description_set(matrix<Tp, 3, 3> M) 
      : matrix_(M) {}

  description_set(std::initializer_list<Tp> il) {
    if (il.size() == n_items) {
      Tp arr[n_items];
      std::copy(il.begin(), il.end(), arr); // can't directly copy in vector class, TODO
      items_ = vector<Tp, n_items>(arr);

    } else {
      throw("List initialiser length does not match matrix dimensions.");
    }
  }

  template<typename Tp2, size_t n2_items>
  description_set(description_set<Tp2, n2_items> M) 
      : matrix_(M.matrix()) {}

 public:
  // dke (function)
  // Differential Kinematic Equation.
  // Pure virtual function. Indicate concrete descriptions sets must specify
  // a matrix for mapping angular velocity onto their parameter rates.
  virtual matrix<Tp, n_items, n_items> dke() = 0;

  // matrix (function)
  // Returns the directional cosine matrix equivalent of the represented rotation.
  matrix<Tp, 3, 3> matrix() {
    if (update_dcm_) { dcm_from_parameters_(); }
    return matrix_;
  }

  // reverse (function)
  // Returns reverse rotation of current. auto keyword used so return DCM can be 
  // to any attitude representation.
  auto reverse() { return matrix().transpose(); }

  // propagate (function)
  // Given angular velocity, propagate current representation over timestep dt.
  void propagate(Tp dt, vector<Tp, 3> omega) {
    vector<Tp, n_items> p_rate = dke() * omega;
    for (int i = 0; i < n_items; ++i) { items_[i] = get_(i) + dt * dp[i]; }
  }

  // -------------------- DCM-based Addition/Subtraction. --------------------
  // TODO: right now, cannot instantiate abstract class to 'set' - but how can I define 
  // DCM based addition operations without writing it in every subclass?
  /*template <typename Tp2, size_t n2_items>
  auto operator+ (description_set<Tp2, n2_items> set) {
    return matrix() * set.matrix();
  }

  template <typename Tp2, size_t n2_items>
  auto operator- (description_set<Tp2, n2_items> set) {
    return matrix() / set.matrix();
  }

  template <typename Tp2, size_t n2_items>
  description_set<Tp, n_items> * operator+= (description_set<Tp2, n2_items> set) {
    matrix_ = matrix() * set.matrix();
    parameters_from_dcm_();
    return this;
  }

  template <typename Tp2, size_t n2_items>
  description_set<Tp, n_items> * operator-= (description_set<Tp2, n2_items> set) {
    matrix_ = matrix() / set.matrix();
    parameters_from_dcm_();
    return this;
  }*/

  // -------------------- Comparison. --------------------
  template<size_t n2_items>
  bool operator== (description_set<Tp, n2_items> set) { return matrix() == set.matrix(); }
  // bool operator== (::matrix<Tp, 3, 3> & rhs) { return matrix() == rhs; } // TODO: no effect?

  // -------------------- Indexing. --------------------
  Tp & operator[] (int i) { update_dcm_ = true; return items_[i]; }
};


// display (function)
// STDOUT display for description sets. Items are display inline enclosed in curly braces.
template <typename Tp, size_t n_items>
void display(attitude::description_set<Tp, n_items> & p) {
  printf("{ %.8f", p[0]);
  for (int i = 1; i < n_items; ++i) {
    printf(", %.8f", p[i]);
  }
  std::cout << " }" << std::endl;
}
  

} // namespace attitude
#endif // ATT_BASE_H_
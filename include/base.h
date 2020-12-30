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

#include "matrix.h"

namespace attitude {


template<typename Tp, size_t n_items>
class description_set {
 protected:
  vector<Tp, n_items> items_;
  
  Tp get_(int i) { return items_[i]; } 
  void set_(int i, Tp val) { items_[i] = val; } 

  bool update_dcm_ = true;
  virtual void dcm_from_parameters_() = 0;
  virtual void parameters_from_dcm_() = 0;

  matrix<Tp, 3, 3> matrix_;

  description_set() {}
  description_set(matrix<Tp, 3, 3> M) 
      : matrix_(M) {}
  //{ parameters_from_dcm_(); }

  description_set(std::initializer_list<Tp> il) {
    if (il.size() == n_items) {
      std::copy(il.begin(), il.end(), items_);

    } else {
      throw("List initialiser length does not match matrix dimensions.");
    }
  }

  template<typename Tp2, size_t n2_items>
  description_set(description_set<Tp2, n2_items> M) 
      : matrix_(M.matrix()) {}
  //{ parameters_from_dcm_(); }

 public:
  virtual matrix<Tp, n_items, n_items> dke() = 0;

  matrix<Tp, 3, 3> matrix() {
    if (update_dcm_) { dcm_from_parameters_(); }
    return matrix_;
  }

  auto reverse() { return matrix().transpose(); }

  void propagate(Tp dt, vector<Tp, 3> omega) {
    vector<Tp, n_items> p_rate = dke() * omega;
    for (int i = 0; i < n_items; ++i) { items_[i] += dt * dp[i]; }
  }

  // -------------------- DCM-based Addition/Subtraction. --------------------
  // TODO: right now, cannot instantiate abstract class to 'set' - how can I define 
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
  bool operator== (::matrix<Tp, 3, 3> rhs) { return matrix() == rhs; }

  // -------------------- Indexing. --------------------
  Tp & operator[] (int i) { update_dcm_ = true; return items_[i]; }
};
  

} // namespace attitude
#endif // ATT_BASE_H_
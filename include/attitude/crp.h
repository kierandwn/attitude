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
// crp.h
// Contains concrete definition for classical rodriguez attitude
// description parameter set.
//
// File namespace attitude::
//
#ifndef ATT_CRODR_H_
#define ATT_CRODR_H_

#include <cmath>
#include <initializer_list>

#include "attitude/base.h"
#include "attitude/matrix.h"
#include "attitude/shepherds.h"

namespace attitude {
// TODO: a lot of the below relies on the vector operations. More efficient
// to code up 'long-hand' to avoid creating many matrix/vecotr objects unnecessarily?


// crp (class)
// Classical Rodriguez Parameters (attitude representation).
// Decribes direct addition & subtraction methods for cl. rodriguez parameters,
// and defines mapping to/from DCM, and angular velocity to CRP rates.
template <typename Tp>
class crp : public virtual description_set<Tp, 3> 
{
  using description_set<Tp, 3>::get_;
  // using description_set<Tp, 3>::set_;
  using description_set<Tp, 3>::items_;
  using description_set<Tp, 3>::matrix_;
  using description_set<Tp, 3>::update_dcm_;

 public:
  crp() {}

  crp(Tp q1, Tp q2, Tp q3) : description_set{q1, q2, q3} { dcm_from_parameters_(); }
  crp(::mn_matrix<Tp, 3, 3> R) : description_set(R) { parameters_from_dcm_(); }
  crp(::vector<Tp, 3> q) : description_set{q[0], q[1], q[2]} { dcm_from_parameters_(); }

  template <typename Tp2, size_t n2_items>
  crp(description_set<Tp2, n2_items>* set) : description_set(set) {}

  // dke (function)
  // Returns a new 3x3 matrix (type Tp) that maps angular velocity onto classical
  // rodriguez rates.
  ::mn_matrix<Tp, 3, 3> dke() override {
    return ::mn_matrix<Tp, 3, 3>{
      1. + pow(get_(0), 2), get_(0) * get_(1) - get_(2), get_(0) * get_(2) + get_(1),
      get_(0) * get_(1) + get_(2), 1. + pow(get_(1), 2), get_(1) * get_(2) - get_(0),
      get_(0) * get_(2) - get_(1), get_(1) * get_(2) + get_(0), 1. + pow(get_(2), 2)
    } * 0.5;
  }

  // reverse (function)
  // Returns the reverse rotation represented as a CRP.
  crp<Tp> reverse() { return crp<Tp>(items_ * -1); }


  // -------------------- Classical Rodriguez Addition/Subtraction. --------------------
  crp<Tp> operator+(crp<Tp> rhs) {
    vector<Tp, 3> rhs_vec{ rhs[0], rhs[1], rhs[2] };
    return (items_ + rhs_vec - cross(items_, rhs_vec)) / (1 - items_.inner(rhs_vec));
    // TODO: add vector to constructor of description sets.
    // TODO: also add a PRV object: inheriting from vector, description_set
  }

  crp<Tp> operator-(crp<Tp> rhs) { // TODO
    vector<Tp, 3> rhs_vec{rhs[0], rhs[1], rhs[2]};
    return (items_ - rhs_vec + cross(items_, rhs_vec)) / (1 + items_.inner(rhs_vec));
  }

  crp<Tp> * operator+=(crp<Tp> rjs) {
    vector<Tp, 3> rhs_vec{rhs[0], rhs[1], rhs[2]};
    items_ = (items_ + rhs_vec - cross(items_, rhs_vec)) / (1 - items_.inner(rhs_vec));
    return this;
  }

  crp<Tp> * operator-=(crp<Tp> rhs) { // TODO
    vector<Tp, 3> rhs_vec{rhs[0], rhs[1], rhs[2]};
    items_ = (items_ - rhs_vec + cross(items_, rhs_vec)) / (1 + items_.inner(rhs_vec));
    return this;
  }


  // -------------------- Comparison. --------------------
  bool operator==(crp<Tp> q) { return matrix() == q.matrix(); }


 private:
  // -------------------- Description set virtual overrides. --------------------

  // parameters_from_dcm_ (function)
  // Computes the classical rodriguez components from directional cosine matrix.
  void parameters_from_dcm_() override {
    items_ = vector<Tp, 3>{
        matrix_[1][2] - matrix_[2][1],
        matrix_[2][0] - matrix_[0][2],
        matrix_[0][1] - matrix_[1][0]
    } * (1. / (1. + trace(matrix_)));
    update_dcm_ = false; // dcm & parameters are now consistent
  }

  // dcm_from_parameters_ (function)
  // Computes the directional cosine matrix from classical rodriguez components.
  void dcm_from_parameters_() override {
    Tp norm2 = pow(get_(0), 2) + pow(get_(1), 2) + pow(get_(2), 2);

    matrix_ = ::mn_matrix<Tp, 3, 3>{
      // Row 1
      1 + pow(get_(0), 2) - pow(get_(1), 2) - pow(get_(2), 2),
      2. * (get_(0) * get_(1) + get_(2)), 
      2. * (get_(0) * get_(2) - get_(1)),
      // Row 2
      2. * (get_(0) * get_(1) - get_(2)), 
      1 - pow(get_(0), 2) + pow(get_(1), 2) - pow(get_(2), 2),
      2. * (get_(1) * get_(2) + get_(0)),
      // Row 3
      2. * (get_(0) * get_(2) + get_(1)), 
      2. * (get_(1) * get_(2) - get_(0)), 
      1 - pow(get_(0), 2) - pow(get_(1), 2) + pow(get_(2), 2)
    } * (1. / (1. + norm2));
    update_dcm_ = false; // dcm & parameters are now consistent
  }
};


}  // namespace attitude
#endif  // ATT_CRODR_H_
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
// mrp.h
// Contains concrete definition for modified rodriguez attitude
// description parameter set.
//
// File namespace attitude::
//
#ifndef ATT_MRODR_H_
#define ATT_MRODR_H_


#include "base.h"
#include "matrix.h"
#include "shepherds.h"

namespace attitude {
// TODO: a lot of the below relies on the vector operations. More efficient 
// to code up 'long-hand' to avoid creating many matrix/vecotr objects unnecessarily?


// mrp (class)
// Modified Rodriguez Parameters (attitude representation).
// Decribes direct addition & subtraction methods for m. rodriguez parameters,
// and defines mapping to/from DCM, and angular velocity to MRP rates.
template <typename Tp>
class mrp : public virtual description_set<Tp, 3> {
 public:
  mrp() {}

  mrp(Tp s1, Tp s2, Tp s3) : description_set{s1, s2, s3} { dcm_from_parameters_(); }

  mrp(::matrix<Tp, 3, 3> R) : description_set(R) { parameters_from_dcm_(); }
  mrp(::vector<Tp, 3> s) : description_set{s[0], s[1], s[2]} { dcm_from_parameters_(); }

  template <typename Tp2, size_t n2_items>
  mrp(description_set<Tp2, n2_items> * set) : description_set(set) {}

  // dke (function)
  // Returns a new 3x3 matrix (type Tp) that maps angular velocity onto modified
  // rodriguez rates.
  ::matrix<Tp, 3, 3> dke() override {
    Tp norm2 = pow(get_(0), 2) + pow(get_(1), 2) + pow(get_(2), 2);

    return ::matrix<Tp, 3, 3>{
        1. - norm2 + 2. * pow(get_(0), 2), 2. * (get_(0) * get_(1) - get_(2)), 2. * (get_(0) * get_(2) + get_(1)), 
        2. * (get_(0) * get_(1) + get_(2)), 1. - norm2 + 2. * pow(get_(1), 2), 2. * (get_(1) * get_(2) - get_(0)), 
        2. * (get_(0) * get_(2) - get_(1)), 2. * (get_(1) * get_(2) + get_(0)), 1. - norm2 + 2. * pow(get_(2), 2)
    } * 0.25;  
  }

  // reverse (function)
  // Returns the reverse rotation represented as a MRP.
  mrp<Tp> reverse() { return mrp<Tp>(items_ * -1); }


  // -------------------- Classical Rodriguez Addition/Subtraction. --------------------
  mrp<Tp> operator+ (mrp<Tp> rhs) {
    Tp lhs_norm2 = pow(get_(0), 2) + pow(get_(1), 2) + pow(get_(2), 2);
    Tp rhs_norm2 = pow(rhs[0], 2) + pow(rhs[1], 2) + pow(rhs[2], 2);

    vector<Tp, 3> rhs_vec { rhs[0], rhs[1], rhs[2] };

    vector<Tp, 3> result(
        (items_ * (1 - rhs_norm2) + rhs_vec * (1 - lhs_norm2) - cross(items_, rhs_vec) * 2.) /
        (1 + rhs_norm2 * lhs_norm2 - 2 * rhs_vec.inner(items_)));

    return mrp<Tp>{result[0], result[1], result[2]}; 
    // TODO: add vector to constructor of description sets.
    // TODO: also add a PRV object: inheriting from vector, description_set
  }

  mrp<Tp> operator- (mrp<Tp> rhs) {
    Tp lhs_norm2 = pow(get_(0), 2) + pow(get_(1), 2) + pow(get_(2), 2);
    Tp rhs_norm2 = pow(rhs[0], 2) + pow(rhs[1], 2) + pow(rhs[2], 2);

    vector<Tp, 3> rhs_vec{rhs[0], rhs[1], rhs[2]};

    vector<Tp, 3> result(
        (items_ * (1 - rhs_norm2) - rhs_vec * (1 - lhs_norm2) + cross(items_, rhs_vec) * 2.) /
        (1 + rhs_norm2 * lhs_norm2 + 2 * rhs_vec.inner(items_)));

    return mrp<Tp>{result[0], result[1], result[2]};
  }

  mrp<Tp> * operator+= (mrp<Tp> rhs) {
    Tp lhs_norm2 = pow(get_(0), 2) + pow(get_(1), 2) + pow(get_(2), 2);
    Tp rhs_norm2 = pow(rhs[0], 2) + pow(rhs[1], 2) + pow(rhs[2], 2);

    vector<Tp, 3> rhs_vec{rhs[0], rhs[1], rhs[2]};

    vector<Tp, 3> result(
        (items_ * (1 - rhs_norm2) + rhs_vec * (1 - lhs_norm2) - 2 * items_.cross(rhs_vec)) /
        (1 + rhs_norm2 * lhs_norm2 - 2 * rhs_vec.inner(items_)));

    set_([result[0], result[1], result[2]]);
    return this;
  }

  mrp<Tp> * operator-= (mrp<Tp> q) {
    Tp lhs_norm2 = pow(get_(0), 2) + pow(get_(1), 2) + pow(get_(2), 2);
    Tp rhs_norm2 = pow(q[0], 2) + pow(q[1], 2) + pow(q[2], 2);

    vector<Tp, 3> rhs_vec{rhs[0], rhs[1], rhs[2]};

    vector<Tp, 3> result(
        (items_ * (1 - rhs_norm2) - rhs_vec * (1 - lhs_norm2) + 2 * items_.cross(rhs_vec)) /
        (1 + rhs_norm2 * lhs_norm2 + 2 * rhs_vec.inner(items_)));

    set_([ result[0], result[1], result[2] ]);
    return this;
  }


  // -------------------- Comparison. --------------------
  bool operator==(mrp<Tp> q) { return matrix() == q.matrix(); }

 private:
  // -------------------- Description set virtual overrides. --------------------

  // parameters_from_dcm_ (function)
  // Computes the modified rodriguez components from directional cosine matrix.
  void parameters_from_dcm_() override {
    Tp quaternion[4];
    shepherds_rule(matrix_, quaternion);

    Tp sigma[3] {
        quaternion[1] / quaternion[0],
        quaternion[2] / quaternion[0],
        quaternion[3] / quaternion[0]
    };

    set_(sigma);
    update_dcm_ = false; // dcm & parameters are now consistent
  }

  // dcm_from_parameters_ (function)
  // Computes the directional cosine matrix from modified rodriguez components.
  void dcm_from_parameters_() override {
    Tp norm2 = pow(get_(0), 2) + pow(get_(1), 2) + pow(get_(2), 2);

    matrix_ = EYE<Tp, 3>() + ((tilde(items_) * tilde(items_) * 8) - (tilde(items_) * (4 * (1 - norm2)))) /
        pow(1 + norm2, 2);

    update_dcm_ = false; // dcm & parameters are now consistent
  }
};

} // namespace attitude
#endif  // ATT_MRODR_H_
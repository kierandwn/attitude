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
// quaternion.h
// Contains concrete definition for quaternion (euler parameter) attitude 
// description parameter set.
//
// File namespace attitude::
//
#ifndef ATT_QUATERNION_H_
#define ATT_QUATERNION_H_

#include <initializer_list>

#include "attitude/base.h"
#include "attitude/matrix.h"
#include "attitude/shepherds.h"

namespace attitude {


// quaternion (class)
// Quaternion/ Euler Parameters (attitude representation).
// Decribes direct addition & subtraction methods for quaternion,
// defines mapping to/from DCM, and angular velocity to quaternion rates.
template<typename Tp>
class quaternion : public virtual description_set<Tp, 4> {
 public:
  quaternion() {}
  
  quaternion(Tp q0, Tp q1, Tp q2, Tp q3) : description_set{q0, q1, q2, q3} 
  { dcm_from_parameters_(); }
  
  quaternion(::mn_matrix<Tp, 3, 3> R) : description_set(R) 
  { parameters_from_dcm_(); }

  template <typename Tp2, size_t n2_items>
  quaternion(description_set<Tp2, n2_items> set) : description_set(set) {}

  // dke (function)
  // Returns a new 4x4 matrix (type Tp) that maps modified angular velocity
  // vector onto euler rates. (mod. angular velocity vec [0, w1, w2, w3]^T).
  ::mn_matrix<Tp, 4, 4> dke() override {
    return ::mn_matrix<Tp, 4, 4>{
        get_(0), -1 * get_(1), -1 * get_(2), -1 * get_(3),
        get_(1),      get_(0), -1 * get_(3),      get_(2),
        get_(2),      get_(3),      get_(0), -1 * get_(1),
        get_(3), -1 * get_(2),      get_(1),      get_(0)
    } * 0.5;
  }

  // -------------------- Quaternion Addition/Subtraction. --------------------
  quaternion<Tp> operator+ (quaternion<Tp> q) {
    return quaternion<Tp> {
        q[0] * get_(0) + -1. * q[1] * get_(1) + -1. * q[2] * get_(2) + -1. * q[3] * get_(3),
        q[1] * get_(0) +       q[0] * get_(1) +       q[3] * get_(2) +       q[2] * get_(3),
        q[2] * get_(0) + -1. * q[3] * get_(1) +       q[0] * get_(2) +       q[1] * get_(3),
        q[3] * get_(0) +       q[2] * get_(1) + -1. * q[1] * get_(2) +       q[0] * get_(3)
    };
  }

  quaternion<Tp> operator- (quaternion<Tp> q)
  {
    return quaternion<Tp> {
              q[0] * get_(0) +       q[1] * get_(1) +       q[2] * get_(2) +       q[3] * get_(3),
        -1. * q[1] * get_(0) +       q[0] * get_(1) + -1. * q[3] * get_(2) +       q[2] * get_(3),
        -1. * q[2] * get_(0) +       q[3] * get_(1) +       q[0] * get_(2) + -1. * q[1] * get_(3),
        -1. * q[3] * get_(0) + -1. * q[2] * get_(1) +       q[1] * get_(2) +       q[0] * get_(3)
    };
  }

  quaternion<Tp> * operator+= (quaternion<Tp> q)
    {
      set_([
          q[0] * get_(0) + -1. * q[1] * get_(1) + -1. * q[2] * get_(2) + -1. * q[3] * get_(3),
          q[1] * get_(0) +       q[0] * get_(1) +       q[3] * get_(2) +       q[2] * get_(3),
          q[2] * get_(0) + -1. * q[3] * get_(1) +       q[0] * get_(2) +       q[1] * get_(3),
          q[3] * get_(0) +       q[2] * get_(1) + -1. * q[1] * get_(2) +       q[0] * get_(3)
      ]);
      return this;
    }

    quaternion<Tp> * operator-= (quaternion<Tp> q)
    {
      set_([
                  q[0] * get_(0) +       q[1] * get_(1) +       q[2] * get_(2) +       q[3] * get_(3),
            -1. * q[1] * get_(0) +       q[0] * get_(1) + -1. * q[3] * get_(2) +       q[2] * get_(3),
            -1. * q[2] * get_(0) +       q[3] * get_(1) +       q[0] * get_(2) + -1. * q[1] * get_(3),
            -1. * q[3] * get_(0) + -1. * q[2] * get_(1) +       q[1] * get_(2) +       q[0] * get_(3)
      ]);
      return this;
    }

    // -------------------- Comparison. --------------------
    bool operator==(quaternion<Tp> q) { return matrix() == q.matrix(); }
    bool operator==(::mn_matrix<Tp, 3, 3> rhs) { return matrix() == rhs; }

 private:
  // -------------------- Description set virtual overrides. --------------------

  // parameters_from_dcm_ (function)
  // Computes the quaternion components from directional cosine matrix (shepherd's method).
  void parameters_from_dcm_() override {
    Tp quaternion[4];
    shepherds_rule(matrix_, quaternion);

    set_(quaternion);
    update_dcm_ = false; // dcm & parameters are now consistent
  }

  // dcm_from_parameters_ (function)
  // Computes the directional cosine matrix from quaternion components.
  void dcm_from_parameters_() override {
    // Row 1
    matrix_[0][0] = pow(get_(0), 2) + pow(get_(1), 2) - pow(get_(2), 2) - pow(get_(3), 2.);
    matrix_[0][1] = 2 * (get_(1) * get_(2) + get_(0) * get_(3));
    matrix_[0][2] = 2. * (get_(1) * get_(3) - get_(0) * get_(2));
    // Row 2
    matrix_[1][0] = 2. * (get_(1) * get_(2) - get_(0) * get_(3));
    matrix_[1][1] = pow(get_(0), 2) - pow(get_(1), 2) + pow(get_(2), 2) - pow(get_(3), 2);
    matrix_[1][2] = 2. * (get_(2) * get_(3) + get_(0) * get_(1));
    // Row 3
    matrix_[2][0] = 2. * (get_(1) * get_(3) + get_(0) * get_(2));
    matrix_[2][1] = 2. * (get_(2) * get_(3) - get_(0) * get_(1));
    matrix_[2][2] = pow(get_(0), 2) - pow(get_(1), 2) - pow(get_(2), 2) + pow(get_(3), 2);
    update_dcm_ = false; // dcm & parameters are now consistent
  }
};


} // namespace attitude
#endif  // ATT_QUATERNION_H_
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

#include <cmath>

#include "base.h"
#include "matrix.h"

namespace attitude {


// Returns the index of the maximum value in array of type ElementType
//  
template<typename Tp>
int argmax(Tp * arr, int length)
{
    int max_index = 0;

    for (int k = 1; k < length; k++)
        if (arr[k] > arr[max_index]){ max_index = k; }

    return max_index;
}

template<typename Tp>
class quaternion : public virtual description_set<Tp, 4> {
 public:
  quaternion() {}
  
  quaternion(Tp q0, Tp q1, Tp q2, Tp q3) : description_set{q0, q1, q2, q3} 
  { dcm_from_parameters_(); }
  
  quaternion(::matrix<Tp, 3, 3> R) : description_set(R) 
  { parameters_from_dcm_(); }

  template <typename Tp2, size_t n2_items>
  quaternion(description_set<Tp2, n2_items> set) : description_set(set) {}

  ::matrix<Tp, 4, 4> dke() override {
    return ::matrix<Tp, 4, 4>{
        get_(0), -1 * get_(1), -1 * get_(2), -1 * get_(3),
        get_(1),      get_(0), -1 * get_(3),      get_(2),
        get_(2),      get_(3),      get_(0), -1 * get_(1),
        get_(3), -1 * get_(2),      get_(1),      get_(0)
    } * 0.5;
  }

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

 private:
  void parameters_from_dcm_() override {
    // Shepherd's Rule
    Tp _4q0q0 = 1 + matrix_[0][0] + matrix_[1][1] + matrix_[2][2];
    Tp _4q1q1 = 1 + matrix_[0][0] - matrix_[1][1] - matrix_[2][2];
    Tp _4q2q2 = 1 - matrix_[0][0] + matrix_[1][1] - matrix_[2][2];
    Tp _4q3q3 = 1 - matrix_[0][0] - matrix_[1][1] + matrix_[2][2];

    Tp vals[4]{abs(_4q0q0), abs(_4q1q1), abs(_4q2q2), abs(_4q3q3)};
    int i = argmax(vals, 4);

    switch (i) { 
     case 0: {
      Tp _4q0q1 = matrix_[1][2] - matrix_[2][1];
      Tp _4q0q2 = matrix_[2][0] - matrix_[0][2];
      Tp _4q0q3 = matrix_[0][1] - matrix_[1][0];
      Tp _q0 = sqrt(0.25 * _4q0q0);

      set_(0, _q0);
      set_(1, 0.25 * _4q0q1 / _q0);
      set_(2, 0.25 * _4q0q2 / _q0);
      set_(3, 0.25 * _4q0q3 / _q0);
     } break;

     case 1: {
      Tp _4q0q1 = matrix_[1][2] - matrix_[2][1];
      Tp _4q1q2 = matrix_[0][1] + matrix_[1][0];
      Tp _4q1q3 = matrix_[2][0] + matrix_[0][2];
      Tp _q1 = sqrt(0.25 * _4q1q1);

      set_(0, 0.25 * _4q0q1 / _q1);
      set_(1, _q1);
      set_(2, 0.25 * _4q1q2 / _q1);
      set_(3, 0.25 * _4q1q3 / _q1);
     } break;

     case 2: {
      Tp _4q0q2 = matrix_[2][0] - matrix_[0][2];
      Tp _4q1q2 = matrix_[0][1] + matrix_[1][0];
      Tp _4q2q3 = matrix_[1][2] + matrix_[2][1];
      Tp _q2 = sqrt(0.25 * _4q2q2);

      set_(0, 0.25 * _4q0q2 / _q2);
      set_(1, 0.25 * _4q1q2 / _q2);
      set_(2, _q2);
      set_(3, 0.25 * _4q2q3 / _q2);
     } break;

     case 3: {
      Tp _4q0q3 = matrix_[0][1] - matrix_[1][0];
      Tp _4q1q3 = matrix_[2][0] + matrix_[0][2];
      Tp _4q2q3 = matrix_[1][2] + matrix_[2][1];
      Tp _q3 = sqrt(0.25 * _4q3q3);

      set_(0, 0.25 * _4q0q3 / _q3);
      set_(1, 0.25 * _4q1q3 / _q3);
      set_(2, 0.25 * _4q2q3 / _q3);
      set_(3, _q3);
     } break;
    }
    update_dcm_ = false;
  }

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
    update_dcm_ = false;
  }
};


} // namespace attitude

template<typename T>
void display(attitude::quaternion<T> q)
{
    std::cout
        << "[ "
        << q[0] << ","
        << q[1] << ","
        << q[2] << ","
        << q[3]
        << " ]"
        << std::endl;
}

#endif  // ATT_QUATERNION_H_
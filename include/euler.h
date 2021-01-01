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
// euler.h
// Contains concrete definition for euler angles attitude description parameter 
// set.
//
// File namespace attitude::
//
#ifndef ATT_EULER_H_
#define ATT_EULER_H_

#include <cmath>
#include <initializer_list>

#include "base.h"
#include "matrix.h"

namespace attitude {


// euler (class)
// Euler Angles (attitude representation).
// Decribes direct addition & subtraction methods for euler angles, of any
// order, and defines mapping to/from DCM, and angular velocity to MRP rates.
template <typename Tp>
class euler : public virtual description_set<Tp, 3> 
{
private:
  const uint16_t ijk_;

public:
  euler(uint16_t ijk=123) : ijk_(ijk) {}

  euler(Tp t1, Tp t2, Tp t3, uint16_t ijk=123)
      : description_set{ t1, t2, t3 },
        ijk_(ijk) 
  { dcm_from_parameters_(); }

  template <typename Tp2, size_t n2_items>
  euler(description_set<Tp2, n2_items> set) : description_set(set) 
  { parameters_from_dcm_(); }

  euler(::matrix<Tp, 3, 3> R, uint16_t ijk=123) 
      : description_set(R),
        ijk_(ijk)
  { parameters_from_dcm_(); }

  // dke (function)
  // Returns a new 3x3 matrix (type Tp) that maps angular velocity
  // onto euler rates.
  ::matrix<Tp, 3, 3> dke() override {
    switch (order()) {
      case 123:
        return
        {
                 Tp(cos(get_(2)) / cos(get_(1))), Tp(-1 * sin(get_(2)) / cos(get_(1))), Tp(0.), 
                                Tp(sin(get_(2))),                     Tp(cos(get_(2))), Tp(0.),
            Tp(-1 * cos(get_(2)) * tan(get_(1))),      Tp(sin(get_(2)) * tan(get_(1))), Tp(1.)
        };
      
      case 321:
        return
        {
            Tp(0.),      Tp(sin(get_(2)) / cos(get_(1))), Tp(cos(get_(2)) / cos(get_(1))),
            Tp(0.),                     Tp(cos(get_(2))),           Tp(-1 * sin(get_(2))),
            Tp(1.), Tp(-1 * sin(get_(2)) * tan(get_(1))), Tp(cos(get_(2)) * tan(get_(1)))
        };
      
      default:
        return EYE<Tp, 3>(); // identity
    }
  }

  // order (function)
  // Returns a eulr angle order (16bit uint).
  uint16_t order() { return ijk_; }


  // -------------------- Euler Addition/Subtraction. --------------------
  // No direct optimisation - has to go through rotation matrices.

  euler<Tp> operator+ (euler<Tp> theta) {
    return euler<Tp>(matrix() * theta.matrix(), order());
  }

  euler<Tp> operator- (euler<Tp> theta) {
    return euler<Tp>(matrix() / theta.matrix(), order());
  }

  euler<Tp> * operator+= (euler<Tp> theta) {
    matrix_ *= theta.matrix();
    parameters_from_dcm_();
    return this;
  }

  euler<Tp> * operator-= (euler<Tp> theta) {
    matrix_ /= theta.matrix();
    parameters_from_dcm_();
    return this;
  }

  // -------------------- Comparison. --------------------
  bool operator== (euler<Tp> theta) { return matrix() == theta.matrix(); }

 private:
  // -------------------- Description set virtual overrides. --------------------
  
  // parameters_from_dcm_ (function)
  // Computes the euler angle parameters from directional cosine matrix.
  void parameters_from_dcm_ () override {
    switch (order()) {
      case 121: {
        items_[0] = atan2(-matrix_[0][1], matrix_[0][2]);
        items_[1] = acos(matrix_[0][0]);
        items_[2] = atan2(matrix_[1][0], matrix_[2][0]);
      } break;
      case 123: {
        items_[0] = atan2(-matrix_[2][1], matrix_[2][2]);
        items_[1] = asin(matrix_[2][0]);
        items_[2] = atan2(-matrix_[1][0], matrix_[0][0]);
      } break;
      case 131: {
        items_[0] = atan2(matrix_[0][2], matrix_[0][1]);
        items_[1] = acos(matrix_[0][0]);
        items_[2] = atan2(matrix_[2][0], -matrix_[1][0]);
      } break;
      case 132: {
        items_[0] = atan2(matrix_[1][2], matrix_[1][1]);
        items_[1] = asin(-matrix_[1][0]);
        items_[2] = atan2(matrix_[2][0], matrix_[0][0]);
      } break;
      case 212: {
        items_[0] = atan2(matrix_[1][0], matrix_[1][2]);
        items_[1] = acos(matrix_[1][1]);
        items_[2] = atan2(matrix_[0][1], -matrix_[2][1]);
      } break;
      case 213: {
        items_[0] = atan2(matrix_[2][0], matrix_[2][2]);
        items_[1] = asin(-matrix_[2][1]);
        items_[2] = atan2(matrix_[0][1], matrix_[1][1]);
      } break;
      case 231: {
        items_[0] = atan2(-matrix_[0][2], matrix_[0][0]);
        items_[1] = asin(matrix_[0][1]);
        items_[2] = atan2(-matrix_[2][1], matrix_[1][1]);
      } break;
      case 232: {
        items_[0] = atan2(matrix_[1][2], -matrix_[1][0]);
        items_[1] = acos(matrix_[1][1]);
        items_[2] = atan2(matrix_[2][1], matrix_[0][1]);
      } break;
      case 312: {
        items_[0] = atan2(-matrix_[1][0], matrix_[1][1]);
        items_[1] = asin(matrix_[1][2]);
        items_[2] = atan2(-matrix_[0][2], matrix_[2][2]);
      } break;
      case 313: {
        items_[0] = atan2(matrix_[2][0], -matrix_[2][1]);
        items_[1] = acos(matrix_[2][2]);
        items_[2] = atan2(matrix_[0][2], matrix_[1][2]);
      } break;
      case 321: {
        items_[0] = atan2(matrix_[0][1], matrix_[0][0]);
        items_[1] = asin(-matrix_[0][2]);
        items_[2] = atan2(matrix_[1][2], matrix_[2][2]);
      } break;
      case 323: {
        items_[0] = atan2(matrix_[2][0], -matrix_[2][1]);
        items_[1] = acos(matrix_[2][2]);
        items_[2] = atan2(matrix_[0][2], matrix_[1][2]);
      } break;
    }
    update_dcm_ = false; // dcm & parameters are now consistent
  }

  // dcm_from_parameters_ (function)
  // Computes the directional cosine matrix from euler angle parameters.
  void dcm_from_parameters_ () override {
    uint8_t i = ijk_ / 100;
    uint8_t j = (ijk_ - (i * 100)) / 10;
    uint8_t k = (ijk_ - (i * 100) - (j * 10));

    matrix_  = dcm::AXIS(k, items_[2]) *
               dcm::AXIS(j, items_[1]) *
               dcm::AXIS(i, items_[0]);
    update_dcm_ = false; // dcm & parameters are now consistent
  }
};

} // namespace attitude
#endif // ATT_EULER_H_
#ifndef EULER_H
#define EULER_H

#include <cmath>

#include "dcm.h"

namespace attitude
{
    template <typename T>
    class euler_
    {
    private:
        dcm_<T> _R;
        const uint16_t _ijk = 123;

        T _theta[3] = { 0 };
        T _get(int i) { return _theta[i]; }

    public:
        euler_() : _theta{ T(0.), T(0.), T(0.) },
                   _R(to_dcm(*this))
        {}

        euler_(T t1, T t2, T t3) : _theta{ t1, t2, t3 },
                                   _R(to_dcm(*this))
        {}

        euler_(T t1, T t2, T t3, char ijk[4]) : _theta{ t1, t2, t3 },
                                                _ijk(ijk),
                                                _R(to_dcm(*this))
        {}

        euler_(dcm_<T> R, uint16_t ijk) : // _theta{ extract_euler_from_dcm(R, ijk) },
                                          _ijk(ijk),
                                          _R(R)
        {
            extract_euler_from_dcm(R, ijk);
        }

        uint16_t order() { return _ijk; }
        dcm_<T> dcm() { return _R; }
        euler_<T> reverse() { return euler_<T>( dcm().reverse(), order() ); }
        
        euler_<T> operator+(euler_<T> theta)
        {
            dcm_<T> R = dcm() + theta.dcm();
            return euler_<T>(R, order());
        }

        euler_<T> operator-(euler_<T> theta)
        {
            dcm_<T> R = dcm() - theta.dcm();
            return euler_<T>(R, order());
        }

        euler_<T>* operator+=(euler_<T> theta)
        {
            _R += theta.dcm();
            return this;
        }

        T operator[](int i) { return _get(i); }

        // Comparison
        bool operator==(euler_<T> theta) { return dcm() == theta.dcm(); }
        bool operator==(dcm_<T> R) { return R == dcm(); }

        void extract_euler_from_dcm(dcm_<T> R, uint16_t ijk)
        {
            if (ijk == 121)
            {
                _theta[0] = atan2(-R[0][1], R[0][2]);
                _theta[1] = acos(R[0][0]);
                _theta[2] = atan2(R[1][0], R[2][0]);
            }
            else if (ijk == 123)
            {
                _theta[0] = atan2(-R[2][1], R[2][2]);
                _theta[1] = asin(R[2][0]);
                _theta[2] = atan2(-R[1][0], R[0][0]);
            }
            else if (ijk == 131)
            {
                _theta[0] = atan2(R[0][2], R[0][1]);
                _theta[1] = acos(R[0][0]);
                _theta[2] = atan2(R[2][0], -R[1][0]);
            }
            else if (ijk == 132)
            {
                _theta[0] = atan2(R[1][2], R[1][1]);
                _theta[1] = asin(-R[1][0]);
                _theta[2] = atan2(R[2][0], R[0][0]);
            }
            else if (ijk == 212)
            {
                _theta[0] = atan2(R[1][0], R[1][2]);
                _theta[1] = acos(R[1][1]);
                _theta[2] = atan2(R[0][1], -R[2][1]);
            }
            else if (ijk == 213)
            {
                _theta[0] = atan2(R[2][0], R[2][2]);
                _theta[1] = asin(-R[2][1]);
                _theta[2] = atan2(R[0][1], R[1][1]);
            }
            else if (ijk == 231)
            {
                _theta[0] = atan2(-R[0][2], R[0][0]);
                _theta[1] = asin(R[0][1]);
                _theta[2] = atan2(-R[2][1], R[1][1]);
            }
            else if (ijk == 232)
            {
                _theta[0] = atan2(R[1][2], -R[1][0]);
                _theta[1] = acos(R[1][1]);
                _theta[2] = atan2(R[2][1], R[0][1]);
            }
            else if (ijk == 312)
            {
                _theta[0] = atan2(-R[1][0], R[1][1]);
                _theta[1] = asin(R[1][2]);
                _theta[2] = atan2(-R[0][2], R[2][2]);
            }
            else if (ijk == 313)
            {
                _theta[0] = atan2(R[2][0], -R[2][1]);
                _theta[1] = acos(R[2][2]);
                _theta[2] = atan2(R[0][2], R[1][2]);
            }
            else if (ijk == 321)
            {
                _theta[0] = atan2(R[0][1], R[0][0]);
                _theta[1] = asin(-R[0][2]);
                _theta[2] = atan2(R[1][2], R[2][2]);
            }
            else if (ijk == 323)
            {
                _theta[0] = atan2(R[2][0], -R[2][1]);
                _theta[1] = acos(R[2][2]);
                _theta[2] = atan2(R[0][2], R[1][2]);
            }
        }
    };

    // Conversions
    template<typename T>
    dcm_<T> to_dcm(euler_<T> theta)
    {
        uint16_t ijk = theta.order();

        uint8_t i = ijk / 100;
        uint8_t j = (ijk - (i * 100)) / 10;
        uint8_t k = (ijk - (i * 100) - (j * 10));

        dcm_<T> R = AXIS(i, theta[0]);
        R += AXIS(j, theta[1]);
        R += AXIS(k, theta[2]);
        return R;
    }
}

template<typename T>
void display(attitude::euler_<T> theta)
{
    std::cout
        << "[ "
        << theta[0] << ","
        << theta[1] << ","
        << theta[2]
        << " ]"
        << std::endl;
}

#endif // EULER_H
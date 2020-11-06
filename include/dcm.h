#ifndef DCM_H
#define DCM_H

#include "common.h"

// Compiler instructions - option for lightwight trig implementation?
#include <array>
#include <cmath>

namespace attitude
{
    template<typename T>
    class dcm_
    {
    private:
        std::array<T, 9> _R;

    public:
        dcm_() : _R{ 1, 0, 0,
                     0, 1, 0,
                     0, 0, 1 }
        {}

        dcm_(std::array<T, 9> R) : _R(R)
        {}

        dcm_<T> get_reverse()
        {
            return dcm_<T>(math::matrix_transpose(_R));
        }

        std::array<T, 9> get_matrix()
        {
            return _R;
        }

        dcm_<T> operator+(dcm_<T> R)
        {
            return dcm_<T>(math::matrix_multiply(_R, R.get_matrix()));
        }

        dcm_<T> operator-(dcm_<T> R)
        {
            dcm_<T> R3 = operator+(R.get_reverse());
            return R3;
        }

        bool operator==(dcm_<T> R)
        {
            std::array<T, 9> diff = math::array_diff(_R, R.get_matrix());
            return math::matrix_is_value(diff, T(0));
        }

        bool operator==(std::array<T, 9> R)
        {
            return _R == R;
        }

        std::array<T, 3> operator[](int i)
        {
            int lowest_index = i * 3;
            return { _R[lowest_index], _R[lowest_index + 1], _R[lowest_index + 2] };
        }
    };

    template<typename T>
    dcm_<T> ZERO()
    {
        return dcm_<T>();
    }

    template<typename T>
    dcm_<T> R1(T theta)
    {
        return dcm_<T>({
            T(1.), T(0.), T(0.),
            T(0.),  cos(theta), sin(theta),
            T(0.), -sin(theta), cos(theta)
            });
    }

    template<typename T>
    dcm_<T> R2(T theta)
    {
        return dcm_<T>({
            cos(theta), T(0.), -sin(theta),
            T(0.), T(1.), T(0.),
            sin(theta), T(0.),  cos(theta)
            });
    }

    template<typename T>
    dcm_<T> R3(T theta)
    {
        return dcm_<T>({
            cos(theta), sin(theta), T(0.),
           -sin(theta), cos(theta), T(0.),
            T(0.), T(0.), T(1.)
            });
    }
}

template<typename T>
void display(attitude::dcm_<T> R)
{
    std::array<T, 9> _R = R.get_matrix();

    for (int i = 0; i < 3; i++)
    {
        int lowest_index = i * 3;
        T row[3]{ _R[lowest_index], _R[lowest_index + 1], _R[lowest_index + 2] };

        std::cout
            << "[ "
            << row[0] << ","
            << row[1] << ","
            << row[2]
            << " ]"
            << std::endl;
    }
}

#endif // DCM_H
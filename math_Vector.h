#ifndef MATH_VECTOR_H
#define MATH_VECTOR_H

#include <cmath>

namespace math
{
    struct Vector
    {
        double x, y, z;

        friend Vector operator+ (const Vector& u, const Vector& v)
        {
            return Vector{ u.x + v.x, u.y + v.y, u.z + v.z };
        }

        friend Vector operator- (const Vector& u, const Vector& v)
        {
            return Vector{ u.x - v.x, u.y - v.y, u.z - v.z };
        }

        friend Vector operator* (double w, const Vector& v)
        {
            return Vector{ w * v.x, w * v.y, w * v.z };
        }

        friend double dot(const Vector& u, const Vector& v)
        {
            return u.x * v.x + u.y * v.y + u.z * v.z;
        }

        friend double norm(const Vector& v)
        {
            return std::sqrt(dot(v, v));
        }

        friend Vector cross(const Vector& u, const Vector& v)
        {
            return Vector{
                u.y * v.z - u.z * v.y,
                u.z * v.x - u.x * v.z,
                u.x * v.y - u.y * v.x };
        }
    };
} // namespace math

#endif // MATH_VECTOR_H

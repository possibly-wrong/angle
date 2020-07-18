#ifndef MATH_ANGLE_H
#define MATH_ANGLE_H

#include "math_Rational.h"
#include "math_Vector.h"

namespace math
{
    // Return sin(x) Taylor series to order x^(2n+1).
    Rational sin(const Rational& x, Unsigned n)
    {
        Rational x2 = -x * x;
        Rational term = x;
        Rational series = term;
        for (Unsigned k = 1; k <= n; ++k)
        {
            term *= (x2 / (2 * k * (2 * k + 1)));
            series += term;
        }
        return series;
    }

    // Return cos(x) Taylor series to order x^(2n).
    Rational cos(const Rational& x, Unsigned n)
    {
        Rational x2 = -x * x;
        Rational term = 1;
        Rational series = term;
        for (Unsigned k = 1; k <= n; ++k)
        {
            term *= (x2 / (2 * k * (2 * k - 1)));
            series += term;
        }
        return series;
    }

    // Return 1/sqrt(x) with ~200-bit accuracy.
    Rational inv_sqrt(const Rational& x)
    {
        // Iterate Newton's method with f(y) = y^-2 - x.
        Rational y = std::sqrt((1 / x).to_double());
        for (int i = 0; i < 2; ++i)
        {
            y *= ((3 - x * y * y) / 2);
        }
        return y;
    }

    const Rational PI = Rational(
        "3.14159265358979323846264338327950288419716939937510582097494");

    // Return atan(sqrt(x2)) with ~200-bit accuracy.
    Rational atan_sqrt(Rational x2)
    {
        // Assume 0 < x2 <= 1 using atan(x) == pi/2 - atan(1/x).
        if (x2 == 0)
        {
            return 0;
        }
        bool invert = (x2 > 1);
        if (invert)
        {
            x2 = 1 / x2;
        }

        // Compute Euler's series for atan(x) to n=200.
        Rational y = x2 / (1 + x2);
        Rational term = y;
        Rational series = term;
        for (Unsigned n = 1; n <= 200; ++n)
        {
            term *= (y * (2 * n) / (2 * n + 1));
            series += term;
        }
        y = inv_sqrt(x2) * series;

        // Convert back if we used reciprocal identity.
        if (invert)
        {
            y = PI / 2 - y;
        }
        return y;
    }

    // Return "exact" angle(u, v) with ~200-bit accuracy.
    Rational angle_exact(const Vector& u, const Vector& v)
    {
        Rational ux = u.x, uy = u.y, uz = u.z;
        Rational vx = v.x, vy = v.y, vz = v.z;
        Rational d = ux * vx + uy * vy + uz * vz; // dot(u, v)
        Rational
            cx = uy * vz - uz * vy, // cross(u, v)
            cy = uz * vx - ux * vz,
            cz = ux * vy - uy * vx;
        if (d == 0)
        {
            return PI / 2;
        }
        Rational angle = atan_sqrt((cx * cx + cy * cy + cz * cz) / (d * d));
        if (d < 0)
        {
            angle = PI - angle;
        }
        return angle;
    }

    // Return angle between vectors using same "cross/dot product" formula.
    double angle(const Vector& u, const Vector& v)
    {
        return std::atan2(norm(cross(u, v)), dot(u, v));
    }

    // Return angle between vectors using law of cosines.
    double angle0(const Vector& u, const Vector& v)
    {
        double u2 = u.x * u.x + u.y * u.y + u.z * u.z;
        double v2 = v.x * v.x + v.y * v.y + v.z * v.z;
        Vector d = u - v;
        return std::acos(std::min(1.0, std::max(-1.0,
            (u2 + v2 - (d.x * d.x + d.y * d.y + d.z * d.z)) /
            (2 * std::sqrt(u2) * std::sqrt(v2)))));
    }

    // Following are three alternative formulas discussed in:
    // Kahan, W., How Futile are Mindless Assessments of Roundoff in
    // Floating-Point Computation?
    // (http://people.eecs.berkeley.edu/~wkahan/Mindless.pdf)

    double angle1(const Vector& u, const Vector& v)
    {
        return std::acos(std::min(1.0, std::max(-1.0,
            dot(u, v) / (norm(u) * norm(v)))));
    }

    double angle2(const Vector& u, const Vector& v)
    {
        double angle = std::asin(std::min(1.0,
            norm(cross(u, v)) / (norm(u) * norm(v))));
        if (dot(u, v) < 0)
        {
            angle = 3.141592653589793 - angle;
        }
        return angle;
    }

    double angle3(const Vector& u, const Vector& v)
    {
        double nu = norm(u);
        double nv = norm(v);
        return 2 * std::atan2(norm(nv * u - nu * v), norm(nv * u + nu * v));
    }
} // namespace math

#endif // MATH_ANGLE_H

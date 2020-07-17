#include "math_angle.h"
#include <random>
#include <limits>
#include <iostream>
using namespace math;

// Compute and display angle between vectors u and v.
bool show(int region, double log10_offset, int rotate,
    const Vector& u, const Vector& v)
{
    Rational truth = angle_exact(u, v);
    if (std::isfinite(log10_offset) && truth == 0)
    {
        return true;
    }
    std::cout << region << " " << log10_offset << " " << rotate << std::endl <<
        "    " << Rational(u.x) << std::endl <<
        "    " << Rational(u.y) << std::endl <<
        "    " << Rational(u.z) << std::endl <<
        "    " << Rational(v.x) << std::endl <<
        "    " << Rational(v.y) << std::endl <<
        "    " << Rational(v.z) << std::endl <<
        "    " << truth.to_string(80) << std::endl <<
        "    " << Rational(angle(u, v)) << std::endl <<
        "    " << Rational(angle0(u, v)) << std::endl <<
        "    " << Rational(angle1(u, v)) << std::endl <<
        "    " << Rational(angle2(u, v)) << std::endl <<
        "    " << Rational(angle3(u, v)) << std::endl;
    return false;
}

int main()
{
    std::mt19937 rng;
    rng.seed(5489);
    std::normal_distribution<double> randn(0.0, 1.0);

    // Evaluate a range of magnitudes of "small" angle offsets over
    // [0, 10^(-18:0.125:-0.5)] radians.
    std::vector<double> log10_offsets;
    log10_offsets.push_back(-std::numeric_limits<double>::infinity());
    for (double log10_offset = -18; log10_offset <= -0.5; log10_offset += 0.125)
    {
        log10_offsets.push_back(log10_offset);
    }

    // Evaluate angles in the regions
    // (0+offset, pi/4-offset, pi/2-offset, pi-offset).
    for (int region = 0; region < 4; ++region)
    {
        for (auto&& log10_offset : log10_offsets)
        {
            double offset = std::pow(10.0, log10_offset);
            Rational vx = cos(offset, 20); // ~200 bits
            Rational vy = sin(offset, 20);
            if (region == 1)
            {
                vx += vy; // near pi/4
                vy = vx - 2 * vy;
            }
            else if (region == 2)
            {
                std::swap(vx, vy); // near pi/2
            }
            else if (region == 3)
            {
                vx = -vx; // near pi
            }

            // Evaluate "2D" problem u=(1,0,0) and v=(vx,vy,0).
            Vector u{1.0, 0.0, 0.0};
            Vector v{vx.to_double(), vy.to_double(), 0.0};
            show(region, log10_offset, 0, u, v);

            // Randomly rotate u and v, repeatedly as needed until angle > 0.
            bool searching = true;
            while (searching)
            {
                // Generate uniform random quaternion rotation.
                Rational w = randn(rng),
                    x = randn(rng), y = randn(rng), z = randn(rng);
                Rational inv_mag =
                    inv_sqrt(w*w + x*x + y*y + z*z);
                w *= inv_mag;
                x *= inv_mag;
                y *= inv_mag;
                z *= inv_mag;

                // Rotate vectors and evaluate.
                u.x = (w*w + x*x - y*y - z*z).to_double();
                u.y = (2 * (x*y + w*z)).to_double();
                u.z = (2 * (x*z - w*y)).to_double();
                v.x = ((w*w + x*x - y*y - z*z)*vx + 2 * (x*y - w*z)*vy).to_double();
                v.y = ((w*w - x*x + y*y - z*z)*vy + 2 * (x*y + w*z)*vx).to_double();
                v.z = (2 * ((w*x + y*z)*vy - (w*y - x*z)*vx)).to_double();
                searching = show(region, log10_offset, 1, u, v);
            }
        }
    }
}

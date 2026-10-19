#pragma once

#include <cmath>
#include <ostream>

namespace nbody {

// Minimal 3-vector with value semantics. Kept as a plain aggregate so the
// compiler can keep components in registers inside the force loop.
struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;

    constexpr Vec3() = default;
    constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    constexpr Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    constexpr Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    constexpr Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
};

constexpr Vec3 operator+(Vec3 a, const Vec3& b) { return a += b; }
constexpr Vec3 operator-(Vec3 a, const Vec3& b) { return a -= b; }
constexpr Vec3 operator-(const Vec3& a) { return {-a.x, -a.y, -a.z}; }
constexpr Vec3 operator*(Vec3 a, double s) { return a *= s; }
constexpr Vec3 operator*(double s, Vec3 a) { return a *= s; }

constexpr double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

constexpr Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

inline std::ostream& operator<<(std::ostream& os, const Vec3& v) {
    return os << '(' << v.x << ", " << v.y << ", " << v.z << ')';
}

}  // namespace nbody

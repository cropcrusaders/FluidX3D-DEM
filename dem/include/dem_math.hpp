#pragma once
// DEM Math Utilities - vec3, quaternion, mat3 for particle simulation
// Minimal, self-contained vector math (no Eigen dependency)

#include <cmath>
#include <cstdint>
#include <algorithm>
#include <limits>

namespace dem {

struct vec3 {
    double x, y, z;
    vec3() : x(0), y(0), z(0) {}
    vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    explicit vec3(double s) : x(s), y(s), z(s) {}

    vec3 operator+(const vec3& b) const { return {x+b.x, y+b.y, z+b.z}; }
    vec3 operator-(const vec3& b) const { return {x-b.x, y-b.y, z-b.z}; }
    vec3 operator*(double s) const { return {x*s, y*s, z*s}; }
    vec3 operator/(double s) const { double inv = 1.0/s; return {x*inv, y*inv, z*inv}; }
    vec3 operator-() const { return {-x, -y, -z}; }
    vec3& operator+=(const vec3& b) { x+=b.x; y+=b.y; z+=b.z; return *this; }
    vec3& operator-=(const vec3& b) { x-=b.x; y-=b.y; z-=b.z; return *this; }
    vec3& operator*=(double s) { x*=s; y*=s; z*=s; return *this; }

    double dot(const vec3& b) const { return x*b.x + y*b.y + z*b.z; }
    vec3 cross(const vec3& b) const {
        return {y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x};
    }
    double length2() const { return x*x + y*y + z*z; }
    double length() const { return std::sqrt(length2()); }
    vec3 normalized() const {
        double len = length();
        if (len < 1e-30) return {0,0,0};
        return *this / len;
    }

    double& operator[](int i) { return (&x)[i]; }
    double operator[](int i) const { return (&x)[i]; }
};

inline vec3 operator*(double s, const vec3& v) { return v*s; }
inline double dot(const vec3& a, const vec3& b) { return a.dot(b); }
inline vec3 cross(const vec3& a, const vec3& b) { return a.cross(b); }

struct mat3 {
    double m[3][3] = {};

    mat3() = default;
    static mat3 identity() {
        mat3 r;
        r.m[0][0] = r.m[1][1] = r.m[2][2] = 1.0;
        return r;
    }
    static mat3 diag(double a, double b, double c) {
        mat3 r;
        r.m[0][0] = a; r.m[1][1] = b; r.m[2][2] = c;
        return r;
    }

    vec3 operator*(const vec3& v) const {
        return {
            m[0][0]*v.x + m[0][1]*v.y + m[0][2]*v.z,
            m[1][0]*v.x + m[1][1]*v.y + m[1][2]*v.z,
            m[2][0]*v.x + m[2][1]*v.y + m[2][2]*v.z
        };
    }
    mat3 operator*(const mat3& b) const {
        mat3 r;
        for (int i=0; i<3; i++)
            for (int j=0; j<3; j++)
                for (int k=0; k<3; k++)
                    r.m[i][j] += m[i][k] * b.m[k][j];
        return r;
    }
    mat3 transposed() const {
        mat3 r;
        for (int i=0; i<3; i++)
            for (int j=0; j<3; j++)
                r.m[i][j] = m[j][i];
        return r;
    }
};

struct quat {
    double w, x, y, z;
    quat() : w(1), x(0), y(0), z(0) {}
    quat(double w_, double x_, double y_, double z_) : w(w_), x(x_), y(y_), z(z_) {}

    quat operator*(const quat& q) const {
        return {
            w*q.w - x*q.x - y*q.y - z*q.z,
            w*q.x + x*q.w + y*q.z - z*q.y,
            w*q.y - x*q.z + y*q.w + z*q.x,
            w*q.z + x*q.y - y*q.x + z*q.w
        };
    }
    quat conjugate() const { return {w, -x, -y, -z}; }
    double norm2() const { return w*w + x*x + y*y + z*z; }
    quat normalized() const {
        double n = std::sqrt(norm2());
        if (n < 1e-30) return {1,0,0,0};
        return {w/n, x/n, y/n, z/n};
    }

    vec3 rotate(const vec3& v) const {
        // q * v * q^-1
        quat p(0, v.x, v.y, v.z);
        quat r = (*this) * p * conjugate();
        return {r.x, r.y, r.z};
    }

    mat3 to_matrix() const {
        double xx=x*x, yy=y*y, zz=z*z;
        double xy=x*y, xz=x*z, yz=y*z;
        double wx=w*x, wy=w*y, wz=w*z;
        mat3 r;
        r.m[0][0] = 1-2*(yy+zz); r.m[0][1] = 2*(xy-wz);   r.m[0][2] = 2*(xz+wy);
        r.m[1][0] = 2*(xy+wz);   r.m[1][1] = 1-2*(xx+zz); r.m[1][2] = 2*(yz-wx);
        r.m[2][0] = 2*(xz-wy);   r.m[2][1] = 2*(yz+wx);   r.m[2][2] = 1-2*(xx+yy);
        return r;
    }

    // Integrate angular velocity: q' = q + 0.5*dt*omega*q
    static quat integrate(const quat& q, const vec3& omega, double dt) {
        quat omega_q(0, omega.x, omega.y, omega.z);
        quat dq;
        dq.w = q.w + 0.5*dt*( -omega.x*q.x - omega.y*q.y - omega.z*q.z );
        dq.x = q.x + 0.5*dt*(  omega.x*q.w + omega.z*q.y - omega.y*q.z );
        dq.y = q.y + 0.5*dt*(  omega.y*q.w - omega.z*q.x + omega.x*q.z );
        dq.z = q.z + 0.5*dt*(  omega.z*q.w + omega.y*q.x - omega.x*q.y );
        return dq.normalized();
    }
};

struct AABB {
    vec3 lo, hi;
    AABB() : lo(vec3(std::numeric_limits<double>::max())),
             hi(vec3(-std::numeric_limits<double>::max())) {}
    AABB(const vec3& lo_, const vec3& hi_) : lo(lo_), hi(hi_) {}

    void expand(const vec3& p) {
        lo.x = std::min(lo.x, p.x); lo.y = std::min(lo.y, p.y); lo.z = std::min(lo.z, p.z);
        hi.x = std::max(hi.x, p.x); hi.y = std::max(hi.y, p.y); hi.z = std::max(hi.z, p.z);
    }
    void expand(const AABB& b) {
        lo.x = std::min(lo.x, b.lo.x); lo.y = std::min(lo.y, b.lo.y); lo.z = std::min(lo.z, b.lo.z);
        hi.x = std::max(hi.x, b.hi.x); hi.y = std::max(hi.y, b.hi.y); hi.z = std::max(hi.z, b.hi.z);
    }
    void pad(double r) {
        lo.x -= r; lo.y -= r; lo.z -= r;
        hi.x += r; hi.y += r; hi.z += r;
    }
    vec3 center() const { return (lo + hi) * 0.5; }
    vec3 size() const { return hi - lo; }
    double surface_area() const {
        vec3 s = size();
        return 2.0*(s.x*s.y + s.y*s.z + s.z*s.x);
    }

    bool overlaps(const AABB& b) const {
        return lo.x <= b.hi.x && hi.x >= b.lo.x
            && lo.y <= b.hi.y && hi.y >= b.lo.y
            && lo.z <= b.hi.z && hi.z >= b.lo.z;
    }
    bool contains(const vec3& p) const {
        return p.x >= lo.x && p.x <= hi.x
            && p.y >= lo.y && p.y <= hi.y
            && p.z >= lo.z && p.z <= hi.z;
    }
};

// Random number generator (deterministic, xoshiro256**)
struct RNG {
    uint64_t s[4];

    explicit RNG(uint64_t seed = 12345) {
        // SplitMix64 to initialize state
        s[0] = splitmix(seed);
        s[1] = splitmix(s[0]);
        s[2] = splitmix(s[1]);
        s[3] = splitmix(s[2]);
    }

    uint64_t next() {
        uint64_t result = rotl(s[1] * 5, 7) * 9;
        uint64_t t = s[1] << 17;
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 45);
        return result;
    }

    double uniform() { return (next() >> 11) * 0x1.0p-53; }
    double uniform(double lo, double hi) { return lo + (hi - lo) * uniform(); }

    // Box-Muller normal
    double normal(double mean = 0.0, double stddev = 1.0) {
        double u1 = uniform(), u2 = uniform();
        return mean + stddev * std::sqrt(-2.0 * std::log(u1 + 1e-30)) * std::cos(2.0 * M_PI * u2);
    }

private:
    static uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
    static uint64_t splitmix(uint64_t& state) {
        uint64_t z = (state += 0x9e3779b97f4a7c15);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
        z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
        return z ^ (z >> 31);
    }
};

} // namespace dem

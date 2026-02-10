#pragma once
// Geometry system: triangle mesh (STL/OBJ) with BVH, and analytic primitives

#include "dem_math.hpp"
#include <vector>
#include <string>
#include <memory>
#include <functional>

namespace dem {

// Triangle mesh
struct Triangle {
    vec3 v0, v1, v2;
    vec3 normal() const { return (v1 - v0).cross(v2 - v0).normalized(); }
    AABB bounds() const {
        AABB b;
        b.expand(v0); b.expand(v1); b.expand(v2);
        return b;
    }
};

// Closest point on triangle to point p, returns distance and contact info
struct TriContactResult {
    vec3 closest_point;
    vec3 normal;   // outward normal at contact
    double distance;
    bool valid = false;
};

TriContactResult closest_point_on_triangle(const vec3& p, const Triangle& tri);

// BVH node for triangle mesh
struct BVHNode {
    AABB bounds;
    int left = -1, right = -1;  // child indices (-1 = leaf)
    int tri_start = 0, tri_count = 0;  // leaf triangle range
    bool is_leaf() const { return left == -1 && right == -1; }
};

class TriangleMesh {
public:
    std::vector<Triangle> triangles;
    std::vector<BVHNode> bvh_nodes;
    std::vector<int> tri_indices;  // reordered triangle indices

    bool load_stl(const std::string& path);
    bool load_obj(const std::string& path);
    void build_bvh(int max_leaf_tris = 4);

    // Query: find closest triangle to sphere at pos with radius r
    // Returns penetration depth (positive = overlap) and contact info
    struct SphereQueryResult {
        vec3 contact_point;
        vec3 contact_normal;
        double penetration;
        bool hit = false;
    };
    SphereQueryResult query_sphere(const vec3& pos, double radius) const;

    AABB total_bounds() const;

private:
    int build_bvh_recursive(int start, int count, int depth);
    void query_sphere_recursive(int node_idx, const vec3& pos, double radius,
                                SphereQueryResult& best) const;
};

// Analytic geometry primitives
enum class PrimitiveType { CYLINDER, CONE, BENT_CYLINDER, SPIRAL_INSERT, PLANE };

struct AnalyticPrimitive {
    PrimitiveType type = PrimitiveType::CYLINDER;

    // Common params
    vec3 center = {0, 0, 0};
    vec3 axis = {0, 0, -1};  // tube axis (default: downward)
    double length = 0.3;      // m
    double radius = 0.025;    // m (inner radius)

    // Cone/frustum
    double radius_top = 0.025;
    double radius_bottom = 0.020;

    // Bent cylinder
    double bend_radius = 0.1;  // radius of curvature
    double bend_angle = 0.5;   // radians

    // Spiral insert
    double spiral_pitch = 0.05;
    double spiral_height = 0.02;
    int spiral_turns = 2;

    // Plane
    vec3 plane_normal = {0, 0, 1};
    double plane_offset = 0.0;

    // Query: sphere collision with this primitive
    struct CollisionResult {
        vec3 contact_point;
        vec3 contact_normal;
        double penetration;
        bool hit = false;
    };
    CollisionResult collide_sphere(const vec3& pos, double radius) const;

private:
    CollisionResult collide_cylinder(const vec3& pos, double r) const;
    CollisionResult collide_cone(const vec3& pos, double r) const;
    CollisionResult collide_bent_cylinder(const vec3& pos, double r) const;
    CollisionResult collide_spiral(const vec3& pos, double r) const;
    CollisionResult collide_plane(const vec3& pos, double r) const;
};

// Complete geometry: a collection of mesh + primitives
class Geometry {
public:
    std::vector<TriangleMesh> meshes;
    std::vector<AnalyticPrimitive> primitives;

    struct ContactResult {
        vec3 contact_point;
        vec3 contact_normal;
        double penetration;
        bool hit = false;
    };

    ContactResult query_sphere(const vec3& pos, double radius) const;
    bool load_stl(const std::string& path);
    void add_primitive(const AnalyticPrimitive& prim);
};

} // namespace dem

// Geometry implementation: STL/OBJ loading, BVH, analytic primitives

#include "geometry.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstring>
#include <cmath>

namespace dem {

// ─── Closest point on triangle ───────────────────────────────────────────────

TriContactResult closest_point_on_triangle(const vec3& p, const Triangle& tri) {
    vec3 ab = tri.v1 - tri.v0;
    vec3 ac = tri.v2 - tri.v0;
    vec3 ap = p - tri.v0;

    double d1 = dot(ab, ap);
    double d2 = dot(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) {
        TriContactResult r;
        r.closest_point = tri.v0;
        r.normal = tri.normal();
        r.distance = (p - tri.v0).length();
        r.valid = true;
        return r;
    }

    vec3 bp = p - tri.v1;
    double d3 = dot(ab, bp);
    double d4 = dot(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) {
        TriContactResult r;
        r.closest_point = tri.v1;
        r.normal = tri.normal();
        r.distance = (p - tri.v1).length();
        r.valid = true;
        return r;
    }

    double vc = d1*d4 - d3*d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        double v = d1 / (d1 - d3);
        TriContactResult r;
        r.closest_point = tri.v0 + ab * v;
        r.normal = tri.normal();
        r.distance = (p - r.closest_point).length();
        r.valid = true;
        return r;
    }

    vec3 cp = p - tri.v2;
    double d5 = dot(ab, cp);
    double d6 = dot(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) {
        TriContactResult r;
        r.closest_point = tri.v2;
        r.normal = tri.normal();
        r.distance = (p - tri.v2).length();
        r.valid = true;
        return r;
    }

    double vb = d5*d2 - d1*d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        double w = d2 / (d2 - d6);
        TriContactResult r;
        r.closest_point = tri.v0 + ac * w;
        r.normal = tri.normal();
        r.distance = (p - r.closest_point).length();
        r.valid = true;
        return r;
    }

    double va = d3*d6 - d5*d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        TriContactResult r;
        r.closest_point = tri.v1 + (tri.v2 - tri.v1) * w;
        r.normal = tri.normal();
        r.distance = (p - r.closest_point).length();
        r.valid = true;
        return r;
    }

    double denom = 1.0 / (va + vb + vc);
    double v = vb * denom;
    double w = vc * denom;
    TriContactResult r;
    r.closest_point = tri.v0 + ab * v + ac * w;
    r.normal = tri.normal();
    r.distance = (p - r.closest_point).length();
    r.valid = true;
    return r;
}

// ─── STL loader ──────────────────────────────────────────────────────────────

bool TriangleMesh::load_stl(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) return false;

    // Read header
    char header[80];
    file.read(header, 80);
    if (!file) return false;

    uint32_t num_tris;
    file.read(reinterpret_cast<char*>(&num_tris), 4);
    if (!file || num_tris > 50000000) return false;

    triangles.resize(num_tris);
    for (uint32_t i = 0; i < num_tris; i++) {
        float n[3], v0[3], v1[3], v2[3];
        file.read(reinterpret_cast<char*>(n), 12);
        file.read(reinterpret_cast<char*>(v0), 12);
        file.read(reinterpret_cast<char*>(v1), 12);
        file.read(reinterpret_cast<char*>(v2), 12);
        uint16_t attr;
        file.read(reinterpret_cast<char*>(&attr), 2);
        if (!file) return false;

        triangles[i].v0 = {v0[0], v0[1], v0[2]};
        triangles[i].v1 = {v1[0], v1[1], v1[2]};
        triangles[i].v2 = {v2[0], v2[1], v2[2]};
    }

    build_bvh();
    return true;
}

bool TriangleMesh::load_obj(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) return false;

    std::vector<vec3> vertices;
    std::string line;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string prefix;
        iss >> prefix;

        if (prefix == "v") {
            double x, y, z;
            iss >> x >> y >> z;
            vertices.push_back({x, y, z});
        } else if (prefix == "f") {
            // Parse face (handles v, v/vt, v/vt/vn, v//vn formats)
            std::vector<int> indices;
            std::string token;
            while (iss >> token) {
                int idx = std::stoi(token.substr(0, token.find('/')));
                if (idx < 0) idx = (int)vertices.size() + idx + 1;
                indices.push_back(idx - 1);  // OBJ is 1-indexed
            }
            // Triangulate face (fan)
            for (size_t i = 2; i < indices.size(); i++) {
                Triangle tri;
                tri.v0 = vertices[indices[0]];
                tri.v1 = vertices[indices[i-1]];
                tri.v2 = vertices[indices[i]];
                triangles.push_back(tri);
            }
        }
    }

    if (triangles.empty()) return false;
    build_bvh();
    return true;
}

// ─── BVH construction ────────────────────────────────────────────────────────

void TriangleMesh::build_bvh(int max_leaf_tris) {
    if (triangles.empty()) return;

    bvh_nodes.clear();
    tri_indices.resize(triangles.size());
    for (size_t i = 0; i < triangles.size(); i++) tri_indices[i] = (int)i;

    bvh_nodes.reserve(2 * triangles.size());
    build_bvh_recursive(0, (int)triangles.size(), 0);
}

int TriangleMesh::build_bvh_recursive(int start, int count, int depth) {
    BVHNode node;
    node.bounds = AABB();
    for (int i = start; i < start + count; i++) {
        node.bounds.expand(triangles[tri_indices[i]].bounds());
    }

    if (count <= 4 || depth > 30) {
        // Leaf
        node.tri_start = start;
        node.tri_count = count;
        int idx = (int)bvh_nodes.size();
        bvh_nodes.push_back(node);
        return idx;
    }

    // Find split axis (longest)
    vec3 size = node.bounds.size();
    int axis = 0;
    if (size.y > size.x) axis = 1;
    if (size.z > size[axis]) axis = 2;

    // Sort by centroid along axis
    std::sort(tri_indices.begin() + start, tri_indices.begin() + start + count,
        [&](int a, int b) {
            return triangles[a].bounds().center()[axis] < triangles[b].bounds().center()[axis];
        });

    int mid = count / 2;
    int idx = (int)bvh_nodes.size();
    bvh_nodes.push_back(node);  // placeholder

    int left = build_bvh_recursive(start, mid, depth + 1);
    int right = build_bvh_recursive(start + mid, count - mid, depth + 1);

    bvh_nodes[idx].left = left;
    bvh_nodes[idx].right = right;
    return idx;
}

// ─── BVH sphere query ────────────────────────────────────────────────────────

TriangleMesh::SphereQueryResult TriangleMesh::query_sphere(const vec3& pos, double radius) const {
    SphereQueryResult best;
    best.penetration = -1e30;
    best.hit = false;

    if (bvh_nodes.empty()) return best;
    query_sphere_recursive(0, pos, radius, best);
    return best;
}

void TriangleMesh::query_sphere_recursive(int node_idx, const vec3& pos, double radius,
                                           SphereQueryResult& best) const {
    const BVHNode& node = bvh_nodes[node_idx];

    // Expand AABB by radius for sphere test
    AABB expanded = node.bounds;
    expanded.pad(radius);
    if (!expanded.contains(pos)) return;

    if (node.is_leaf()) {
        for (int i = node.tri_start; i < node.tri_start + node.tri_count; i++) {
            const Triangle& tri = triangles[tri_indices[i]];
            auto cr = closest_point_on_triangle(pos, tri);
            if (!cr.valid) continue;

            double pen = radius - cr.distance;
            if (pen > best.penetration) {
                best.penetration = pen;
                best.contact_point = cr.closest_point;
                // Normal: from triangle surface toward sphere center
                vec3 diff = pos - cr.closest_point;
                double d = diff.length();
                if (d > 1e-12) {
                    best.contact_normal = diff / d;
                } else {
                    best.contact_normal = cr.normal;
                }
                best.hit = (pen > 0);
            }
        }
    } else {
        if (node.left >= 0) query_sphere_recursive(node.left, pos, radius, best);
        if (node.right >= 0) query_sphere_recursive(node.right, pos, radius, best);
    }
}

AABB TriangleMesh::total_bounds() const {
    if (!bvh_nodes.empty()) return bvh_nodes[0].bounds;
    AABB b;
    for (auto& t : triangles) b.expand(t.bounds());
    return b;
}

// ─── Analytic primitives ─────────────────────────────────────────────────────

AnalyticPrimitive::CollisionResult AnalyticPrimitive::collide_sphere(const vec3& pos, double r) const {
    switch (type) {
        case PrimitiveType::CYLINDER:     return collide_cylinder(pos, r);
        case PrimitiveType::CONE:         return collide_cone(pos, r);
        case PrimitiveType::BENT_CYLINDER: return collide_bent_cylinder(pos, r);
        case PrimitiveType::SPIRAL_INSERT: return collide_spiral(pos, r);
        case PrimitiveType::PLANE:        return collide_plane(pos, r);
    }
    return {};
}

AnalyticPrimitive::CollisionResult AnalyticPrimitive::collide_cylinder(const vec3& pos, double r) const {
    // Infinite cylinder along axis, then clamp to length
    vec3 d = pos - center;
    double along = dot(d, axis);

    // Check if within cylinder length
    if (along < 0 || along > length) return {};

    // Radial distance from axis
    vec3 radial = d - axis * along;
    double radial_dist = radial.length();

    // Inside cylinder (tube): contact with inner wall when particle reaches wall
    // Wall is at distance `radius` from axis. Particle surface reaches wall when
    // radial_dist + r > radius, i.e., penetration = (radial_dist + r) - radius > 0
    double penetration = (radial_dist + r) - radius;
    if (penetration > 0 && radial_dist > 1e-12) {
        CollisionResult res;
        res.penetration = penetration;
        // Normal points from wall toward tube center (pushes particle inward)
        res.contact_normal = -radial.normalized();
        res.contact_point = center + axis * along + radial.normalized() * radius;
        res.hit = true;
        return res;
    }
    return {};
}

AnalyticPrimitive::CollisionResult AnalyticPrimitive::collide_cone(const vec3& pos, double r) const {
    vec3 d = pos - center;
    double along = dot(d, axis);
    if (along < 0 || along > length) return {};

    double t = along / length;
    double local_radius = radius_top * (1.0 - t) + radius_bottom * t;

    vec3 radial = d - axis * along;
    double radial_dist = radial.length();
    double penetration = (radial_dist + r) - local_radius;

    if (penetration > 0 && radial_dist > 1e-12) {
        CollisionResult res;
        res.penetration = penetration;
        res.contact_normal = -radial.normalized();
        res.contact_point = center + axis * along + radial.normalized() * local_radius;
        res.hit = true;
        return res;
    }
    return {};
}

AnalyticPrimitive::CollisionResult AnalyticPrimitive::collide_bent_cylinder(const vec3& pos, double r) const {
    // Bent cylinder: circular arc in the XZ plane
    // The arc center is at center + (bend_radius, 0, 0)
    vec3 arc_center = center + vec3(bend_radius, 0, 0);

    // Project point onto the arc plane
    vec3 d = pos - arc_center;
    // Angle along the arc (in XZ plane)
    double theta = std::atan2(-d.z, -d.x);
    if (theta < 0) theta += 2.0 * M_PI;

    if (theta > bend_angle) return {};

    // Point on arc centerline
    vec3 arc_point = arc_center + vec3(-bend_radius * std::cos(theta), 0, -bend_radius * std::sin(theta));

    // Radial distance from arc centerline
    vec3 to_particle = pos - arc_point;
    // Remove the component along the arc tangent
    vec3 tangent = vec3(bend_radius * std::sin(theta), 0, -bend_radius * std::cos(theta)).normalized();
    double along_tangent = dot(to_particle, tangent);
    vec3 radial = to_particle - tangent * along_tangent;

    double radial_dist = radial.length();
    double penetration = (radial_dist + r) - radius;

    if (penetration > 0 && radial_dist > 1e-12) {
        CollisionResult res;
        res.penetration = penetration;
        res.contact_normal = -radial.normalized();
        res.contact_point = arc_point + radial.normalized() * radius;
        res.hit = true;
        return res;
    }
    return {};
}

AnalyticPrimitive::CollisionResult AnalyticPrimitive::collide_spiral(const vec3& pos, double r) const {
    // Spiral insert: helical surface inside a cylinder
    // First check cylinder containment
    auto cyl_result = collide_cylinder(pos, r);

    // Now check spiral fins
    vec3 d = pos - center;
    double along = dot(d, axis);
    if (along < 0 || along > length) return cyl_result;

    // Spiral angle at this height
    double theta_offset = (along / length) * spiral_turns * 2.0 * M_PI;
    vec3 radial = d - axis * along;
    double radial_dist = radial.length();

    if (radial_dist < 1e-10) return cyl_result;

    // Angle of particle in radial plane
    double theta_p = std::atan2(radial.y, radial.x);

    // Check proximity to spiral fin
    double fin_theta = theta_offset;
    double angle_diff = theta_p - fin_theta;
    // Normalize to [-pi, pi]
    while (angle_diff > M_PI) angle_diff -= 2.0 * M_PI;
    while (angle_diff < -M_PI) angle_diff += 2.0 * M_PI;

    // Fin is a thin surface at angle fin_theta
    double angular_distance = std::abs(angle_diff) * radial_dist;
    double fin_gap = angular_distance - r;

    if (fin_gap < 0 && radial_dist > 0.3 * radius) {
        // Contact with spiral fin
        CollisionResult res;
        res.penetration = -fin_gap;
        // Normal perpendicular to fin surface (tangential direction)
        double sign = (angle_diff > 0) ? 1.0 : -1.0;
        vec3 tang_dir = vec3(-std::sin(fin_theta), std::cos(fin_theta), 0) * sign;
        res.contact_normal = tang_dir;
        res.contact_point = pos - tang_dir * angular_distance;
        res.hit = true;

        // Return whichever has more penetration
        if (!cyl_result.hit || res.penetration > cyl_result.penetration) return res;
    }

    return cyl_result;
}

AnalyticPrimitive::CollisionResult AnalyticPrimitive::collide_plane(const vec3& pos, double r) const {
    double dist = dot(pos - center, plane_normal) - plane_offset;
    double pen = r - dist;
    if (pen > 0) {
        CollisionResult res;
        res.penetration = pen;
        res.contact_normal = plane_normal;
        res.contact_point = pos - plane_normal * dist;
        res.hit = true;
        return res;
    }
    return {};
}

// ─── Geometry aggregate ──────────────────────────────────────────────────────

Geometry::ContactResult Geometry::query_sphere(const vec3& pos, double radius) const {
    ContactResult best;
    best.penetration = -1e30;
    best.hit = false;

    for (auto& mesh : meshes) {
        auto r = mesh.query_sphere(pos, radius);
        if (r.hit && r.penetration > best.penetration) {
            best.contact_point = r.contact_point;
            best.contact_normal = r.contact_normal;
            best.penetration = r.penetration;
            best.hit = true;
        }
    }
    for (auto& prim : primitives) {
        auto r = prim.collide_sphere(pos, radius);
        if (r.hit && r.penetration > best.penetration) {
            best.contact_point = r.contact_point;
            best.contact_normal = r.contact_normal;
            best.penetration = r.penetration;
            best.hit = true;
        }
    }
    return best;
}

bool Geometry::load_stl(const std::string& path) {
    TriangleMesh mesh;
    if (!mesh.load_stl(path)) return false;
    meshes.push_back(std::move(mesh));
    return true;
}

void Geometry::add_primitive(const AnalyticPrimitive& prim) {
    primitives.push_back(prim);
}

} // namespace dem

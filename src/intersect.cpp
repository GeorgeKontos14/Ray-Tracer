#include "intersect.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <limits>


bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p) {
    // TODO: implement this function.
    glm::vec3 nb = glm::cross(v0 - v2, p - v2);
    glm::vec3 nc = glm::cross(v1 - v0, p - v0);
    float beta = glm::dot(n, nb) / glm::dot(n, n);
    if (beta < 0)
        return false;
    float gamma = glm::dot(n, nc) / glm::dot(n, n);
    if (gamma < 0)
        return false;
    if (beta + gamma <= 1)
        return true;
    return false;
}

bool intersectRayWithPlane(const Plane& plane, Ray& ray)
{
    if (glm::dot(plane.normal, ray.direction) == 0) {
        return false;
    }

    float newT = (plane.D - glm::dot(plane.normal, ray.origin)) / glm::dot(plane.normal, ray.direction);
    if (newT <= 0 || newT >= ray.t)
        return false;
    ray.t = std::min(ray.t, newT);
    return true;
}

Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
    Plane plane;
    plane.normal = glm::cross(v1 - v0, v2 - v0);
    plane.D = glm::dot(plane.normal, v0);
    return plane;
}

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, Ray& ray, HitInfo& hitInfo)
{
    Plane plane = trianglePlane(p.v0.position, p.v1.position, p.v2.position);
    Ray copyRay;
    copyRay.direction = ray.direction;
    copyRay.origin = ray.origin;
    copyRay.t = ray.t;
    bool intersectPlane = intersectRayWithPlane(plane, copyRay);
    if (!intersectPlane)
        return false;
    glm::vec3 point = copyRay.origin + copyRay.t * copyRay.direction;
    bool isInsideTriangle = pointInTriangle(p.v0.position, p.v1.position, p.v2.position, plane.normal, point);
    if (isInsideTriangle) {
        ray.t = std::min(copyRay.t, ray.t);
        return true;
    }
    return false;
}

/// Input: a sphere with the following attributes: sphere.radius, sphere.center
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const Sphere& sphere, Ray& ray, HitInfo& hitInfo)
{
    // TODO: implement this function.
    return false;
}

/// Input: an axis-aligned bounding box with the following parameters: minimum coordinates box.lower and maximum coordinates box.upper
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const AxisAlignedBox& box, Ray& ray)
{
    float xmin = (box.lower.x - ray.origin.x) / ray.direction.x;
    float ymin = (box.lower.y - ray.origin.y) / ray.direction.y;
    float zmin = (box.lower.z - ray.origin.z) / ray.direction.z;
    float xmax = (box.upper.x - ray.origin.x) / ray.direction.x;
    float ymax = (box.upper.y - ray.origin.y) / ray.direction.y;
    float zmax = (box.upper.z - ray.origin.z) / ray.direction.z;

    float inx;
    float iny;
    float inz;
    float outx;
    float outy;
    float outz;

    inx = glm::min(xmin, xmax);
    iny = glm::min(ymin, ymax);
    inz = glm::min(zmin, zmax);
    outx = glm::max(xmax, xmin);
    outy = glm::max(ymax, ymin);
    outz = glm::max(zmax, zmin);

    float in = glm::max(inx, glm::max(iny, inz));
    float out = glm::min(outx, glm::min(outy, outz));

    if (in > out || out < 0)
        return false;

    if (in < ray.t) {
        // ray.t = in;
        return true;
    }
    return false;
    // Time Complexity: O(1)
    // Space Complexity: O(1)
}

#include "bvh.h"
#include "draw.h"
#include "interpolate.h"
#include "intersect.h"
#include "render.h"
#include "scene.h"
#include "extra.h"
#include "texture.h"
#include <algorithm>
#include <bit>
#include <chrono>
#include <framework/opengl_includes.h>
#include <iostream>
#include <stack>

// Helper method to fill in hitInfo object. This can be safely ignored (or extended).
// Note: many of the functions in this helper tie in to standard/extra features you will have
// to implement separately, see interpolate.h/.cpp for these parts of the project
void updateHitInfo(RenderState& state, const BVHInterface::Primitive& primitive, const Ray& ray, HitInfo& hitInfo)
{
    const auto& [v0, v1, v2] = std::tie(primitive.v0, primitive.v1, primitive.v2);
    const auto& mesh = state.scene.meshes[primitive.meshID];
    const auto n = glm::normalize(glm::cross(v1.position - v0.position, v2.position - v0.position));
    const auto p = ray.origin + ray.t * ray.direction;

    // First, fill in default data, unrelated to separate features
    hitInfo.material = mesh.material;
    hitInfo.normal = n;
    hitInfo.barycentricCoord = computeBarycentricCoord(v0.position, v1.position, v2.position, p);

    // Next, if `features.enableNormalMapping` is true, generate smoothly interpolated vertex normals
    if (state.features.enableNormalInterp) {
        hitInfo.normal = interpolateNormal(v0.normal, v1.normal, v2.normal, hitInfo.barycentricCoord);
    }

    // Next, if `features.enableTextureMapping` is true, generate smoothly interpolated vertex uvs
    if (state.features.enableTextureMapping) {
        hitInfo.texCoord = interpolateTexCoord(v0.texCoord, v1.texCoord, v2.texCoord, hitInfo.barycentricCoord);
    }

    // Finally, catch flipped normals
    if (glm::dot(ray.direction, n) > 0.0f) {
        hitInfo.normal = -hitInfo.normal;
    }
}

// BVH constructor; can be safely ignored. You should not have to touch this
// NOTE: this constructor is tested, so do not change the function signature.
BVH::BVH(const Scene& scene, const Features& features)
{
#ifndef NDEBUG
    // Store start of bvh build for timing
    using clock = std::chrono::high_resolution_clock;
    const auto start = clock::now();
#endif

    // Count the total nr. of triangles in the scene
    size_t numTriangles = 0;
    for (const auto& mesh : scene.meshes)
        numTriangles += mesh.triangles.size();

    // Given the input scene, gather all triangles over which to build the BVH as a list of Primitives
    std::vector<Primitive> primitives;
    primitives.reserve(numTriangles);
    for (uint32_t meshID = 0; meshID < scene.meshes.size(); meshID++) {
        const auto& mesh = scene.meshes[meshID];
        for (const auto& triangle : mesh.triangles) {
            primitives.push_back(Primitive {
                .meshID = meshID,
                .v0 = mesh.vertices[triangle.x],
                .v1 = mesh.vertices[triangle.y],
                .v2 = mesh.vertices[triangle.z] });
        }
    }

    // Tell underlying vectors how large they should approximately be
    m_primitives.reserve(numTriangles);
    m_nodes.reserve(numTriangles + 1);

    // Recursively build BVH structure; this is where your implementation comes in
    m_nodes.emplace_back(); // Create root node
    m_nodes.emplace_back(); // Create dummy node s.t. children are allocated on the same cache line
    buildRecursive(scene, features, primitives, RootIndex);

    // Fill in boilerplate data
    buildNumLevels();
    buildNumLeaves();

#ifndef NDEBUG
    // Output end of bvh build for timing
    const auto end = clock::now();
    std::cout << "BVH construction time: " << std::chrono::duration<double, std::milli>(end - start).count() << "ms" << std::endl;
#endif
}

// BVH helper method; allocates a new node and returns its index
// You should not have to touch this
uint32_t BVH::nextNodeIdx()
{
    const auto idx = static_cast<uint32_t>(m_nodes.size());
    m_nodes.emplace_back();
    return idx;
}

//Helper method that calculates the minimum coordinates for a primitive
// - primitive; a single triangle
// - return; the minimum x,y,z coordinates of the primitive's vertices
glm::vec3 primitiveMin(const BVHInterface::Primitive primitive)
{
    float minX = std::min(primitive.v0.position.x, std::min(primitive.v1.position.x, primitive.v2.position.x));
    float minY = std::min(primitive.v0.position.y, std::min(primitive.v1.position.y, primitive.v2.position.y));
    float minZ = std::min(primitive.v0.position.z, std::min(primitive.v1.position.z, primitive.v2.position.z));
    return glm::vec3 { minX, minY, minZ };
    // Time complexity: O(1)
    // Space complexity: O(1)
}

// Helper method that calculates the maximum coordinates for a primitive
//  - primitive; a single triangle
//  - return; the maximum x,y,z coordinates of the primitive's vertices
glm::vec3 primitiveMax(const BVHInterface::Primitive primitive)
{
    float maxX = std::max(primitive.v0.position.x, std::max(primitive.v1.position.x, primitive.v2.position.x));
    float maxY = std::max(primitive.v0.position.y, std::max(primitive.v1.position.y, primitive.v2.position.y));
    float maxZ = std::max(primitive.v0.position.z, std::max(primitive.v1.position.z, primitive.v2.position.z));
    return glm::vec3 { maxX, maxY, maxZ };
    // Time complexity: O(1)
    // Space complexity: O(1)
}

// TODO: Standard feature
// Given a BVH triangle, compute an axis-aligned bounding box around the primitive
// - primitive; a single triangle to be stored in the BVH
// - return;    an axis-aligned bounding box around the triangle
// This method is unit-tested, so do not change the function signature.
AxisAlignedBox computePrimitiveAABB(const BVHInterface::Primitive primitive)
{
    //In order for the triangle to be inside the bounding box, we use the minimum and maximum coordinates of the
    //of the triangle in order to specify the lower and the upper corner respectively.

    //Calculating the lower corner, using a helper method
    glm::vec3 lower = primitiveMin(primitive);

    //Calculating the upper corner, using a helper methid
    glm::vec3 upper = primitiveMax(primitive);

    return { .lower = lower, .upper = upper };
    // Time complexity: O(1)
    // Space complexity: O(1)
}

// TODO: Standard feature
// Given a range of BVH triangles, compute an axis-aligned bounding box around the range.
// - primitive; a contiguous range of triangles to be stored in the BVH
// - return;    a single axis-aligned bounding box around the entire set of triangles
// This method is unit-tested, so do not change the function signature.
AxisAlignedBox computeSpanAABB(std::span<const BVHInterface::Primitive> primitives)
{
    //In order to compute a bounding box that wraps around all triangles in the given span, we can set the lower corner to the minimum of the 
    //coordinates of all the triangle vertices and the upper corner to the respective maximum.

    //First, we initialize minimum and maxmum variables for all three coordinates
    float minX = std::numeric_limits<float>::max();
    float minY = std::numeric_limits<float>::max();
    float minZ = std::numeric_limits<float>::max();
    float maxX = std::numeric_limits<float>::min();
    float maxY = std::numeric_limits<float>::min();
    float maxZ = std::numeric_limits<float>::min();

    //Then, we iterate through all primitives in the span, calculate its minimum and maximum coefficients using the same
    //helper methods as for computePrimitiveAABB. We compare these to the current extrema and update them if necessary
    glm::vec3 currentMin;
    glm::vec3 currentMax;
    for (BVHInterface::Primitive p : primitives) {
        currentMin = primitiveMin(p);
        currentMax = primitiveMax(p);
        minX = std::min(currentMin.x, minX);
        minY = std::min(currentMin.y, minY);
        minZ = std::min(currentMin.z, minZ);
        maxX = std::max(currentMax.x, maxX);
        maxY = std::max(currentMax.y, maxY);
        maxZ = std::max(currentMax.z, maxZ);
    }
    glm::vec3 lower = { minX, minY, minZ };
    glm::vec3 upper = { maxX, maxY, maxZ };
    return { .lower = lower, .upper = upper };
    // Time Complexity: O(n), where n = length of primitives span
    // Space complexity: O(1)
}

// TODO: Standard feature
// Given a BVH triangle, compute the geometric centroid of the triangle
// - primitive; a single triangle to be stored in the BVH
// - return;    the geometric centroid of the triangle's vertices
// This method is unit-tested, so do not change the function signature.
glm::vec3 computePrimitiveCentroid(const BVHInterface::Primitive primitive)
{
    //The centroid is the barycenter of the triangle; thus, its coefficients
    //are the average of the coefficients of the triangle
    float avgX = (primitive.v0.position.x + primitive.v1.position.x + primitive.v2.position.x) / 3;
    float avgY = (primitive.v0.position.y + primitive.v1.position.y + primitive.v2.position.y) / 3;
    float avgZ = (primitive.v0.position.z + primitive.v1.position.z + primitive.v2.position.z) / 3;
    return glm::vec3 { avgX, avgY, avgZ };
    // Time Complexity: O(1)
    // Space complexity: O(1)
}

// TODO: Standard feature
// Given an axis-aligned bounding box, compute the longest axis; x = 0, y = 1, z = 2.
// - aabb;   the input axis-aligned bounding box
// - return; 0 for the x-axis, 1 for the y-axis, 2 for the z-axis
//           if several axes are equal in length, simply return the first of these
// This method is unit-tested, so do not change the function signature.
uint32_t computeAABBLongestAxis(const AxisAlignedBox& aabb)
{
    //We compute the difference between the maximum and the minimum of the x, y, z coefficients to fing the 
    //length of the x, y and z axis respectively. Afterwards, we compare those to find the longest axis
    float xAxis = aabb.upper.x - aabb.lower.x;
    float yAxis = aabb.upper.y - aabb.lower.y;
    float zAxis = aabb.upper.z - aabb.lower.z;
    if (xAxis >= yAxis) {
        if (xAxis >= zAxis)
            return 0;
        return 2;
    } else {
        if (yAxis >= zAxis)
            return 1;
        return 2;
    }
    return 0;
    // Time Complexity: O(1)
    // Space Complexity: O(1)
}

// Helper method that sorts the primitives along the specified axis
// - axis; 0, 1, 2 for x, y, z axis respectively
// - data; vector containing all the primitives and their centroids
// The sorting algorithm used is merge sort
void sortPrimitiveData(uint32_t axis, std::vector<PrimitiveData>& data) 
{
    if (data.size() < 2)
        return;

    //Split the vector in two vectors, left and right
    std::vector<PrimitiveData> left;
    std::vector<PrimitiveData> right;
    for (int i = 0; i < data.size(); i++) {
        if (i < data.size() / 2)
            left.push_back(data[i]);
        else
            right.push_back(data[i]);
    }

    //Sort these two vectors recursively
    sortPrimitiveData(axis, left);
    sortPrimitiveData(axis, right);

    //Merge left and right
    data = mergeVectors(axis, left, right);
    return;
    // Time Complexity: O(nlog(n)), where n = data.size()
    // Space Complexity: O(n)
}

// Merges two sorted vectors of primtives
// - axis; 0, 1, 2 for x, y, z respectively
// - v1; the first sorted vector
// - v2; the second sorted vector
// - return; the sorted merged vector
std::vector<PrimitiveData> mergeVectors(uint32_t axis, std::vector<PrimitiveData> v1, std::vector<PrimitiveData> v2) 
{
    std::vector<PrimitiveData> merged;
    int i = 0;
    int j = 0;
    int res;
    PrimitiveData p1;
    PrimitiveData p2;
    while (i < v1.size() && j < v2.size()) {
        p1 = v1[i];
        p2 = v2[j];
        res = comparePrimitives(axis, p1, p2);
        if (res >= 0) {
            merged.push_back(p2);
            j++;
        } else {
            merged.push_back(p1);
            i++;
        }
    }
    while (i < v1.size()) {
        merged.push_back(v1[i]);
        i++;
    }
    while (j < v2.size()) {
        merged.push_back(v2[j]);
        j++;
    }
    return merged;
    // Time Complexity: O(max(m, n)), where m = v1.size(), n = v2.size()
    // Space Complexity: O(m+n)
}

// Helper method for testing
bool equalVectors(std::vector<PrimitiveData> v1, std::vector<PrimitiveData> v2) {
    if (v1.size() != v2.size())
        return false;
    for (int i = 0; i < v1.size(); i++)
        if (v1[i].primitive != v2[i].primitive)
            return false;
    return true;
}

bool equalSpans(std::span<BVHInterface::Primitive> v1, std::span<BVHInterface::Primitive> v2) {
    if (v1.size() != v2.size())
        return false;
    for (int i = 0; i < v1.size(); i++)
        if (v1[i] != v2[i])
            return false;
    return true;
}

// Compares two primitves along the specified axis
// - axis; 0, 1, 2 for x, y, z axis respectively
// - p1; the data of the first primitive
// - p2; the data of the second primitive
// - return; 1 if p1 > p2; 0 if p1 == p2; -1 if p1 < p2
int comparePrimitives(uint32_t axis, PrimitiveData p1, PrimitiveData p2)
{
    float c1; //coefficient of p1 used in the comparison
    float c2; //coefficient of p2 used in the comparison
    if (axis == 0) {
        c1 = p1.centroid.x;
        c2 = p2.centroid.x;
    } else if (axis == 1) {
        c1 = p1.centroid.y;
        c2 = p2.centroid.y;
    } else {
        c1 = p1.centroid.z;
        c2 = p2.centroid.z;
    }
    if (c1 > c2)
        return 1;
    else if (c1 < c2)
        return -1;
    
    return 0;
    // Time complexity: O(1)
    // Space complexity: O(1)
}

// TODO: Standard feature
// Given a range of BVH triangles, sort these along a specified axis based on their geometric centroid.
// Then, find and return the split index in the range, such that the subrange containing the first element 
// of the list is at least as big as the other, and both differ at most by one element in size.
// Hint: you should probably reuse `computePrimitiveCentroid()`
// - aabb;       the axis-aligned bounding box around the given triangle range
// - axis;       0, 1, or 2, determining on which axis (x, y, or z) the split must happen
// - primitives; the modifiable range of triangles that requires sorting/splitting along an axis
// - return;     the split position of the modified range of triangles
// This method is unit-tested, so do not change the function signature.
size_t splitPrimitivesByMedian(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVHInterface::Primitive> primitives)
{
    using Primitive = BVHInterface::Primitive;
    //At first, we create a vector of PrimitiveData objects for all the primitives;
    //This way, we can sort the vector and then just extract the sorted primitives back to the primitives array
    std::vector<PrimitiveData> data;
    PrimitiveData current;
    for (Primitive p : primitives) {
        current.primitive = p;
        current.centroid = computePrimitiveCentroid(p);
        data.push_back(current);
    }

    //Then, we call a helper method to sort this array
    sortPrimitiveData(axis, data);

    //We now modify the primitives span, so that the primitives are sorted along the given axis
    for (int i = 0; i < data.size(); i++)
        primitives[i] = data[i].primitive;

    //All we need to do now is find which index to return;
    //if the span has an even size then the split index should be size/2; this way,
    //both the first and the second subspans will have size of size/2. If, however,
    //the span has odd size, then the split index should be size/2+1, so that the 
    //first subspan has one element more than the second subspan
    int n = primitives.size();
    if (n % 2 == 1)
        return n / 2 + 1;
    return n / 2;
    // Time Complexity: O(nlog(n)), where n = primitives.size()
    // Space Complexity: O(n)
}

bool interBVH(const AxisAlignedBox& box, Ray& ray) {
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
        return true;
    }
    return false;
}


// TODO: Standard feature
// Hierarchy traversal routine; called by the BVH's intersect(),
// you must implement this method and implement it carefully!
//
// If `features.enableAccelStructure` is not enabled, the method should just iterate the BVH's
// underlying primitives (or the scene's geometry). The default imlpementation already does this.
// You will have to implement the part which actually traverses the BVH for a faster intersect,
// given that `features.enableAccelStructure` is enabled.
//
// This method returns `true` if geometry was hit, and `false` otherwise. On first/closest hit, the
// distance `t` in the `ray` object is updated, and information is updated in the `hitInfo` object.
//
// - state;    the active scene, and a user-specified feature config object, encapsulated
// - bvh;      the actual bvh which should be traversed for faster intersection
// - ray;      the ray intersecting the scene's geometry
// - hitInfo;  the return object, with info regarding the hit geometry
// - return;   boolean, if geometry was hit or not
//
// This method is unit-tested, so do not change the function signature.
bool intersectRayWithBVH(RenderState& state, const BVHInterface& bvh, Ray& ray, HitInfo& hitInfo)
{
    // Relevant data in the constructed BVH
    std::span<const BVHInterface::Node> nodes = bvh.nodes();
    std::span<const BVHInterface::Primitive> primitives = bvh.primitives();

    // Return value
    bool is_hit = false;

    if (state.features.enableAccelStructure || state.features.extra.enableBvhSahBinning) {
        // TODO: implement here your (probably stack-based) BVH traversal.
        //
        // Some hints (refer to bvh_interface.h either way). BVH nodes are packed, so the
        // data is not easily extracted. Helper methods are available, however:
        // - For a given node, you can test if the node is a leaf with `node.isLeaf()`.
        // - If the node is not a leaf, you can obtain the left/right children with `node.leftChild()` etc.
        // - If the node is a leaf, you can obtain the offset to and nr. of primitives in the bvh's list
        //   of underlying primitives with `node.primitiveOffset()` and `node.primitiveCount()`
        //
        // In short, you will have to step down the bvh, node by node, and intersect your ray
        // with the node's AABB. If this intersection passes, you should:
        // - if the node is a leaf, intersect with the leaf's primitives
        // - if the node is not a leaf, test the left and right children as well!
        //
        // Note that it is entirely possible for a ray to hit a leaf node, but not its primitives,
        // and it is likewise possible for a ray to hit both children of a node.
          
        BVHInterface::Primitive prim;
        std::stack<int> inds;
        inds.push(0);
        BVHInterface::Node n;
        int primInd = 0;
        while (!inds.empty()) {
            n = nodes[inds.top()];
            inds.pop();
            if (interBVH(n.aabb, ray)) {
                if (n.isLeaf()) {
                    //is_hit |= intersectLeaf(state, n, ray, primInd, hitInfo, primitives);
                    for (int i = n.primitiveOffset(); i < n.primitiveOffset() + n.primitiveCount(); i++) {
                        prim = primitives[i];
                        const auto& [v0, v1, v2] = std::tie(prim.v0, prim.v1, prim.v2);
                        if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                            updateHitInfo(state, primitives[i], ray, hitInfo);
                            is_hit = true;
                        }
                    }
                } else {
                    inds.push(n.leftChild());
                    inds.push(n.rightChild());
                }
            }
        }
        return is_hit;
    } else {
        // Naive implementation; simply iterates over all primitives
        for (const auto& prim : primitives) {
            const auto& [v0, v1, v2] = std::tie(prim.v0, prim.v1, prim.v2);
            if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                updateHitInfo(state, prim, ray, hitInfo);
                is_hit = true;
            }
        }
    }

    // Intersect with spheres.
    for (const auto& sphere : state.scene.spheres)
        is_hit |= intersectRayWithShape(sphere, ray, hitInfo);

    return is_hit;
}

// TODO: Standard feature
// Leaf construction routine; you should reuse this in in `buildRecursive()`
// Given an axis-aligned bounding box, and a range of triangles, generate a valid leaf object
// and store the triangles in the `m_primitives` vector.
// You are free to modify this function's signature, as long as the constructor builds a BVH
// - scene;      the active scene
// - features;   the user-specified features object
// - aabb;       the axis-aligned bounding box around the primitives beneath this leaf
// - primitives; the range of triangles to be stored for this leaf
BVH::Node BVH::buildLeafData(const Scene& scene, const Features& features, const AxisAlignedBox& aabb, std::span<Primitive> primitives)
{
    Node node;
    // TODO fill in the leaf's data; refer to `bvh_interface.h` for details
    
    //The bounding box of the node is the bounding box around all of its primitives
    node.aabb = aabb;
    //By having 1 as the first bit, we clarify that the node is a leaf, as intended in the file 'bvh_interface.h'. The rest of the value is the id of the 
    //first primitive of the span
    node.data[0] = node.LeafBit + m_primitives.size();
    //Since the node is a leaf, the second element of data should indicate the amount of primitives in the node, ie, the size of the primitives span
    node.data[1] = primitives.size();
    // Copy the current set of primitives to the back of the primitives vector
    std::copy(primitives.begin(), primitives.end(), std::back_inserter(m_primitives));

    return node;
    // Time Complexity: O(1)
    // Space Complexity: O(1)
}

// TODO: Standard feature
// Node construction routine; you should reuse this in in `buildRecursive()`
// Given an axis-aligned bounding box, and left/right child indices, generate a valid node object.
// You are free to modify this function's signature, as long as the constructor builds a BVH
// - scene;           the active scene
// - features;        the user-specified features object
// - aabb;            the axis-aligned bounding box around the primitives beneath this node
// - leftChildIndex;  the index of the node's left child in `m_nodes`
// - rightChildIndex; the index of the node's right child in `m_nodes`
BVH::Node BVH::buildNodeData(const Scene& scene, const Features& features, const AxisAlignedBox& aabb, uint32_t leftChildIndex, uint32_t rightChildIndex)
{
    Node node;
    // TODO fill in the node's data; refer to `bvh_interface.h` for details
    
    //The bounding box of the node is the bounding box around all of its primitives
    node.aabb = aabb;
    //Since the node is an inner node, the first element of data should be 0 followed by the index of the left child,
    //and the second element should be the index of the rightChild
    node.data[0] = leftChildIndex;
    node.data[1] = rightChildIndex;

    return node;
    // Time Complexity: O(1)
    // Space Complexity: O(1)
}

// TODO: Standard feature
// Hierarchy construction routine; called by the BVH's constructor,
// you must implement this method and implement it carefully!
//
// You should implement the other BVH standard features first, and this feature last, as you can reuse
// most of the other methods to assemble this part. There are detailed instructions inside the
// method which we recommend you follow.
//
// Arguments:
// - scene;      the active scene
// - features;   the user-specified features object
// - primitives; a range of triangles to be stored in the BVH
// - nodeIndex;  index of the node you are currently working on, this is already allocated
//
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildRecursive(const Scene& scene, const Features& features, std::span<Primitive> primitives, uint32_t nodeIndex)
{
    // WARNING: always use nodeIndex to index into the m_nodes array. never hold a reference/pointer,
    // because a push/emplace (in ANY recursive calls) might grow vectors, invalidating the pointers.

    // Compute the AABB of the current node.
    AxisAlignedBox aabb = computeSpanAABB(primitives);

    // As a starting point, we provide an implementation which creates a single leaf, and stores
    // all triangles inside it. You should remove or comment this, and work on your own recursive
    // construction algorithm that implements the following steps. Make sure to reuse the methods
    // you have previously implemented to simplify this process.
    //
    // 1. Determine if the node should be a leaf, when the nr. of triangles is less or equal to 4
    //    (hint; use the `LeafSize` constant)
    // 2. If it is a leaf, fill in the leaf's data, and store its range of triangles in `m_primitives`
    // 3. If it is a node:
    //    3a. Split the range of triangles along the longest axis into left and right subspans,
    //        using either median or SAH-Binning based on the `Features` object
    //    3b. Allocate left/right child nodes
    //        (hint: use `nextNodeIdx()`)
    //    3c. Fill in the current node's data; aabb, left/right child indices
    //    3d. Recursively build left/right child nodes over their respective triangles
    //        (hint; use `std::span::subspan()` to split into left/right ranges)

    // Just configure the current node as a giant leaf for now
    //m_nodes[nodeIndex] = buildLeafData(scene, features, aabb, primitives);

    //Count the number of primitives to determine whether the node is a leaf or not
    bool isLeafNode = primitives.size() <= LeafSize;
    //If it is a leaf, we call the buildLeafData method and the leaf's primitives to m_primitives
    if (isLeafNode) {
        BVH::Node leaf = buildLeafData(scene, features, aabb, primitives);
        m_nodes[nodeIndex] = leaf;
        return;
    }

    //If it is an inner node, we split the span of primitives
    //If the extra feature is disabled, we split by calculating the median
    uint32_t longest = computeAABBLongestAxis(aabb);
    int splitInd = 0;
    if (features.extra.enableBvhSahBinning) 
        splitInd = splitPrimitivesBySAHBin(aabb, longest, primitives);
    else 
        splitInd = splitPrimitivesByMedian(aabb, longest, primitives);
    if (splitInd < 0 || splitInd > primitives.size()) {
        BVH::Node leaf = buildLeafData(scene, features, aabb, primitives);
        m_nodes[nodeIndex] = leaf;
        return;
    }
    uint32_t leftInd = nextNodeIdx();
    uint32_t rightInd = nextNodeIdx();
    // We store the data of the current node in the m_nodes vector
    m_nodes[nodeIndex] = buildNodeData(scene, features, aabb, leftInd, rightInd);
    // Recursively build child nodes until we reach leaves
    buildRecursive(scene, features, primitives.subspan(0, splitInd), leftInd);
    buildRecursive(scene, features, primitives.subspan(splitInd, primitives.size() - splitInd), rightInd);
    // Time Complexity: O(n(log(n))^2), where n = primitives.size()
    // Space Complexity: O(nlog(n))
}

// TODO: Standard feature, or part of it
// Compute the nr. of levels in your hierarchy after construction; useful for `debugDrawLevel()`
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildNumLevels()
{
    // #levels = 1 if n <= 2; ceil(log2(n))-1 else
    // where n = #primitives
    int n = m_primitives.size();
    if (n <= 2) {
        m_numLevels = 1;
        return;
    }
    int ceil = std::_Ceiling_of_log_2(n);
    m_numLevels = ceil - 1;
    return;
    // Time Complexity: O(1)
    // Space Complexity: O(1)
}

// Compute the nr. of leaves in your hierarchy after construction; useful for `debugDrawLeaf()`
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildNumLeaves()
{
    // #leaves = 2^(floor(log2(n))-2) + n - 2^(floor(log2(n))) if 2^floor(log2(n)) <= n <= 2^floor(log2(n))+2^(floor(log2(n))-2);
    // #leaves = 2^(floor(log2(n))-1) else
    // where n = # primitives
    int n = m_primitives.size();
    if (n <= 4) {
        m_numLeaves = 1;
        return;
    } else if (n <= 8) {
        m_numLeaves = 2;
        return;
    }
    int floor = std::_Floor_of_log_2(n);
    if (n >= pow(2, floor) && n <= pow(2, floor) + pow(2, floor - 2)) {
        m_numLeaves = pow(2, floor - 2) + n - pow(2, floor);
        return;
    }
    m_numLeaves = pow(2, floor - 1);
    return;
    // Time Complexity: O(1)
    // Space Complexity: O(1)
}

// Helper method that finds the next level of nodes in our BVH
// - currentLevelIndices; the indices of the nodes of the current level
// - nodes; the vector containing all the nodes of the BVH
// - return; a vector containing the indices of the nodes of the lower level in our BVH
std::vector<int> nextLevelIndices(std::vector<int> currentLevelIndices, std::vector<BVHInterface::Node> nodes) {
    std::vector<int> res;
    BVHInterface::Node current;
    for (int ind : currentLevelIndices) {
        current = nodes[ind];
        if (!current.isLeaf()) {
            res.push_back(current.leftChild());
            res.push_back(current.rightChild());
        }
    }
    return res;
    // Time Complexity: O(n), where n = #nodes on the current level
    // Space Complexity: O(2n)
}

// Draw the bounding boxes of the nodes at the selected level. Use this function to visualize nodes
// for debugging. You may wish to implement `buildNumLevels()` first. We suggest drawing the AABB
// of all nodes on the selected level.
// You are free to modify this function's signature.
void BVH::debugDrawLevel(int level)
{
    // Example showing how to draw an AABB as a (white) wireframe box.
    // Hint: use draw functions (see `draw.h`) to draw the contained boxes with different
    // colors, transparencies, etc.
    //AxisAlignedBox aabb { .lower = glm::vec3(0.0f), .upper = glm::vec3(0.0f, 1.05f, 1.05f) };
    //drawAABB(aabb, DrawMode::Wireframe, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
    if (level == 0) {
        BVHInterface::Node root = m_nodes[RootIndex];
        drawAABB(root.aabb, DrawMode::Wireframe, glm::vec3(0, 1.0f, 0), 1.0f);
        return;
    }
    std::vector<int> currentLevel;
    currentLevel.push_back(RootIndex);
    for (int i = 0; i < level; i++)
        currentLevel = nextLevelIndices(currentLevel, m_nodes);
    for (int ind : currentLevel)
        drawAABB(m_nodes[ind].aabb, DrawMode::Wireframe, glm::vec3(0, 1.0f, 0), 1.0f);
    return;
    // Time Complexity: O(2^level)
    // Space Complexity: O(2^level)
}

// Draw data of the leaf at the selected index. Use this function to visualize leaf nodes
// for debugging. You may wish to implement `buildNumLeaves()` first. We suggest drawing the AABB
// of the selected leaf, and then its underlying primitives with different colors.
// - leafIndex; index of the selected leaf.
//              (Hint: not the index of the i-th node, but of the i-th leaf!)
// You are free to modify this function's signature.
void BVH::debugDrawLeaf(int leafIndex)
{
    // Example showing how to draw an AABB as a (white) wireframe box.
    // Hint: use drawTriangle (see `draw.h`) to draw the contained primitives
    //AxisAlignedBox aabb { .lower = glm::vec3(0.0f), .upper = glm::vec3(0.0f, 1.05f, 1.05f) };
    //drawAABB(aabb, DrawMode::Wireframe, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
    if (leafIndex == 0)
        return;
    
    BVHInterface::Node current;
    int leafCount = 0;
    for (BVHInterface::Node n: m_nodes) {
        if (n.isLeaf())
            leafCount++;
        if (leafCount == leafIndex) {
            current = n;
            break;
        }
    }
    drawAABB(current.aabb, DrawMode::Wireframe, glm::vec3(1.0f, 0, 0), 1.0f);
    int count = current.primitiveCount();
    BVHInterface::Primitive p;
    for (int i = current.primitiveOffset(); i < current.primitiveOffset() + count; i++) {
        p = m_primitives[i];
        drawTriangle(p.v0, p.v1, p.v2, glm::vec3{0, 0, 1.0f});
    }
    // Time Complexity: O(n), where n = #nodes
    // Space Complexity: O(1)
}

void BVH::debugDrawSplit(float splitLine, int nodeIdx, uint32_t axis) {
    BVHInterface::Node n = m_nodes[nodeIdx];
    drawAABB(n.aabb, DrawMode::Wireframe, glm::vec3(0, 0, 1.0f), 1.0f);
    float length;
    float splitPos;
    AxisAlignedBox splitPlane {
        n.aabb.lower,
        n.aabb.upper
    };
    if (axis == 0) {
        length = n.aabb.upper.x - n.aabb.lower.x;
        splitPos = n.aabb.lower.x + splitLine * length;
        splitPlane.lower.x = splitPos;
        splitPlane.upper.x = splitPos;
        drawAABB(splitPlane, DrawMode::Filled, glm::vec3(0.8f, 0, 0.8f), 0.4f);
    } else if (axis == 1) {
        length = n.aabb.upper.y - n.aabb.lower.y;
        splitPos = n.aabb.lower.y + splitLine * length;
        splitPlane.lower.y = splitPos;
        splitPlane.upper.y = splitPos;
        drawAABB(splitPlane, DrawMode::Filled, glm::vec3(0.8f, 0, 0.8f), 0.4f);
    } else {
        length = n.aabb.upper.z - n.aabb.lower.z;
        splitPos = n.aabb.lower.z + splitLine * length;
        splitPlane.lower.z = splitPos;
        splitPlane.upper.z = splitPos;
        drawAABB(splitPlane, DrawMode::Filled, glm::vec3(0.8f, 0, 0.8f), 0.4f);
    }
    
}

void BVH::debugDrawOptimalSplit() {
    using Primitive = BVH::Primitive;

    std::span<Primitive> prims(m_primitives);
    AxisAlignedBox box = computeSpanAABB(prims);
    drawAABB(box, DrawMode::Wireframe, glm::vec3(0, 0, 1.0f), 1.0f);
    uint32_t axis = computeAABBLongestAxis(box);
    int numBins = 20;
    std::vector<float> centroidCoord;
    float step = 0;
    float lowerCoord = 0;
    if (axis == 0) {
        lowerCoord = box.lower.x;
        step = (box.upper.x - box.lower.x) / numBins;
        for (Primitive p : m_primitives)
            centroidCoord.push_back(computePrimitiveCentroid(p).x);
    } else if (axis == 1) {
        lowerCoord = box.lower.y;
        step = (box.upper.y - box.lower.y) / numBins;
        for (Primitive p : m_primitives)
            centroidCoord.push_back(computePrimitiveCentroid(p).y);
    } else {
        lowerCoord = box.lower.z;
        step = (box.upper.z - box.lower.z) / numBins;
        for (Primitive p : m_primitives)
            centroidCoord.push_back(computePrimitiveCentroid(p).z);
    }
    int nA;
    int nB;
    float currentCost;
    float minCost = std::numeric_limits<float>::max();
    float minSplitLine;
    int minNA = 0;
    for (int i = 1; i < numBins; i++) {
        nA = 0;
        nB = 0;
        for (float c : centroidCoord) {
            if (c <= lowerCoord + i * step)
                nA++;
            else
                nB++;
        }
        currentCost = costOfSplit(box, i * step, axis, nA, nB);
        if (currentCost < minCost) {
            minCost = currentCost;
            minSplitLine = i * step;
            minNA = nA;
        }
    }
    AxisAlignedBox aabb {
        box.lower,
        box.upper
    };
    if (minNA >= m_primitives.size()) {
        if (axis == 0)
            aabb.upper.x = aabb.lower.x;
        else if (axis == 1)
            aabb.upper.y = aabb.lower.y;
        else
            aabb.upper.z = aabb.lower.z;
    } else if (minNA == 0) {
        if (axis == 0)
            aabb.lower.x = aabb.upper.x;
        else if (axis == 1)
            aabb.lower.y = aabb.upper.y;
        else
            aabb.lower.z = aabb.upper.z;
    } else {
        if (axis == 0) {
            aabb.lower.x = lowerCoord + minSplitLine;
            aabb.upper.x = lowerCoord + minSplitLine;
        } else if (axis == 0) {
            aabb.lower.y = lowerCoord + minSplitLine;
            aabb.upper.y = lowerCoord + minSplitLine;
        } else {
            aabb.lower.z = lowerCoord + minSplitLine;
            aabb.upper.z = lowerCoord + minSplitLine;
        }
    }
    drawAABB(aabb, DrawMode::Filled, glm::vec3 { 1.0f, 1.0f, 0 }, 0.4f);
}
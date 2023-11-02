#include "extra.h"
#include "bvh.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include <framework/trackball.h>

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of Depth of Field. Here, you generate camera rays s.t. a focus point and a thin lens camera model
// are in play, allowing objects to be in and out of focus.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithDepthOfField(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableDepthOfField) {
        return;
    }

    // ...
}

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of motion blur. Here, you integrate over a time domain, and not just the pixel's image domain,
// to give objects the appearance of "fast movement".
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithMotionBlur(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableMotionBlur) {
        return;
    }

}

// TODO; Extra feature
// Given a rendered image, compute and apply a bloom post-processing effect to increase bright areas.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& image)
{
    if (!features.extra.enableBloomEffect)
        return;
    int n = 7;
    std::vector<float> filter;
    for (int i = 0; i < n; i++)
        filter.push_back(combinations(n, i) / (pow(2, n) - 1));
    Screen thres = Screen(image.resolution(), true);
    glm::vec3 pixel;
    for (int j = 0; j < image.pixels().size(); j++) {
        pixel = image.pixels()[j];
        if (pixel.x >= 0.9f || pixel.y >= 0.9f || pixel.z >= 0.9f)
            thres.pixels()[j] = pixel;
    }

    Screen filtered = Screen(image.resolution(), true);
    for (int j = 0; j < thres.resolution().y; j++) {
        for (int i = 0; i < thres.resolution().x; i++) {
            pixel = thres.pixels()[thres.indexAt(i, j)];
            if (pixel.x >= 0.9 || pixel.y >= 0.9 || pixel.z >= 0.9)
                applyFilter(thres, filter, i, j, true, filtered);
        }
    }


    for (int j = 0; j < thres.resolution().y; j++) {
        for (int i = 0; i < thres.resolution().x; i++) {
            pixel = thres.pixels()[thres.indexAt(i, j)];
            if (pixel.x >= 0.9 || pixel.y >= 0.9 || pixel.z >= 0.9)
                applyFilter(thres, filter, i, j, false, filtered);
        }
    }

    for (int i = 0; i < image.pixels().size(); i++) {
        image.pixels()[i] = image.pixels()[i] + filtered.pixels()[i];  
    }
}

void applyFilter(Screen& image, std::vector<float> filter, int x, int y, bool horizontal, Screen& after) {
    glm::vec3 res;
    if (horizontal) {
        if (x == 0)
            res = image.pixels()[image.indexAt(x, y)] * filter[3] + image.pixels()[image.indexAt(x + 1, y)] * filter[4] + image.pixels()[image.indexAt(x + 2, y)] * filter[5] + image.pixels()[image.indexAt(x + 3, y)] * filter[6]; 
        else if (x == 1)
            res = image.pixels()[image.indexAt(x - 1, y)] * filter[2] + image.pixels()[image.indexAt(x, y)] * filter[3] + image.pixels()[image.indexAt(x + 1, y)] * filter[4] + image.pixels()[image.indexAt(x + 2, y)] * filter[5] + image.pixels()[image.indexAt(x + 3, y)] * filter[6];
        else if (x == 2)
            res = image.pixels()[image.indexAt(x - 2, y)] * filter[1] + image.pixels()[image.indexAt(x - 1, y)] * filter[2] + image.pixels()[image.indexAt(x, y)] * filter[3] + image.pixels()[image.indexAt(x + 1, y)] * filter[4] + image.pixels()[image.indexAt(x + 2, y)] * filter[5] + image.pixels()[image.indexAt(x + 3, y)] * filter[6];
        else if (x == image.resolution().x-1)
            res = image.pixels()[image.indexAt(x - 3, y)] * filter[0] + image.pixels()[image.indexAt(x - 2, y)] * filter[1] + image.pixels()[image.indexAt(x - 1, y)] * filter[2] + image.pixels()[image.indexAt(x, y)] * filter[3];
        else if (x == image.resolution().x-2)
            res = image.pixels()[image.indexAt(x - 3, y)] * filter[0] + image.pixels()[image.indexAt(x - 2, y)] * filter[1] + image.pixels()[image.indexAt(x - 1, y)] * filter[2] + image.pixels()[image.indexAt(x, y)] * filter[3] + image.pixels()[image.indexAt(x + 1, y)] * filter[4];
        else if (x == image.resolution().x-3)
            res = image.pixels()[image.indexAt(x - 3, y)] * filter[0] + image.pixels()[image.indexAt(x - 2, y)] * filter[1] + image.pixels()[image.indexAt(x - 1, y)] * filter[2] + image.pixels()[image.indexAt(x, y)] * filter[3] + image.pixels()[image.indexAt(x + 1, y)] * filter[4] + image.pixels()[image.indexAt(x + 2, y)] * filter[5];
        else
            res = image.pixels()[image.indexAt(x - 3, y)] * filter[0] +image.pixels()[image.indexAt(x - 2, y)] * filter[1] + image.pixels()[image.indexAt(x - 1, y)] * filter[2] + image.pixels()[image.indexAt(x, y)] * filter[3] + image.pixels()[image.indexAt(x + 1, y)] * filter[4] + image.pixels()[image.indexAt(x + 2, y)] * filter[5] + image.pixels()[image.indexAt(x + 3, y)] * filter[6];
    } else {
        if (y == 0)
            res = image.pixels()[image.indexAt(x, y)] * filter[3] + image.pixels()[image.indexAt(x, y+1)] * filter[4] + image.pixels()[image.indexAt(x, y+2)] * filter[5] + image.pixels()[image.indexAt(x, y+3)] * filter[6];
        else if (y == 1)
            res = image.pixels()[image.indexAt(x, y-1)] * filter[2] + image.pixels()[image.indexAt(x, y)] * filter[3] + image.pixels()[image.indexAt(x, y+1)] * filter[4] + image.pixels()[image.indexAt(x, y+2)] * filter[5] + image.pixels()[image.indexAt(x, y+3)] * filter[6];
        else if (y == 2)
            res = image.pixels()[image.indexAt(x, y-2)] * filter[1] + image.pixels()[image.indexAt(x, y-1)] * filter[2] + image.pixels()[image.indexAt(x, y)] * filter[3] + image.pixels()[image.indexAt(x, y+1)] * filter[4] + image.pixels()[image.indexAt(x, y+2)] * filter[5] + image.pixels()[image.indexAt(x, y+3)] * filter[6];
        else if (y == image.resolution().y - 1)
            res = image.pixels()[image.indexAt(x, y-3)] * filter[0] + image.pixels()[image.indexAt(x, y-2)] * filter[1] + image.pixels()[image.indexAt(x, y-1)] * filter[2] + image.pixels()[image.indexAt(x, y)] * filter[3];
        else if (y == image.resolution().y - 2)
            res = image.pixels()[image.indexAt(x, y-3)] * filter[0] + image.pixels()[image.indexAt(x, y-2)] * filter[1] + image.pixels()[image.indexAt(x, y-1)] * filter[2] + image.pixels()[image.indexAt(x, y)] * filter[3] + image.pixels()[image.indexAt(x, y+1)] * filter[4];
        else if (y == image.resolution().y - 3)
            res = image.pixels()[image.indexAt(x, y-3)] * filter[0] + image.pixels()[image.indexAt(x, y-2)] * filter[1] + image.pixels()[image.indexAt(x, y-1)] * filter[2] + image.pixels()[image.indexAt(x, y)] * filter[3] + image.pixels()[image.indexAt(x, y+1)] * filter[4] + image.pixels()[image.indexAt(x, y+2)];
        else
            res = image.pixels()[image.indexAt(x, y-3)] * filter[0] + image.pixels()[image.indexAt(x, y-2)] * filter[1] + image.pixels()[image.indexAt(x, y-1)] * filter[2] + image.pixels()[image.indexAt(x, y)] * filter[3] + image.pixels()[image.indexAt(x, y+1)] * filter[4] + image.pixels()[image.indexAt(x, y+2)] * filter[5] + image.pixels()[image.indexAt(x, y+3)] * filter[6];
    }
    after.pixels()[image.indexAt(x, y)] = res;

}

unsigned long long int combinations(uint32_t n, uint32_t k)
{
    if (k == 0 || k == n)
        return 1;
    unsigned long long int res = n;
    for (uint32_t i = 1; i < k; i++)
        res *= n - i;
    res = res / factorial(k);
    return res;
}

unsigned long long int factorial(uint32_t k) {
    if (k <= 1)
        return 1;
    return k * factorial(k - 1);
}




// TODO; Extra feature
// Given a camera ray (or reflected camera ray) and an intersection, evaluates the contribution of a set of
// glossy reflective rays, recursively evaluating renderRay(..., depth + 1) along each ray, and adding the
// results times material.ks to the current intersection's hit color.
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderRayGlossyComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    // Generate an initial specular ray, and base secondary glossies on this ray
    // auto numSamples = state.features.extra.numGlossySamples;
    // ...
}

// TODO; Extra feature
// Given a camera ray (or reflected camera ray) that does not intersect the scene, evaluates the contribution
// along the ray, originating from an environment map. You will have to add support for environment textures
// to the Scene object, and provide a scene with the right data to supply this.
// - state; the active scene, feature config, bvh, and sampler
// - ray;   ray object
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
glm::vec3 sampleEnvironmentMap(RenderState& state, Ray ray)
{
    if (state.features.extra.enableEnvironmentMap) {
        // Part of your implementation should go here
        return glm::vec3(0.f);
    } else {
        return glm::vec3(0.f);
    }
}


// TODO: Extra feature
// As an alternative to `splitPrimitivesByMedian`, use a SAH+binning splitting criterion. Refer to
// the `Data Structures` lecture for details on this metric.
// - aabb;       the axis-aligned bounding box around the given triangle set
// - axis;       0, 1, or 2, determining on which axis (x, y, or z) the split must happen
// - primitives; the modifiable range of triangles that requires splitting
// - return;     the split position of the modified range of triangles
// This method is unit-tested, so do not change the function signature.
size_t splitPrimitivesBySAHBin(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVH::Primitive> primitives)
{
    using Primitive = BVH::Primitive;
    int numBins = 20;
    std::vector<float> centroidCoord;
    float step = 0;
    float lowerCoord = 0;
    if (axis == 0) {
        lowerCoord = aabb.lower.x;
        step = (aabb.upper.x - aabb.lower.x) / numBins;
        for (Primitive p : primitives)
            centroidCoord.push_back(computePrimitiveCentroid(p).x);
    } else if (axis == 1) {
        lowerCoord = aabb.lower.y;
        step = (aabb.upper.y - aabb.lower.y) / numBins;
        for (Primitive p : primitives)
            centroidCoord.push_back(computePrimitiveCentroid(p).y);
    } else {
        lowerCoord = aabb.lower.z;
        step = (aabb.upper.z - aabb.lower.z) / numBins;
        for (Primitive p : primitives)
            centroidCoord.push_back(computePrimitiveCentroid(p).z);
    }
    int nA;
    int nB;
    float currentCost;
    float minCost = std::numeric_limits<float>::max();
    float minSplitLine;
    int minNA;
    for (int i = 1; i < numBins; i++) {
        nA = 0;
        nB = 0;
        for (float c : centroidCoord) {
            if (c <= lowerCoord + i * step)
                nA++;
            else
                nB++;
        }
        currentCost = costOfSplit(aabb, i * step, axis, nA, nB);
        if (currentCost < minCost) {
            minCost = currentCost;
            minSplitLine = i * step;
            minNA = nA;
        }
    }
    if (minNA >= primitives.size() || minNA == 0)
        return -1;
    Primitive temp;
    float tempC;
    int i = 0;
    int j = primitives.size() - 1;
    while(i < primitives.size()) {
        if (i == minNA)
            break;
        if (centroidCoord[i] > lowerCoord + minSplitLine) {
            temp = primitives[i];
            primitives[i] = primitives[j];
            primitives[j] = temp;
            tempC = centroidCoord[i];
            centroidCoord[i] = centroidCoord[j];
            centroidCoord[j] = tempC;
            j--;
        } else {
            i++;
        }
    }
    return minNA;
}


// Helper method that calculates the cost of a given split
// - aabb; the AABB around all the primitives
// - split; offset from the lower value on the given axis for split line
// - axis; 0, 1, 2 for x, y, z respectively
// - nA; number of primitives in first partition
// - nB; number of primitives in second partition
// - return; the cost of the given split
// The cost of traversing the BVH is constant for all splits and thus is excluded from the calculations
float costOfSplit(const AxisAlignedBox& aabb, float split, uint32_t axis, int nA, int nB) {
    float S = 2*((aabb.upper.x-aabb.lower.x)*(aabb.upper.y-aabb.lower.y) + (aabb.upper.y - aabb.lower.y)*(aabb.upper.z-aabb.lower.z) +(aabb.upper.z-aabb.lower.z)* (aabb.upper.x-aabb.lower.x));
    float SA;
    float SB;
    if (axis == 0) { 
        SA = 2 * (split * (aabb.upper.y - aabb.lower.y) + (aabb.upper.y - aabb.lower.y) * (aabb.upper.z - aabb.lower.z) + (aabb.upper.z - aabb.lower.z) * split);
        SB = 2 * ((aabb.upper.x - aabb.lower.x - split) * (aabb.upper.y - aabb.lower.y) + (aabb.upper.y - aabb.lower.y) * (aabb.upper.z - aabb.lower.z) + (aabb.upper.z - aabb.lower.z) * (aabb.upper.x - aabb.lower.x - split));
    } else if (axis == 1) {
        SA = 2 * ((aabb.upper.x - aabb.lower.x) * split + split * (aabb.upper.z - aabb.lower.z) + (aabb.upper.z - aabb.lower.z) * (aabb.upper.x - aabb.lower.x));
        SB = 2 * ((aabb.upper.x - aabb.lower.x) * (aabb.upper.y - aabb.lower.y - split) + (aabb.upper.y - aabb.lower.y - split) * (aabb.upper.z - aabb.lower.z) + (aabb.upper.z - aabb.lower.z) * (aabb.upper.x - aabb.lower.x));
    } else {
        SA = 2 * ((aabb.upper.x - aabb.lower.x) * (aabb.upper.y - aabb.lower.y) + (aabb.upper.y - aabb.lower.y) * split + split * (aabb.upper.x - aabb.lower.x));
        SB = 2 * ((aabb.upper.x - aabb.lower.x) * (aabb.upper.y - aabb.lower.y) + (aabb.upper.y - aabb.lower.y) * (aabb.upper.z - aabb.lower.z - split) + (aabb.upper.z - aabb.lower.z - split) * (aabb.upper.x - aabb.lower.x));
    }
    float pA = SA / S;
    float pB = SB / S;

    return pA*nA+pB*nB;
}
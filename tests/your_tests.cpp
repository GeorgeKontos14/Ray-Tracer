// Put your includes here
#include "bvh.h"
#include "render.h"
#include "sampler.h"
#include "scene.h"
#include "shading.h"
#include "extra.h"
#include <limits>

// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <catch2/catch_all.hpp>
#include <glm/glm.hpp>
DISABLE_WARNINGS_POP()

// In this file you can add your own unit tests using the Catch2 library.
// You can find the documentation of Catch2 at the following link:
// https://github.com/catchorg/Catch2/blob/devel/docs/assertions.md
//
// These tests are only to help you verify that your code is correct.
// You don't have to hand them in; we will not consider them when grading.
//

// Add your tests here, if you want :D
TEST_CASE("StudentTest")
{
    // Add your own tests here...
}

TEST_CASE("BVHTesting")
{
    Features features = {
        .enableAccelStructure = true
    };
    Scene scene = loadScenePrebuilt(SceneType::CornellBox, DATA_DIR);
    BVH bvh(scene, features);
    RenderState state = { .scene = scene, .features = features, .bvh = bvh, .sampler = {} };
    Vertex v0 = {
        glm::vec3 { 0, 0, 0 },
        glm::vec3 { 1, 0, 0 },
        glm::vec3 { 0, 0, 0 }
    };
    Vertex v1 = {
        glm::vec3 { 1, 0, 0 },
        glm::vec3 { 1, 0, 0 },
        glm::vec3 { 0, 0, 0 }
    };
    Vertex v2 = {
        glm::vec3 { 0, 1, 0 },
        glm::vec3 { 1, 0, 0 },
        glm::vec3 { 0, 0, 0 }
    };
    Vertex v3 = {
        glm::vec3 { 2, 2, 2 },
        glm::vec3 { 1, 0, 0 },
        glm::vec3 { 0, 0, 0 }
    };
    Vertex v4 = {
        glm::vec3 { 3, 2, 2 },
        glm::vec3 { 1, 0, 0 },
        glm::vec3 { 0, 0, 0 }
    };
    Vertex v5 = {
        glm::vec3 { 3, 4, 1 },
        glm::vec3 { 1, 0, 0 },
        glm::vec3 { 0, 0, 0 }
    };
    Vertex v6 = {
        glm::vec3 { -2, 1, 7 },
        glm::vec3 { 1, 0, 0 },
        glm::vec3 { 0, 0, 0 }
    };
    Vertex v7 = {
        glm::vec3 { 4, -5, 3 },
        glm::vec3 { 1, 0, 0 },
        glm::vec3 { 0, 0, 0 }
    };
    Vertex v8 = {
        glm::vec3 { 0, 1, -8 },
        glm::vec3 { 1, 0, 0 },
        glm::vec3 { 0, 0, 0 }
    };
    Vertex v9 = {
        glm::vec3 { -4, 0, 6 },
        glm::vec3 { 1, 0, 0 },
        glm::vec3 { 0, 0, 0 }
    };
    BVHInterface::Primitive p0 = {
        0,
        v0,
        v1,
        v2
    };
    BVHInterface::Primitive p1 = {
        1,
        v0,
        v1,
        v3
    };
    BVHInterface::Primitive p2 = {
        2,
        v4,
        v1,
        v2
    };
    BVHInterface::Primitive p3 = {
        3,
        v0,
        v1,
        v4
    };
    BVHInterface::Primitive p4 = {
        4,
        v3,
        v4,
        v5
    };
    BVHInterface::Primitive p5 = {
        5,
        v4,
        v5,
        v6
    };
    BVHInterface::Primitive p6 = {
        6,
        v0,
        v1,
        v6
    };
    BVHInterface::Primitive p7 = {
        7,
        v7,
        v8,
        v9
    };
    PrimitiveData d0 = {
        p0,
        computePrimitiveCentroid(p0)
    };
    PrimitiveData d1 = {
        p1,
        computePrimitiveCentroid(p1)
    };
    PrimitiveData d2 = {
        p2,
        computePrimitiveCentroid(p2)
    };
    PrimitiveData d3 = {
        p3,
        computePrimitiveCentroid(p3)
    };
    PrimitiveData d4 = {
        p4,
        computePrimitiveCentroid(p4)
    };
    PrimitiveData d5 = {
        p5,
        computePrimitiveCentroid(p5)
    };
    PrimitiveData d6 = {
        p6,
        computePrimitiveCentroid(p6)
    };
    PrimitiveData d7 = {
        p7,
        computePrimitiveCentroid(p7)
    };
    SECTION("Primitive AABB")
    {
        AxisAlignedBox aabb1 = computePrimitiveAABB(p1);
        AxisAlignedBox aabb3 = computePrimitiveAABB(p3);
        AxisAlignedBox aabb7 = computePrimitiveAABB(p7);
        CHECK((aabb1.lower == glm::vec3 { 0, 0, 0 } && aabb1.upper == glm::vec3 { 2, 2, 2 }));
        CHECK((aabb3.lower == glm::vec3 { 0, 0, 0 } && aabb3.upper == glm::vec3 { 3, 2, 2 }));
        CHECK((aabb7.lower == glm::vec3 { -4, -5, -8 } && aabb7.upper == glm::vec3 { 4, 1, 6 }));
    }

    SECTION("Span AABB")
    {
        BVHInterface::Primitive arr012[] = { p0, p1, p2 };
        AxisAlignedBox aabb012 = computeSpanAABB(arr012);
        BVHInterface::Primitive arr456[] = { p4, p5, p6 };
        AxisAlignedBox aabb456 = computeSpanAABB(arr456);
        BVHInterface::Primitive arr12[] = { p1, p2 };
        AxisAlignedBox aabb12 = computeSpanAABB(arr12);
        BVHInterface::Primitive arr0123[] = { p0, p1, p2, p3 };
        AxisAlignedBox aabb0123 = computeSpanAABB(arr0123);
        BVHInterface::Primitive arr234567[] = { p2, p3, p4, p5, p6, p7 };
        AxisAlignedBox aabb234567 = computeSpanAABB(arr234567);
        CHECK((aabb012.lower == glm::vec3 { 0, 0, 0 } && aabb012.upper == glm::vec3 { 3, 2, 2 }));
        CHECK((aabb456.lower == glm::vec3 { -2, 0, 0 } && aabb456.upper == glm::vec3 { 3, 4, 7 }));
        CHECK((aabb12.lower == glm::vec3 { 0, 0, 0 } && aabb12.upper == glm::vec3 { 3, 2, 2 }));
        CHECK((aabb0123.lower == glm::vec3 { 0, 0, 0 } && aabb0123.upper == glm::vec3 { 3, 2, 2 }));
        CHECK((aabb234567.lower == glm::vec3 { -4, -5, -8 } && aabb234567.upper == glm::vec3 { 4, 4, 7 }));
    }

    SECTION("Primitive Centroid")
    {
        CHECK(computePrimitiveCentroid(p0) == glm::vec3 { (float)1 / 3, (float)1 / 3, 0 });
        CHECK(computePrimitiveCentroid(p2) == glm::vec3 { (float)4 / 3, 1, (float)2 / 3 });
        CHECK(computePrimitiveCentroid(p5) == glm::vec3 { (float)4 / 3, (float)7 / 3, (float)10 / 3 });
        CHECK(computePrimitiveCentroid(p7) == glm::vec3 { 0, (float)-4 / 3, (float)1 / 3 });
    }

    SECTION("AABB Longest Axis")
    {
        AxisAlignedBox aabb0 = {
            .lower = { 0, 0, 0 },
            .upper = { 3, 5, 7 }
        };
        AxisAlignedBox aabb1 = {
            .lower = { -1, 2, 8 },
            .upper = { 3, 4, 9 }
        };
        AxisAlignedBox aabb2 = {
            .lower = { 1, 2, 3 },
            .upper = { 3, 9, 6 }
        };
        AxisAlignedBox aabb3 = {
            .lower = { 1, 1, 1 },
            .upper = { 3, 3, 3 }
        };
        AxisAlignedBox aabb4 = {
            .lower = { 1, 1, 2 },
            .upper = { 5, 3, 6 }
        };
        AxisAlignedBox aabb5 = {
            .lower = { 0, 1, 2 },
            .upper = { 4, 5, 4 }
        };
        AxisAlignedBox aabb6 = {
            .lower = { 0, 0, 0 },
            .upper = { 0.5, 11.3, 11.3 }
        };
        CHECK(computeAABBLongestAxis(aabb0) == 2);
        CHECK(computeAABBLongestAxis(aabb1) == 0);
        CHECK(computeAABBLongestAxis(aabb2) == 1);
        CHECK(computeAABBLongestAxis(aabb3) == 0);
        CHECK(computeAABBLongestAxis(aabb4) == 0);
        CHECK(computeAABBLongestAxis(aabb5) == 0);
        CHECK(computeAABBLongestAxis(aabb6) == 1);
    }

    SECTION("Compare Primitives")
    {
        CHECK(comparePrimitives(0, d0, d0) == 0);
        CHECK(comparePrimitives(0, d2, d3) == 0);
        CHECK(comparePrimitives(1, d1, d3) == 0);
        CHECK(comparePrimitives(2, d1, d2) == 0);
        CHECK(comparePrimitives(0, d0, d6) == 1);
        CHECK(comparePrimitives(0, d0, d1) == -1);
        CHECK(comparePrimitives(1, d1, d7) == 1);
        CHECK(comparePrimitives(1, d2, d4) == -1);
        CHECK(comparePrimitives(2, d4, d0) == 1);
        CHECK(comparePrimitives(2, d2, d5) == -1);
    }

    SECTION("Merge Vectors")
    {
        std::vector<PrimitiveData> vec0x = { d6, d0, d4 };
        std::vector<PrimitiveData> vec1x = { d7, d1, d2 };
        std::vector<PrimitiveData> mergedx = { d6, d7, d0, d1, d2, d4 };
        std::vector<PrimitiveData> resx = mergeVectors(0, vec0x, vec1x);
        CHECK(equalVectors(mergedx, resx));
        std::vector<PrimitiveData> vec0y = { d0, d3, d5 };
        std::vector<PrimitiveData> vec1y = { d7, d6, d1, d2, d4 };
        std::vector<PrimitiveData> mergedy = { d7, d6, d0, d1, d3, d2, d5, d4 };
        std::vector<PrimitiveData> resy = mergeVectors(1, vec0y, vec1y);
        CHECK(equalVectors(mergedy, resy));
        std::vector<PrimitiveData> vec0z = { d0, d2, d5 };
        std::vector<PrimitiveData> vec1z = { d7, d3, d1, d4, d6 };
        std::vector<PrimitiveData> mergedz = { d0, d7, d3, d1, d2, d4, d6, d5 };
        std::vector<PrimitiveData> resz = mergeVectors(2, vec0z, vec1z);
        CHECK(equalVectors(mergedz, resz));
    }

    SECTION("Sort Primitive Data")
    {
        std::vector<PrimitiveData> data = { d0, d1, d2, d3, d4, d5, d6, d7 };
        std::vector<PrimitiveData> sortedX = { d6, d7, d0, d1, d5, d3, d2, d4 };
        std::vector<PrimitiveData> sortedY = { d7, d0, d6, d3, d1, d2, d5, d4 };
        std::vector<PrimitiveData> sortedZ = { d0, d7, d2, d1, d3, d4, d6, d5 };
        sortPrimitiveData(0, data);
        CHECK(equalVectors(data, sortedX));
        sortPrimitiveData(1, data);
        CHECK(equalVectors(data, sortedY));
        sortPrimitiveData(2, data);
        CHECK(equalVectors(data, sortedZ));
    }

    SECTION("Split Primitives By Median")
    {
        std::vector<BVHInterface::Primitive> primitivesVec = { p0, p1, p2, p3, p4, p5, p6, p7 };
        std::vector<BVHInterface::Primitive> sortedXVec = { p6, p7, p0, p1, p5, p3, p2, p4 };
        std::span<BVHInterface::Primitive> primitives(primitivesVec);
        std::span<BVHInterface::Primitive> sortedX(sortedXVec);
        AxisAlignedBox aabb = computeSpanAABB(primitives);
        int splitInd = splitPrimitivesByMedian(aabb, 0, primitives);
        CHECK(equalSpans(primitives, sortedX));
        CHECK(splitInd == 4);
        std::vector<BVHInterface::Primitive> primitivesVec1 = { p0, p1, p2, p3, p4, p5, p6 };
        std::vector<BVHInterface::Primitive> sortedXVec1 = { p6, p0, p1, p5, p3, p2, p4 };
        std::span<BVHInterface::Primitive> primitives1(primitivesVec1);
        std::span<BVHInterface::Primitive> sortedX1(sortedXVec1);
        AxisAlignedBox aabb1 = computeSpanAABB(primitives1);
        int splitInd1 = splitPrimitivesByMedian(aabb1, 0, primitives1);
        CHECK(equalSpans(primitives1, sortedX1));
        CHECK(splitInd1 == 4);
    }

    SECTION("Split Primitives SAH")
    {
        std::vector<BVHInterface::Primitive> primitivesVec = { p0, p1, p2, p3, p4, p5, p6, p7 };
        std::span<BVHInterface::Primitive> primitives(primitivesVec);
        uint32_t ind = splitPrimitivesBySAHBin(computeSpanAABB(primitives), 0, primitives);
        CHECK(ind != 2);
    }
}

TEST_CASE("BloomFilterTesting")
{
    SECTION("Factorial")
    {
        CHECK(factorial(0) == 1);
        CHECK(factorial(1) == 1);
        CHECK(factorial(-7777) == 1);
        CHECK(factorial(2) == 2);
        CHECK(factorial(5) == 120);
    }
     
    SECTION()
    // The below tests are not "good" unit tests. They don't actually test correctness.
// They simply exist for demonstrative purposes. As they interact with the interfaces
// (scene, bvh_interface, etc), they allow you to verify that you haven't broken
// our grading interface. They should compile without changes. If they do
// not compile, neither will our grading tests!
TEST_CASE("InterfaceTest")
{
    // Setup a RenderState object with some defaults
    Features features = {
        .enableShading = true,
        .enableAccelStructure = false, // BVH is not actually active r.n.
        .shadingModel = ShadingModel::Lambertian
    };
    Scene scene = loadScenePrebuilt(SceneType::CornellBox, DATA_DIR);
    BVH bvh(scene, features);
    RenderState state = { .scene = scene, .features = features, .bvh = bvh, .sampler = {} };

    SECTION("BVH generation")
    {
        // There's something in here?
        CHECK(!state.bvh.primitives().empty());
    }

    SECTION("BVH traversal")
    {
        Ray ray = { .origin = glm::vec3(0), .direction = glm::vec3(1) };
        HitInfo hitInfo;

        // Hit something?
        CHECK(state.bvh.intersect(state, ray, hitInfo));
        CHECK(ray.t != std::numeric_limits<float>::max());
    }

    SECTION("Hit shading")
    {
        Ray ray = { .origin = glm::vec3(0), .direction = glm::vec3(1) };
        HitInfo hitInfo;
        state.bvh.intersect(state, ray, hitInfo);

        // Shaded something?
        glm::vec3 Lo = computeShading(state, ray.direction, -ray.direction, glm::vec3(1), hitInfo);
        CHECK(glm::any(glm::notEqual(Lo, glm::vec3(0))));
    }
}

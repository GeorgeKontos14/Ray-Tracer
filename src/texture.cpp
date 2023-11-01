#include "texture.h"
#include "render.h"
#include <framework/image.h>

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// the nearest texel to the coordinates is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the nearest corresponding texel
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureNearest(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)
    int width = image.width;
    int height = image.height;
    float x = texCoord.x;
    float y = texCoord.y;
    float j = glm::floor((width * x));
    float i = glm::floor((height * y));
    int j1 = glm::floor(j);
    int i1 = glm::floor(i);
    //they have to be between bounds
    j1 = glm::clamp(j1, 0, width - 1);
    i1 = glm::clamp(i1, 0, height - 1);
    //j*width+i
    int index = j1 * width + i1;
    
    return image.pixels[index];
}

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// a bilinearly interpolated texel is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the filter of the corresponding texels
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureBilinear(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)
    int width = image.width;
    int height = image.height;
    float x = texCoord.x;
    float y = texCoord.y;

    float i = y * height + 0.5f;
    float j = width * x + 0.5f;

    //we need the fractional and integer parts for the formula
    
    int iint = glm::floor(i);
    float ifractional = i - iint;
    int jint = glm::floor(j);
    float jfractional = j - jint;
    float widt4 = (1 - jfractional) * (1 - ifractional);
    float widt3 = (1 - jfractional) * ifractional;
    float widt2 = jfractional * (1 - ifractional);
    float widt1 = jfractional * ifractional;
    glm::vec3 t4 = image.pixels[width * (jint-1)+iint-1];
    glm::vec3 t3 = image.pixels[width * (jint - 1) + iint];
    glm::vec3 t2 = image.pixels[width * jint + iint - 1];
    glm::vec3 t1 = image.pixels[jint * width + iint];
    glm::vec3 res = widt4*t4+widt3*t3+widt2*t2+widt1*t1;

    //bilinear interpolation
    return res;

}
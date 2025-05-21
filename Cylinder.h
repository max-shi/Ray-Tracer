/*----------------------------------------------------------
* COSC363  Ray Tracer Assignment
*
*  Cylinder class
*  This is a subclass of SceneObject, and hence implements the
*  methods intersect() and normal().
-------------------------------------------------------------*/

#ifndef H_CYLINDER
#define H_CYLINDER

#include <glm/glm.hpp>
#include "SceneObject.h"
#include "TextureBMP.h"

/**
 * Class that defines cylinder object with caps
 */
class Cylinder : public SceneObject
{

private:
    glm::vec3 center;      // Center of the base of the cylinder
    float radius;          // Radius of the cylinder
    float height;          // Height of the cylinder
    bool hasCap;           // Whether the cylinder has a cap
    bool isTextured_;      // Whether the cylinder is textured
    TextureBMP* texture_;  // Texture for the cylinder

public:
    Cylinder() : center(glm::vec3(0)), radius(1), height(1), hasCap(true), isTextured_(false), texture_(nullptr) {}

    Cylinder(glm::vec3 cen, float rad, float hgt, glm::vec3 col, bool cap = true)
        : center(cen), radius(rad), height(hgt), hasCap(cap), isTextured_(false), texture_(nullptr)
    {
        setColor(col);
    }

    ~Cylinder() {
        // No need to delete texture_ as it's managed elsewhere
    }

    float intersect(glm::vec3 pos, glm::vec3 dir);

    glm::vec3 normal(glm::vec3 p);
    
    // Method to set texture
    void setTexture(TextureBMP* tex) {
        texture_ = tex;
        isTextured_ = (tex != nullptr);
    }
    
    // Method to get texture coordinates for a point on the cylinder
    glm::vec2 getTexCoord(glm::vec3 p);
    
    // Method to get color at a point (either texture or material color)
    glm::vec3 getColorAt(glm::vec3 p);
};

#endif //!H_CYLINDER

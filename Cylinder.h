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
    glm::vec3 center;
    float radius;
    float height;
    bool hasCap;
    bool isTextured_;
    TextureBMP* texture_;

public:
    Cylinder() : center(glm::vec3(0)), radius(1), height(1), hasCap(true), isTextured_(false), texture_(nullptr) {}

    Cylinder(glm::vec3 cen, float rad, float hgt, glm::vec3 col, bool cap = true)
        : center(cen), radius(rad), height(hgt), hasCap(cap), isTextured_(false), texture_(nullptr)
    {
        setColor(col);
    }

    ~Cylinder() {
    }

    float intersect(glm::vec3 pos, glm::vec3 dir);

    glm::vec3 normal(glm::vec3 p);
    
    // Method to set texture
    void setTexture(TextureBMP* tex) {
        texture_ = tex;
        isTextured_ = (tex != nullptr);
    }
    glm::vec2 getTexCoord(glm::vec3 p);
    glm::vec3 getColorAt(glm::vec3 p);
};

#endif //!H_CYLINDER

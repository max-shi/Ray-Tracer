#ifndef TORUS_H
#define TORUS_H

#include "SceneObject.h"
#include <complex>
#include <vector>
#include <glm/gtc/matrix_transform.hpp>

class Torus : public SceneObject {
private:
    float majorRadius;
    float minorRadius;
    glm::vec3 center;
    glm::mat4 transform;    // Transformation matrix
    glm::mat4 invTransform; // Inverse transformation matrix

    std::vector<float> solveQuartic(const std::vector<double>& coeffs) const;

public:
    Torus(glm::vec3 center, float majorRadius, float minorRadius);

    // Transformation methods
    void rotate(float angle, glm::vec3 axis);
    void translate(glm::vec3 translation);
    void scale(glm::vec3 scale);
    
    float intersect(glm::vec3 p0, glm::vec3 dir) override;
    glm::vec3 normal(glm::vec3 p) override;
};

#endif
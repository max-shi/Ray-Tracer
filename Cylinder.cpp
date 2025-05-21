/*----------------------------------------------------------
* COSC363  Ray Tracer Assignment
*
*  Cylinder class
*  This is a subclass of SceneObject, and hence implements the
*  methods intersect() and normal().
-------------------------------------------------------------*/

#include "Cylinder.h"
#include <math.h>
#include <algorithm>

/**
* Implementing the intersection check equation for a cylinder with caps
*/
float Cylinder::intersect(glm::vec3 pos, glm::vec3 dir) {
    float t = -1.0;  // Default return value (no intersection)
    
    // Intersection with the cylindrical surface
    glm::vec3 d = pos - center;
    
    // Coefficients for the quadratic equation
    float a = powf(dir.x, 2) + powf(dir.z, 2);
    float b = 2*(d.x * dir.x + d.z * dir.z);
    float c = powf(d.x, 2) + powf(d.z, 2) - powf(radius, 2);
    float delta = powf(b, 2) - 4*(a * c);
    
    // Check if ray intersects the cylindrical surface
    if(delta >= 0.0) {
        float t1 = (-b - sqrt(delta))/(2*a);
        float t2 = (-b + sqrt(delta))/(2*a);
        
        // Ensure t1 <= t2
        if (t1 > t2) {
            std::swap(t1, t2);
        }
        
        // Check if the intersection points are within the cylinder's height
        float y1 = pos.y + dir.y * t1;
        float y2 = pos.y + dir.y * t2;
        
        // Valid intersection points must be within the cylinder's height range
        bool valid1 = (t1 > 0.01) && (y1 >= center.y) && (y1 <= center.y + height);
        bool valid2 = (t2 > 0.01) && (y2 >= center.y) && (y2 <= center.y + height);
        
        if (valid1) t = t1;
        else if (valid2) t = t2;
    }
    
    // If hasCap is true, check for intersection with the caps
    if (hasCap) {
        // Bottom cap (at center.y)
        if (fabs(dir.y) > 1e-6) {  // Avoid division by zero
            float tBottom = (center.y - pos.y) / dir.y;
            if (tBottom > 0.01) {
                glm::vec3 p = pos + tBottom * dir;
                float dist = powf(p.x - center.x, 2) + powf(p.z - center.z, 2);
                if (dist <= powf(radius, 2)) {
                    if (t < 0 || tBottom < t) t = tBottom;
                }
            }
        }
        
        // Top cap (at center.y + height)
        if (fabs(dir.y) > 1e-6) {  // Avoid division by zero
            float tTop = (center.y + height - pos.y) / dir.y;
            if (tTop > 0.01) {
                glm::vec3 p = pos + tTop * dir;
                float dist = powf(p.x - center.x, 2) + powf(p.z - center.z, 2);
                if (dist <= powf(radius, 2)) {
                    if (t < 0 || tTop < t) t = tTop;
                }
            }
        }
    }
    
    return t;
}

/**
* Implementing normal calculation for cylinder with caps
*/
glm::vec3 Cylinder::normal(glm::vec3 p) {
    // Check if the point is on one of the caps
    if (hasCap) {
        // Bottom cap
        if (fabs(p.y - center.y) < 1e-4) {
            return glm::vec3(0, -1, 0);
        }
        // Top cap
        if (fabs(p.y - (center.y + height)) < 1e-4) {
            return glm::vec3(0, 1, 0);
        }
    }
    
    // Point is on the cylindrical surface
    glm::vec3 n(p.x - center.x, 0, p.z - center.z);
    return glm::normalize(n);
}

/**
* Calculate texture coordinates for a point on the cylinder
* s = angle around the cylinder (0 to 1)
* t = height along the cylinder (0 to 1)
*/
glm::vec2 Cylinder::getTexCoord(glm::vec3 p) {
    // For points on the cylindrical surface
    float s, t;
    
    // Calculate the angle around the cylinder (for s coordinate)
    float dx = p.x - center.x;
    float dz = p.z - center.z;
    float angle = atan2(dz, dx);
    
    // Convert angle to [0,1] range for s coordinate
    s = angle / (2.0f * M_PI);
    if (s < 0) s += 1.0f;  // Ensure s is in [0,1]
    
    // Calculate t coordinate based on height
    t = (p.y - center.y) / height;
    
    // Ensure t is in [0,1] range
    t = glm::clamp(t, 0.0f, 1.0f);
    
    return glm::vec2(s, t);
}

/**
* Get color at a point on the cylinder (either from texture or material color)
*/
glm::vec3 Cylinder::getColorAt(glm::vec3 p) {
    if (!isTextured_ || texture_ == nullptr) {
        return color_; // Return the material color if not textured
    }
    
    // Get texture coordinates
    glm::vec2 texCoord = getTexCoord(p);
    
    // Get color from texture
    return texture_->getColorAt(texCoord.x, texCoord.y);
}

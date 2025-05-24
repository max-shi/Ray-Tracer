/*==================================================================================
* COSC 363  Computer Graphics
* Ray Tracer for Cosc363 Assignment 2 (msh254)
*===================================================================================
*/
#include <iostream>
#include <cmath>
#include "Plane.h"
#include <vector>
#include <glm/glm.hpp>
#include "Sphere.h"
#include "SceneObject.h"
#include "Ray.h"
#include <GL/freeglut.h>
#include "TextureBMP.h"
#include "Torus.h"
#include "Cylinder.h"
#include <random>
#include <chrono>
using namespace std;

TextureBMP texture;
TextureBMP cylinderTexture("../fabric-pattern-polyhaven.bmp");

// These below parameters affect the speed at which the ray tracer loads
const int NUMDIV = 500;
const int MAX_STEPS = 10;
bool antiAliasingEnabled = false;
bool stochasticSamplingEnabled = false;
int SAMPLES_PER_PIXEL = 4;
const float EDIST = 40;
const float XMIN = -20.0;
const float XMAX = 20.0;
const float YMIN = -20.0;
const float YMAX = 20.0;
const int MAX_ADAPTIVE_DEPTH = 1;
float COLOR_THRESHOLD = 0.1f;
float LIGHT_RADIUS = 3.0f;  // Radius of the area light for soft shadows

// Random number generator for stochastic sampling
std::mt19937 rng(std::chrono::high_resolution_clock::now().time_since_epoch().count());
std::uniform_real_distribution<float> dist(0.0f, 1.0f);
vector<SceneObject*> sceneObjects;

//-------------------- Trace Function --------------------------
glm::vec3 trace(Ray ray, int step) {
    glm::vec3 backgroundCol(0);
    glm::vec3 lightPos(0, 40, -90);
    glm::vec3 ambientLight(0.2, 0.2, 0.2);
    glm::vec3 color(0);
    SceneObject* obj;

    ray.closestPt(sceneObjects);
    // The below line should never be actually called as we have a bounding box, but good for a check i guess..
    if(ray.index == -1) return backgroundCol;
    obj = sceneObjects[ray.index];
    
    // Checkered floor computations. As the floor is the first object in the scene, we can check for it with ray.index == 0
    if (ray.index == 0) {
        int squareSize = 10;
        int ix = floor(ray.hit.x / squareSize);
        int iz = floor(ray.hit.z / squareSize);
        if ((ix + iz) % 2 == 0) {
            color = glm::vec3(1, 1, 1);
        } else {
            color = glm::vec3(0, 0, 0);
        }
        obj->setColor(color);
    }
    
    // Method only for cylinder, as it has a texture mapped to it
    Cylinder* cylinderObj = dynamic_cast<Cylinder*>(obj);
    if (cylinderObj != nullptr) {
        color = cylinderObj->getColorAt(ray.hit);
        obj->setColor(color);
    }

    glm::vec3 LightVec = lightPos - ray.hit;
    float lightDist = glm::length(LightVec);

    /* see this https://learnwebgl.brown37.net/09_lights/lights_attenuation.html#:~:text=In%20the%20physical%20world%20the,be%20proportional%20to%201%2Fd.
     * Main points : in real life light attenuation is 1/d^2 : d being distance
     * However this causes light dropoff very fast, which is why it is common to  use 1/d, but in context of a ray tracer
     *  " In the original OpenGL lighting model, the equation 1.0/(c1 + c2*d + c3*d^2) was used to give programmers control over attenuation. " (from source)
     * hence the attenuation.
     * Note: the values for c2, c3 are small because of the scale of the room.. (the room is 160x160 units) : which was an oversight when creating the program.
     */

    float attenuation = 1.0f / (1.0f + 0.002f * lightDist + 0.0004f * lightDist * lightDist);
    float lightIntensity = 1.2f;

    //-------------------- Stochastic Sampling --------------------------
    if (stochasticSamplingEnabled) { // Check if stochastic sampling is enabled for soft shadows
        // see https://www.youtube.com/watch?v=NCptEJ1Uevg&ab_channel=OGLDEV
        // https://cg.ivd.kit.edu/publications/2015/sssm/StochasticSoftShadows.pdf
        // Soft shadows with area light sampling

        // Setup; normalise light direction and find u,v vectors
        glm::vec3 lightDir = glm::normalize(LightVec);
        glm::vec3 u = glm::normalize(glm::cross(lightDir, glm::vec3(0, 1, 0)));
        if (glm::length(u) < 0.1f) u = glm::normalize(glm::cross(lightDir, glm::vec3(1, 0, 0)));
        glm::vec3 v = glm::cross(u, lightDir);
        
        // Calculate lighting with soft shadows
        glm::vec3 diffuseColor(0);
        glm::vec3 specularColor(0);
        // Counter for samples that are in shadow
        int shadowCount = 0;

        // Loop through multiple samples for the area light
        for (int i = 0; i < SAMPLES_PER_PIXEL; i++) {
            // Generate random point on area light
            float rx = (dist(rng) * 2.0f - 1.0f) * LIGHT_RADIUS;
            float ry = (dist(rng) * 2.0f - 1.0f) * LIGHT_RADIUS;
            // Calculate position on area light using the random offsets
            glm::vec3 randomLightPos = lightPos + rx * u + ry * v;
            
            // Calculate lighting from this sample point
            glm::vec3 sampleLightVec = randomLightPos - ray.hit;
            float sampleLightDist = glm::length(sampleLightVec);
            
            // Check if this light sample is visible (shadow ray)
            Ray shadowRay(ray.hit, sampleLightVec);
            // Find closest intersection with scene objects
            shadowRay.closestPt(sceneObjects);

            if (shadowRay.index > -1 && shadowRay.dist < sampleLightDist) {
                shadowCount++;
                continue;
            }
            
            // Add contribution from this light sample
            glm::vec3 sampleColor = obj->lighting(randomLightPos, -ray.dir, ray.hit);
            diffuseColor += sampleColor;
        }
        
        // Average the results
        // Calculate percentage of samples not in shadow
        float shadowFactor = 1.0f - (float)shadowCount / (float)SAMPLES_PER_PIXEL;
        if (shadowFactor > 0) {
            color = (diffuseColor / (float)SAMPLES_PER_PIXEL) * attenuation * lightIntensity * shadowFactor; // Average color with shadow factor
        } else { // If all samples are in shadow
            color = glm::vec3(0); // Set color to black (full shadow)
        }
    } else {
    //-------------------- Non Stochastic Sampling --------------------------
        // Normal flow; shadow is black
        color = obj->lighting(lightPos, -ray.dir, ray.hit) * attenuation * lightIntensity;
        Ray shadowRay(ray.hit, LightVec);
        shadowRay.closestPt(sceneObjects);
        if (shadowRay.index > -1 && shadowRay.dist < lightDist) {
            color = glm::vec3(0);
        }
    }
    //-------------------- Ambient light Calculations --------------------------
    color += ambientLight * obj->getColor();

    //-------------------- Transparency Ray + Calculations --------------------------
    if (obj->isTransparent() && step < MAX_STEPS) {
        float t = obj->getTransparencyCoeff();
        // Create a secondary ray along the direction of transmission of light
        Ray transmissionRay(ray.hit, ray.dir);
        glm::vec3 transmittedColor = trace(transmissionRay, step + 1);
        // Add the scaled transmitted color to the pixel color (I = IA + t*IC)
        color = color + (t * transmittedColor);
    }
    
    //-------------------- Refractive Ray + Calculations --------------------------
    if (obj->isRefractive() && step < MAX_STEPS) {
        const float ETA = obj->getRefractiveIndex(); // Refractive index
        float transVal = obj->getRefractionCoeff(); // Refraction coefficient
        glm::vec3 normalVec = obj->normal(ray.hit);
        
        // First refraction (entering the object)
        glm::vec3 refractedDir = glm::refract(ray.dir, normalVec, 1.0f/ETA);
        Ray refractedRay(ray.hit, refractedDir);
        refractedRay.closestPt(sceneObjects);
        
        if (refractedRay.index != -1) {
            // Second refraction (exiting the object)
            glm::vec3 exitNormal = -sceneObjects[refractedRay.index]->normal(refractedRay.hit);
            glm::vec3 exitDir = glm::refract(refractedDir, exitNormal, ETA);
            Ray exitRay(refractedRay.hit, exitDir);
            
            glm::vec3 refractedColor = trace(exitRay, step + 1);
            
            // Blend original color with refracted color based on refraction coefficient
            color = color * (1.0f - transVal) + refractedColor * transVal;
        }
    }

    //-------------------- Reflective Ray + Calculations --------------------------
    if (obj->isReflective() && step < MAX_STEPS) {
        float rho = obj->getReflectionCoeff();
        glm::vec3 normalVec = obj->normal(ray.hit);
        
        // Standard perfect mirror reflection
        // TODO rough reflections?
        glm::vec3 reflectedDir = glm::reflect(ray.dir, normalVec);
        Ray reflectedRay(ray.hit, reflectedDir);
        glm::vec3 reflectedColor = trace(reflectedRay, step + 1);
        color = color + (rho * reflectedColor);
    }

    color = glm::clamp(color, 0.0f, 1.0f);
    return color;
}

/**-------------------- Helper Function For Adaptive Sampling --------------------------
 * Determines if two colors are sufficiently different to require subdivision in adaptive sampling.
 * This function calculates the distance between two color vectors and compares
 * it with the COLOR_THRESHOLD to decide if further sampling is needed.
 */
bool needsSubdivision(const glm::vec3& col1, const glm::vec3& col2) {
    float diff = glm::length(col1 - col2);
    return diff > COLOR_THRESHOLD;
}

/**-------------------- Adaptive Sampling --------------------------
 * Recursively samples a pixel region using adaptive subdivision for anti-aliasing.
 * This function implements adaptive supersampling by recursively subdividing pixel regions
 * where color variation is high. It samples four sub-regions of the current region and
 * determines whether further subdivision is needed based on color differences.
 * Note: this is never called directly by the ray tracer: see calculatePixelColor().
 */
glm::vec3 adaptiveSample(float x, float y, float width, float height, const glm::vec3& eye,
                         int maxDepth, int currentDepth, const glm::vec3& parentColor) {
    if (currentDepth >= maxDepth) {
        glm::vec3 dir(x + width/2, y + height/2, -EDIST);
        Ray ray(eye, dir);
        return trace(ray, 1);
    }

    // Sample four sub-regions
    glm::vec3 samples[4];
    samples[0] = adaptiveSample(x, y, width/2, height/2, eye, maxDepth, currentDepth+1, parentColor);
    samples[1] = adaptiveSample(x + width/2, y, width/2, height/2, eye, maxDepth, currentDepth+1, parentColor);
    samples[2] = adaptiveSample(x, y + height/2, width/2, height/2, eye, maxDepth, currentDepth+1, parentColor);
    samples[3] = adaptiveSample(x + width/2, y + height/2, width/2, height/2, eye, maxDepth, currentDepth+1, parentColor);

    // Check if subdivision is needed : comparison with previous parentColor (this does not compare between samples -> only with the previous color)
    bool subdivide = false;
    if (currentDepth > 0) {
        for (int i = 0; i < 4; i++) {
            if (needsSubdivision(samples[i], parentColor)) {
                subdivide = true;
                break;
            }
        }
    }

    if (subdivide && currentDepth < maxDepth-1) {
        // Average of samples
        return (samples[0] + samples[1] + samples[2] + samples[3]) / 4.0f;
    } else {
        // Check if subdivision is needed : comparison with samples
        bool similar = true;
        for (int i = 1; i < 4; i++) {
            if (needsSubdivision(samples[i], samples[0])) {
                similar = false;
                break;
            }
        }

        if (similar) {
            // if all similar, average the colors
            return (samples[0] + samples[1] + samples[2] + samples[3]) / 4.0f;
        } else if (currentDepth < maxDepth-1) {
            // if NOT all similar, recurse again
            // FIXME: currently max adaptive depth is set to 1 meaning this actually will NEVER get called.
            return adaptiveSample(x, y, width, height, eye, maxDepth, currentDepth+1, samples[0]);
        } else {
            // this is the case where we reach the maxDepth (then we just take the average) -> no recursion
            return (samples[0] + samples[1] + samples[2] + samples[3]) / 4.0f;
        }
    }
}

/** -------------------- Calculate Pixel Color --------------------------
* Calculates the color for a pixel based on ray tracing.
* This function either performs a single ray trace for the pixel (when anti-aliasing is disabled)
* or uses adaptive sampling for higher quality anti-aliasing when enabled.
*/
glm::vec3 calculatePixelColor(float xp, float yp, float cellX, float cellY, const glm::vec3& eye) {
    if (!antiAliasingEnabled) {
        glm::vec3 dir(xp + 0.5f * cellX, yp + 0.5f * cellY, -EDIST);
        Ray ray(eye, dir);
        return trace(ray, 1);
    } else {
        glm::vec3 initialColor(0);
        return adaptiveSample(xp, yp, cellX, cellY, eye, MAX_ADAPTIVE_DEPTH, 0, initialColor);
    }
}

void display() {
    float xp, yp;
    float cellX = (XMAX - XMIN) / NUMDIV;
    float cellY = (YMAX - YMIN) / NUMDIV;
    glm::vec3 eye(0., 0., 0.);

    glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glBegin(GL_QUADS);

    for (int i = 0; i < NUMDIV; i++) {
        xp = XMIN + i * cellX;
        for (int j = 0; j < NUMDIV; j++) {
            yp = YMIN + j * cellY;

            glm::vec3 col = calculatePixelColor(xp, yp, cellX, cellY, eye);
            glColor3f(col.r, col.g, col.b);
            glVertex2f(xp, yp);
            glVertex2f(xp + cellX, yp);
            glVertex2f(xp + cellX, yp + cellY);
            glVertex2f(xp, yp + cellY);
        }
    }

    glEnd();
    glFlush();
}

void keyboard(unsigned char key, int x, int y) {
    if (key == 'a' || key == 'A') {
        antiAliasingEnabled = !antiAliasingEnabled;
        printf("Anti-aliasing %s\n", antiAliasingEnabled ? "ENABLED" : "DISABLED");
        glutPostRedisplay();
    }
    else if (key == 's' || key == 'S') {
        stochasticSamplingEnabled = !stochasticSamplingEnabled;
        printf("Stochastic sampling %s\n", stochasticSamplingEnabled ? "ENABLED" : "DISABLED");
        glutPostRedisplay();
    }
    else if (key == '+') {
        COLOR_THRESHOLD += 0.02f;
        printf("Color threshold increased to %.2f\n", COLOR_THRESHOLD);
        glutPostRedisplay();
    }
    else if (key == '-') {
        COLOR_THRESHOLD = std::max(0.02f, COLOR_THRESHOLD - 0.02f);
        printf("Color threshold decreased to %.2f\n", COLOR_THRESHOLD);
        glutPostRedisplay();
    }
    else if (key == 'l' || key == 'L') {
        // Increase/decrease light radius for soft shadows
        if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) {
            LIGHT_RADIUS = std::max(0.5f, LIGHT_RADIUS - 0.5f);
            printf("Light radius decreased to %.1f\n", LIGHT_RADIUS);
        } else {
            LIGHT_RADIUS += 0.5f;
            printf("Light radius increased to %.1f\n", LIGHT_RADIUS);
        }
        glutPostRedisplay();
    }
    // Depth of field and roughness controls removed as these features are disabled
    else if (key == 'n' || key == 'N') {
        // Increase/decrease number of samples
        if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) {
            SAMPLES_PER_PIXEL = std::max(4, SAMPLES_PER_PIXEL / 2);
            printf("Samples per pixel decreased to %d\n", SAMPLES_PER_PIXEL);
        } else {
            SAMPLES_PER_PIXEL = std::min(64, SAMPLES_PER_PIXEL * 2);
            printf("Samples per pixel increased to %d\n", SAMPLES_PER_PIXEL);
        }
        glutPostRedisplay();
    }
}

void initialize() {
    glMatrixMode(GL_PROJECTION);
    gluOrtho2D(XMIN, XMAX, YMIN, YMAX);
    glClearColor(0, 0, 0, 1);

    // Bounding Box

    // Checkered Floor
    // This MUST be first, as we check for intersection with it first (and then set the checkered pattern).
    Plane *floor = new Plane(glm::vec3(-45, -22.5, 0),
                             glm::vec3(45, -22.5, 0),
                             glm::vec3(45, -22.5, -160),
                             glm::vec3(-45, -22.5, -160));
    floor->setColor(glm::vec3(1, 1, 1));
    floor->setSpecularity(false);
    sceneObjects.push_back(floor);

    Plane *leftWall = new Plane(glm::vec3(-45, -22.5, 0),
                               glm::vec3(-45, -22.5, -160),
                               glm::vec3(-45, 45, -160),
                               glm::vec3(-45, 45, 0));
    leftWall->setColor(glm::vec3(1, 0, 0));
    leftWall->setSpecularity(false);
    sceneObjects.push_back(leftWall);

    Plane *rightWall = new Plane(glm::vec3(45, -22.5, 0),
                                glm::vec3(45, 45, 0),
                                glm::vec3(45, 45, -160),
                                glm::vec3(45, -22.5, -160));
    rightWall->setColor(glm::vec3(0.2, 0.8, 0.));
    rightWall->setSpecularity(false);
    sceneObjects.push_back(rightWall);

    Plane *backWall = new Plane(glm::vec3(-45, -22.5, -160),
                               glm::vec3(45, -22.5, -160),
                               glm::vec3(45, 45, -160),
                               glm::vec3(-45, 45, -160));
    backWall->setColor(glm::vec3(0.9, 0.9, 0.9));
    backWall->setSpecularity(true);
    backWall->setShininess(100.0);
    backWall->setReflectivity(true, 0.9);
    sceneObjects.push_back(backWall);

    Plane *frontWall = new Plane(glm::vec3(-45, -22.5, 0),
                                glm::vec3(-45, 45, 0),
                                glm::vec3(45, 45, 0),
                                glm::vec3(45, -22.5, 0));
    frontWall->setColor(glm::vec3(0.68, 0.85, 0.9));
    sceneObjects.push_back(frontWall);

    Plane *ceiling = new Plane(glm::vec3(-45, 45, 0),
                              glm::vec3(-45, 45, -160),
                              glm::vec3(45, 45, -160),
                              glm::vec3(45, 45, 0));
    ceiling->setColor(glm::vec3(0.7, 0.7, 0.7));
    ceiling->setSpecularity(false);
    sceneObjects.push_back(ceiling);

    // Big Sphere on the Left
    Sphere *sphereBig = new Sphere(glm::vec3(0,0,-140), 6.0);
    sphereBig->setColor(glm::vec3(1, 0.5, 1));
    sphereBig->setSpecularity(true);
    sphereBig->setReflectivity(true, 0.1);
    sceneObjects.push_back(sphereBig);


    Cylinder* cylinder = new Cylinder(glm::vec3(0, -22.5, -140), 6.0, 17.0, glm::vec3(0.0, 0.8, 0.5), true);
    cylinder->setColor(glm::vec3(1,1,1));
    cylinder->setTexture(&cylinderTexture);
    sceneObjects.push_back(cylinder);

    // Transparent Sphere
    Sphere *sphere1 = new Sphere(glm::vec3(-18, -9, -75), 6.0);
    sphere1->setColor(glm::vec3(0.7, 0.7, 1.0));
    sphere1->setReflectivity(true, 0.1);
    sphere1->setTransparency(true, 0.7);
    sceneObjects.push_back(sphere1);

    Cylinder* cylinder3 = new Cylinder(glm::vec3(-18,-22.5,-75), 6.0, 8.0, glm::vec3(0.0, 0.8, 0.5), true );
    cylinder3->setColor(glm::vec3(0.6,1,0.6));
    sceneObjects.push_back(cylinder3);

    // Refractive Sphere
    Sphere *sphereRefract = new Sphere(glm::vec3(18, -9, -75), 6.0);
    sphereRefract->setColor(glm::vec3(0.2, 1, 0.2));
    sphereRefract->setRefractivity(true, 1, 1.02);
    sphereRefract->setReflectivity(true, 0.2);
    sceneObjects.push_back(sphereRefract);

    Cylinder* cylinder2 = new Cylinder(glm::vec3(18,-22.5,-75), 6.0, 8.0, glm::vec3(0.0, 0.8, 0.5), true );
    cylinder2->setColor(glm::vec3(0.4,0,0.4));
    sceneObjects.push_back(cylinder2);

    // Blue sphere (in the middle of the torus)
    Sphere *sphere3 = new Sphere(glm::vec3(25, 0, -100), 4.0);
    sphere3->setColor(glm::vec3(0, 1, 1));
    sphere3->setReflectivity(true, 0.8);
    sceneObjects.push_back(sphere3);

    Sphere *sphere4 = new Sphere(glm::vec3(-25, 0, -100), 4.0);
    sphere4->setColor(glm::vec3(0, 1, 1));
    sphere4->setReflectivity(true, 0.8);
    sceneObjects.push_back(sphere4);

    // Torus
    Torus* torus = new Torus(glm::vec3(25, 0, -100), 9.0f, 3.0f);
    torus->rotate(-55.0f, glm::vec3(0, 1, 0));
    torus->setColor(glm::vec3(0.8, 0.2, 0.2));
    torus->setSpecularity(true);
    torus->setShininess(50.0f);
    torus->setReflectivity(true, 0.3f);
    sceneObjects.push_back(torus);

    Torus* torus2 = new Torus(glm::vec3(-25, 0, -100), 9.0f, 3.0f);
    torus2->rotate(55.0f, glm::vec3(0, 1, 0));
    torus2->setColor(glm::vec3(0.2, 0.8, 0.2));
    torus2->setSpecularity(true);
    torus2->setShininess(50.0f);
    torus2->setReflectivity(true, 0.3f);
    sceneObjects.push_back(torus2);

    // Flat pink sphere
    Sphere* flat = new Sphere(glm::vec3(25,-75,-100), 12.0f);
    flat->scale(glm::vec3(1.0f, 0.2f, 1.0f));
    flat->setColor(glm::vec3(0.8,0.2,0.6));
    sceneObjects.push_back(flat);

    Sphere* flat2 = new Sphere(glm::vec3(-25,-75,-100), 12.0f);
    flat2->scale(glm::vec3(1.0f, 0.2f, 1.0f));
    flat2->setColor(glm::vec3(0.2,0.8,0.4));
    sceneObjects.push_back(flat2);
}
int main(int argc, char *argv[]) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB );
	glutInitWindowSize(1000, 1000);
	glutInitWindowPosition(20, 20);
	glutCreateWindow("Raytracing");

	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	initialize();

	glutMainLoop();
	return 0;
}
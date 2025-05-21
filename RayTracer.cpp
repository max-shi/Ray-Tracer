/*==================================================================================
* COSC 363  Computer Graphics
* Department of Computer Science and Software Engineering, University of Canterbury.
*
* A basic ray tracer
* See Lab07.pdf   for details.
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
using namespace std;

TextureBMP texture;
const float EDIST = 40;
const int NUMDIV = 500;
const int MAX_STEPS = 2;
const float XMIN = -20.0;
const float XMAX = 20.0;
const float YMIN = -20.0;
const float YMAX = 20.0;
bool antiAliasingEnabled = true;
const int MAX_ADAPTIVE_DEPTH = 1;
float colorThreshold = 0.1f;

vector<SceneObject*> sceneObjects;

//---The most important function in a ray tracer! ----------------------------------
glm::vec3 trace(Ray ray, int step) {
    glm::vec3 backgroundCol(0);
    glm::vec3 lightPos(0, 40, -90);
    glm::vec3 ambientLight(0.2, 0.2, 0.2);
    glm::vec3 color(0);
    SceneObject* obj;

    ray.closestPt(sceneObjects);
    if(ray.index == -1) return backgroundCol;
    obj = sceneObjects[ray.index];

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

    glm::vec3 LightVec = lightPos - ray.hit;
    float lightDist = glm::length(LightVec);
    float attenuation = 1.0f / (1.0f + 0.0008f * lightDist + 0.0004f * lightDist * lightDist);
    float lightIntensity = 1.2f;
    color = obj->lighting(lightPos, -ray.dir, ray.hit) * attenuation * lightIntensity;
    color += ambientLight * obj->getColor();

    Ray shadowRay(ray.hit, LightVec);
    shadowRay.closestPt(sceneObjects);
    if (shadowRay.index > -1 && shadowRay.dist < lightDist) {
        color = ambientLight * obj->getColor();
    }

    if (obj->isTransparent() && step < MAX_STEPS) {
        float tran = obj->getTransparencyCoeff();
        
        // Create a new ray that continues through the object
        // The Ray constructor already has ray stepping built in (RSTEP = 0.005f)
        Ray transparentRay(ray.hit, ray.dir);
        
        // Find the next object this ray hits
        transparentRay.closestPt(sceneObjects);
        
        // If it hits something, blend the colors based on transparency coefficient
        if (transparentRay.index > -1 && transparentRay.index != ray.index) {
            glm::vec3 transparentColor = trace(transparentRay, step + 1);
            color = (1.0f - tran) * color + tran * transparentColor;
        }
    }

    if (obj->isReflective() && step < MAX_STEPS) {
        float rho = obj->getReflectionCoeff();
        glm::vec3 normalVec = obj->normal(ray.hit);
        glm::vec3 reflectedDir = glm::reflect(ray.dir, normalVec);
        Ray reflectedRay(ray.hit, reflectedDir);
        glm::vec3 reflectedColor = trace(reflectedRay, step + 1);
        color = color + (rho * reflectedColor);
    }

    color = glm::clamp(color, 0.0f, 1.0f);
    return color;
}

// Helper function for adaptive sampling
bool needsSubdivision(const glm::vec3& col1, const glm::vec3& col2) {
    float diff = glm::length(col1 - col2);
    return diff > colorThreshold;
}

// Recursive adaptive sampling function
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

    // Check if subdivision is needed
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
        return (samples[0] + samples[1] + samples[2] + samples[3]) / 4.0f;
    } else {
        bool similar = true;
        for (int i = 1; i < 4; i++) {
            if (needsSubdivision(samples[i], samples[0])) {
                similar = false;
                break;
            }
        }

        if (similar) {
            return (samples[0] + samples[1] + samples[2] + samples[3]) / 4.0f;
        } else if (currentDepth < maxDepth-1) {
            return adaptiveSample(x, y, width, height, eye, maxDepth, currentDepth+1, samples[0]);
        } else {
            return (samples[0] + samples[1] + samples[2] + samples[3]) / 4.0f;
        }
    }
}

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
    else if (key == '+') {
        colorThreshold += 0.02f;
        printf("Color threshold increased to %.2f\n", colorThreshold);
        glutPostRedisplay();
    }
    else if (key == '-') {
        colorThreshold = std::max(0.02f, colorThreshold - 0.02f);
        printf("Color threshold decreased to %.2f\n", colorThreshold);
        glutPostRedisplay();
    }
}

void initialize() {
    glMatrixMode(GL_PROJECTION);
    gluOrtho2D(XMIN, XMAX, YMIN, YMAX);
    glClearColor(0, 0, 0, 1);

    // Bounding Box

    // Checkered Floor
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
    Sphere *sphereBig = new Sphere(glm::vec3(-12,0,-85), 10.0);
    sphereBig->setColor(glm::vec3(1, 1, 1));
    sphereBig->setSpecularity(true);
    sceneObjects.push_back(sphereBig);


    // Transparent Sphere
    Sphere *sphere1 = new Sphere(glm::vec3(15, 0, -50), 3.0);
    sphere1->setColor(glm::vec3(0.7, 0.7, 1.0));
    sphere1->setReflectivity(true, 0.1);
    // sphere1->setRefractivity(true);
    sphere1->setTransparency(true, 0.9);
    sceneObjects.push_back(sphere1);

    // Blue sphere (in the middle of the torus)
    Sphere *sphere3 = new Sphere(glm::vec3(25, 0, -100), 4.0);
    sphere3->setColor(glm::vec3(0, 1, 1));
    sphere3->setReflectivity(true, 0.8);
    sceneObjects.push_back(sphere3);

    // Torus
    Torus* torus = new Torus(glm::vec3(25, 0, -100), 9.0f, 3.0f);
    torus->rotate(-55.0f, glm::vec3(0, 1, 0));
    torus->setColor(glm::vec3(0.8, 0.2, 0.2));
    torus->setSpecularity(true);
    torus->setShininess(50.0f);
    torus->setReflectivity(true, 0.3f);
    sceneObjects.push_back(torus);

    // Flat pink sphere
    Sphere* flat = new Sphere(glm::vec3(25,-75,-100), 12.0f);
    flat->scale(glm::vec3(1.0f, 0.2f, 1.0f));
    flat->setColor(glm::vec3(0.8,0.2,0.6));
    sceneObjects.push_back(flat);
}

int main(int argc, char *argv[]) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB );
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(20, 20);
	glutCreateWindow("Raytracing");

	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	initialize();

	glutMainLoop();
	return 0;
}
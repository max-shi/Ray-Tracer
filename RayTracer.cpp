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
const int NUMDIV = 400;
const int MAX_STEPS = 2;
const float XMIN = -20.0;  // Widened viewing window
const float XMAX = 20.0;   // Widened viewing window
const float YMIN = -20.0;  // Widened viewing window
const float YMAX = 20.0;   // Widened viewing window

vector<SceneObject*> sceneObjects;

//---The most important function in a ray tracer! ---------------------------------- 
//   Computes the colour value obtained by tracing a ray and finding its 
//     closest point of intersection with objects in the scene.
//----------------------------------------------------------------------------------
glm::vec3 trace(Ray ray, int step) {
	glm::vec3 backgroundCol(0);					//Background colour = (0,0,0)
	glm::vec3 lightPos(0, 40, -90);				//Light's position at the top of the Cornell box
	// Add a small ambient light component to ensure all surfaces have some illumination
	glm::vec3 ambientLight(0.2, 0.2, 0.2);	// Slightly increased ambient light
	glm::vec3 color(0);
	SceneObject* obj;

	ray.closestPt(sceneObjects);				//Compare the ray with all objects in the scene
	if(ray.index == -1) return backgroundCol;		//no intersection
	obj = sceneObjects[ray.index];
	
	// Check if the hit object is the floor (index 0) to create checkered pattern
	if (ray.index == 0)
	{
		// Create a checkered pattern for the floor
		int squareSize = 10; // Size of each square in the checkered pattern
		
		// Calculate checkerboard pattern based on x and z coordinates
		int ix = floor(ray.hit.x / squareSize);
		int iz = floor(ray.hit.z / squareSize);
		
		// If (ix+iz) is even, color is white, otherwise black
		if ((ix + iz) % 2 == 0) {
			color = glm::vec3(1, 1, 1); // White square
		} else {
			color = glm::vec3(0, 0, 0); // Black square
		}
		obj->setColor(color);
	}
	//object on which the closest point of intersection is found

	// Calculate distance to light for attenuation
	glm::vec3 LightVec = lightPos - ray.hit;
	float lightDist = glm::length(LightVec);
	
	// Apply light attenuation based on distance (slightly reduced attenuation)
	float attenuation = 1.0f / (1.0f + 0.0008f * lightDist + 0.0004f * lightDist * lightDist);
	
	// Adjust light intensity
	float lightIntensity = 1.2f;
	
	// Calculate base color with lighting and apply attenuation
	color = obj->lighting(lightPos, -ray.dir, ray.hit) * attenuation * lightIntensity;
	
	// Add ambient light to ensure all surfaces have some illumination
	color += ambientLight * obj->getColor();
	
	// Shadow calculation
	Ray shadowRay(ray.hit, LightVec);
	shadowRay.closestPt(sceneObjects);
	if (shadowRay.index > -1 && shadowRay.dist < lightDist) {
		// In shadow - reduce to ambient only
		color = ambientLight * obj->getColor();
	}
	
	// Handle transparency - must be done before reflections
	if (obj->isTransparent() && step < MAX_STEPS)
	{
		float tran = obj->getTransparencyCoeff();
		// Create a ray that continues in the same direction from the hit point
		Ray transparentRay(ray.hit, ray.dir);
		// Find the next object behind the transparent object
		transparentRay.closestPt(sceneObjects);
		
		if (transparentRay.index > -1) {
			// Get the color from the object behind the transparent object
			glm::vec3 transparentColor = trace(transparentRay, step + 1);
			// Blend the current color with the transparent color
			color = (1.0f - tran) * color + tran * transparentColor;
		}
	}
	
	// Handle reflections after transparency
	if (obj->isReflective() && step < MAX_STEPS)
	{
		float rho = obj->getReflectionCoeff();
		glm::vec3 normalVec = obj->normal(ray.hit);
		glm::vec3 reflectedDir = glm::reflect(ray.dir, normalVec);
		Ray reflectedRay(ray.hit, reflectedDir);
		glm::vec3 reflectedColor = trace(reflectedRay, step + 1);
		color = color + (rho * reflectedColor);
	}
	
	// Ensure color values don't exceed 1.0 (moved after reflection)
	color = glm::clamp(color, 0.0f, 1.0f);
	
	return color;
}

//---The main display module -----------------------------------------------------------
// In a ray tracing application, it just displays the ray traced image by drawing
// each cell as a quad.
//---------------------------------------------------------------------------------------
void display() {
	float xp, yp;  //grid point
	float cellX = (XMAX - XMIN) / NUMDIV;  //cell width
	float cellY = (YMAX - YMIN) / NUMDIV;  //cell height
	glm::vec3 eye(0., 0., 0.);

	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glBegin(GL_QUADS);  //Each cell is a tiny quad.

	for (int i = 0; i < NUMDIV; i++) {	//Scan every cell of the image plane
		xp = XMIN + i * cellX;
		for (int j = 0; j < NUMDIV; j++) {
			yp = YMIN + j * cellY;

			glm::vec3 dir(xp + 0.5 * cellX, yp + 0.5 * cellY, -EDIST);	//direction of the primary ray

			Ray ray = Ray(eye, dir);

			glm::vec3 col = trace(ray, 1); //Trace the primary ray and get the colour value
			glColor3f(col.r, col.g, col.b);
			glVertex2f(xp, yp);				//Draw each cell with its color value
			glVertex2f(xp + cellX, yp);
			glVertex2f(xp + cellX, yp + cellY);
			glVertex2f(xp, yp + cellY);
		}
	}

	glEnd();
	glFlush();
}

//---This function initializes the scene ------------------------------------------- 
//   Specifically, it creates scene objects (spheres, planes, cones, cylinders etc)
//     and add them to the list of scene objects.
//   It also initializes the OpenGL 2D orthographc projection matrix for drawing the
//     the ray traced image.
//----------------------------------------------------------------------------------
void initialize() {
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(XMIN, XMAX, YMIN, YMAX);
	glClearColor(0, 0, 0, 1);

	// Create a Cornell-style box with camera positioned closer to the back wall
	// The walls are positioned to make it seem like the camera is closer to the back

	// Floor - checkered pattern will be implemented in the trace function
	Plane *floor = new Plane(glm::vec3(-45, -22.5, 0),      // Point A
		                    glm::vec3(45, -22.5, 0),       // Point B
		                    glm::vec3(45, -22.5, -160),    // Point C
		                    glm::vec3(-45, -22.5, -160));  // Point D
	floor->setColor(glm::vec3(1, 1, 1));                // White base color for checkered pattern
	floor->setSpecularity(false);
	sceneObjects.push_back(floor);

	// Left wall - RED (changed from green)
	Plane *leftWall = new Plane(glm::vec3(-45, -22.5, 0),    // Point A
		                       glm::vec3(-45, -22.5, -160), // Point B
		                       glm::vec3(-45, 45, -160),    // Point C
		                       glm::vec3(-45, 45, 0));      // Point D
	leftWall->setColor(glm::vec3(1, 0, 0));              // Red
	leftWall->setSpecularity(false);
	sceneObjects.push_back(leftWall);

	// Right wall - Silver, Reflective.
	Plane *rightWall = new Plane(glm::vec3(45, -22.5, 0),    // Point A
		                        glm::vec3(45, 45, 0),       // Point B
		                        glm::vec3(45, 45, -160),    // Point C
		                        glm::vec3(45, -22.5, -160));// Point D
	rightWall->setColor(glm::vec3(0.9, 0.9, 0.9));             // Green
	rightWall->setSpecularity(true);                       // Enable specularity for mirror-like shine
	rightWall->setShininess(100.0);                        // High shininess for mirror effect
	rightWall->setReflectivity(true, 0.9);
	sceneObjects.push_back(rightWall);

	// Back wall (behind camera) - Green
	Plane *backWall = new Plane(glm::vec3(-45, -22.5, -160),   // Point A
		                       glm::vec3(45, -22.5, -160),    // Point B
		                       glm::vec3(45, 45, -160),       // Point C
		                       glm::vec3(-45, 45, -160));     // Point D
	backWall->setColor(glm::vec3(0.2039f, 0.9216f, 0.5725f));              // White
	backWall->setSpecularity(false);
	sceneObjects.push_back(backWall);

	// Front wall (where you are facing) - BABY BLUE (swapped with back wall)
	Plane *frontWall = new Plane(glm::vec3(-45, -22.5, 0),   // Point A
		                        glm::vec3(-45, 45, 0),      // Point B
		                        glm::vec3(45, 45, 0),       // Point C
		                        glm::vec3(45, -22.5, 0));   // Point D
	frontWall->setColor(glm::vec3(0.68, 0.85, 0.9));     // Light baby blue
	frontWall->setSpecularity(false);
	sceneObjects.push_back(frontWall);

	// Ceiling - GREY with light source
	// Reordered vertices to make normal point upward
	Plane *ceiling = new Plane(glm::vec3(-45, 45, 0),        // Point A
		                      glm::vec3(-45, 45, -160),     // Point B
		                      glm::vec3(45, 45, -160),      // Point C
		                      glm::vec3(45, 45, 0));        // Point D
	ceiling->setColor(glm::vec3(0.7, 0.7, 0.7));         // Grey
	ceiling->setSpecularity(false);
	sceneObjects.push_back(ceiling);

	// Add four orbs spaced evenly in the scene
	// Orb 1 - top left (transparent sphere)
	Sphere *sphere1 = new Sphere(glm::vec3(-10, 0, -40), 3.0);
	sphere1->setColor(glm::vec3(0.7, 0.7, 1.0));   // Light blue tint for transparency
	sphere1->setReflectivity(true, 0.3);    // Reduced reflection coefficient for transparency
	sphere1->setTransparency(true, 0.7);    // Make this sphere transparent with 70% transparency
	sceneObjects.push_back(sphere1);


	// // Orb 3 - bottom left
	// Sphere *sphere3 = new Sphere(glm::vec3(-5, 0, -80), 13.0);
	// sphere3->setColor(glm::vec3(0, 1, 1));   // Cyan
	// sphere3->setReflectivity(true, 0.8);   // Reduced reflection coefficient
	// sceneObjects.push_back(sphere3);

	// Add a torus to the scene
	Torus* torus = new Torus(glm::vec3(15, 0, -40), 9.0f, 3.0f);
	// torus->translate(glm::vec3(-5, 0, -80)); // Move to position
	torus->rotate(15.0f, glm::vec3(1, 1, 0)); // Rotate around Y
	// torus->rotate(15.0f, glm::vec3(1, 0, 0)); // Additional rotation around X
	torus->setColor(glm::vec3(0.8, 0.2, 0.2));  // Reddish colo
	torus->setSpecularity(true);
	torus->setShininess(50.0f);
	torus->setReflectivity(true, 0.3f);
	sceneObjects.push_back(torus);
}

int main(int argc, char *argv[]) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB );
	glutInitWindowSize(1000, 1000);
	glutInitWindowPosition(20, 20);
	glutCreateWindow("Raytracing");

	glutDisplayFunc(display);
	initialize();

	glutMainLoop();
	return 0;
}

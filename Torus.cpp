#include "Torus.h"
#include <cmath>
#include <complex>
#include <algorithm>

using namespace std;

/**
 * Constructor for the Torus class
 *
 * @param center The center point of the torus in 3D space
 * @param majorRadius The distance from the center of the tube to the center of the torus (outer radius)
 * @param minorRadius The radius of the tube itself (inner radius)
 * 
 * Initializes a torus with the specified center and radii. The transform and invTransform
 * matrices are initialized to identity matrices, representing no transformation.
 */
Torus::Torus(glm::vec3 center, float majorRadius, float minorRadius)
    : center(center), majorRadius(majorRadius), minorRadius(minorRadius) {
    transform = glm::mat4(1.0f);
    invTransform = glm::mat4(1.0f);
}

/**
 * Rotate method for the torus around a specified axis
 * 
 * @param angle The rotation angle in degrees
 * @param axis The axis of rotation as a 3D vector
 * 
 * This method rotates the torus around its center by the specified angle along the given axis.
 * The rotation is performed by:
 * 1. Translating the torus to the origin
 * 2. Applying the rotation
 * 3. Translating back to the original position
 * 
 * The transformation matrix is updated accordingly, and the inverse transformation matrix
 * is recalculated to ensure correct ray transformations during intersection tests.
 */
void Torus::rotate(float angle, glm::vec3 axis) {
    // Create a translation matrix to move to origin
    glm::mat4 toOrigin = glm::translate(glm::mat4(1.0f), -center);
    
    // Create a rotation matrix
    glm::mat4 rotation = glm::rotate(glm::mat4(1.0f), glm::radians(angle), axis);
    
    // Create a translation matrix to move back
    glm::mat4 fromOrigin = glm::translate(glm::mat4(1.0f), center);
    
    // Apply the transformations in sequence: translate to origin, rotate, translate back
    transform = fromOrigin * rotation * toOrigin * transform;
    invTransform = glm::inverse(transform);
}

/**
 * Translates the torus by a specified vector
 * 
 * @param translation The translation vector in 3D space
 * 
 * This method moves the torus by the specified translation vector.
 * The transformation matrix is updated by applying a translation matrix,
 * and the inverse transformation matrix is recalculated to ensure correct
 * ray transformations during intersection tests.
 */
void Torus::translate(glm::vec3 translation) {
    glm::mat4 translationMat = glm::translate(glm::mat4(1.0f), translation);
    transform = translationMat * transform;
    invTransform = glm::inverse(transform);
}

/**
 * Scales the torus by a specified factor along each axis
 * 
 * @param scale The scale factors as a 3D vector (x, y, z)
 * 
 * This method resizes the torus by the specified scale factors along each axis.
 * The transformation matrix is updated by applying a scaling matrix,
 * and the inverse transformation matrix is recalculated to ensure correct
 * ray transformations during intersection tests.
 */
void Torus::scale(glm::vec3 scale) {
    glm::mat4 scaleMat = glm::scale(glm::mat4(1.0f), scale);
    transform = scaleMat * transform;
    invTransform = glm::inverse(transform);
}
/**
 * Tests if a ray intersects with the torus's bounding sphere
 * 
 * @param p0 The origin point of the ray in world space
 * @param dir The direction vector of the ray in world space
 * @param t Output parameter that stores the intersection distance if an intersection occurs
 * @return true if the ray intersects the bounding sphere, false otherwise
 * 
 * This method is an optimization that quickly determines if a ray might intersect the torus
 * by first checking for intersection with a bounding sphere. The bounding sphere has a radius
 * equal to the sum of the major and minor radii of the torus.
 * 
 * The ray is first transformed to the torus's local coordinate system using the inverse
 * transformation matrix. Then a standard ray-sphere intersection test is performed.
 * If an intersection is found, the parameter t is set to the distance to the intersection point.
 */
bool Torus::intersectBoundingSphere(glm::vec3 p0, glm::vec3 dir, float& t) const {
    // Transform ray to torus's local coordinate system
    glm::vec4 localP0_4 = invTransform * glm::vec4(p0, 1.0f);
    glm::vec4 localDir_4 = invTransform * glm::vec4(dir, 0.0f);

    glm::vec3 localP0 = glm::vec3(localP0_4) - center;
    glm::vec3 localDir = glm::vec3(localDir_4);

    // The bounding sphere radius is majorRadius + minorRadius
    float boundRadius = majorRadius + minorRadius;

    // Solve quadratic equation for sphere intersection
    float a = glm::dot(localDir, localDir);
    float b = 2.0f * glm::dot(localP0, localDir);
    float c = glm::dot(localP0, localP0) - boundRadius * boundRadius;

    float discriminant = b * b - 4.0f * a * c;

    if (discriminant < 0.0f) {
        // if discriminant is negative -> then we false
        return false;
    }

    float sqrtDisc = sqrt(discriminant);
    float t1 = (-b - sqrtDisc) / (2.0f * a);
    float t2 = (-b + sqrtDisc) / (2.0f * a);

    t = (t1 > 0.0f) ? t1 : t2;
    return (t > 0.0f);
}

/**
 * Solves a quartic equation using the Durand-Kerner method
 * 
 * @param coeffs A vector of 5 coefficients [a, b, c, d, e] for the quartic equation ax^4 + bx^3 + cx^2 + dx + e = 0
 * @return A vector of real roots of the quartic equation
 * 
 * This method implements the Durand-Kerner numerical algorithm to find all roots of a quartic polynomial.
 * The algorithm works by starting with initial guesses for the roots and iteratively refining them until convergence.
 * 
 * The steps are:
 * 1. Start with initial guesses for the roots (equally spaced around the unit circle in the complex plane)
 * 2. Iteratively refine each root using the Durand-Kerner formula
 * 3. Check for convergence after each iteration
 * 4. Extract and return only the real roots (roots with negligible imaginary parts)
 * 
 * This method is used internally by the intersect method to solve the quartic equation that arises
 * when finding the intersection of a ray with a torus.
 *
 * Sources: https://en.wikipedia.org/wiki/Durand%E2%80%93Kerner_method
 * https://stackoverflow.com/questions/40274198/ray-tracing-quartic-surfaces <- This one is GOOD
 */
vector<float> Torus::solveQuartic(const vector<double>& coeffs) const {
    const int maxIterations = 100;
    const double tolerance = 1e-6;

    // Initial guesses - equally spaced around the unit circle
    // TODO : are these ideal?
    complex<double> roots[4] = {
        complex<double>(0.4, 0.9),
        complex<double>(-0.9, 0.4),
        complex<double>(-0.4, -0.9),
        complex<double>(0.9, -0.4)
    };

    // Durand-Kerner iteration
    for (int iter = 0; iter < maxIterations; ++iter) {
        complex<double> newRoots[4];
        // start with conversion assumed
        bool converged = true;

        // If in here, any time the convergence is false, then we break (at the bottom).
        for (int i = 0; i < 4; ++i) {
            complex<double> numerator =
                coeffs[0] * roots[i] * roots[i] * roots[i] * roots[i] +
                coeffs[1] * roots[i] * roots[i] * roots[i] +
                coeffs[2] * roots[i] * roots[i] +
                coeffs[3] * roots[i] +
                coeffs[4];

            complex<double> denominator = 1.0;
            for (int j = 0; j < 4; ++j) {
                if (i != j) {
                    denominator *= (roots[i] - roots[j]);
                }
            }

            newRoots[i] = roots[i] - (numerator / denominator);

            if (abs(newRoots[i] - roots[i]) > tolerance) {
                // if the val between the root and new root is greater than the tolerance, we set converged to false (and we run the alg. again)
                converged = false;
            }
        }

        if (converged) break;

        for (int i = 0; i < 4; ++i) {
            roots[i] = newRoots[i];
        }
    }

    // Extract real roots
    vector<float> realRoots;
    for (int i = 0; i < 4; ++i) {
        if (abs(roots[i].imag()) < tolerance) {
            realRoots.push_back(static_cast<float>(roots[i].real()));
        }
    }

    return realRoots;
}

/**
 * Calculates the intersection of a ray with the torus
 * 
 * @param p0 The origin point of the ray in world space
 * @param dir The direction vector of the ray in world space
 * @return The distance to the intersection point, or -1 if no intersection occurs
 * 
 * This method implements the ray-torus intersection algorithm. The steps are:
 * 1. First check for intersection with a bounding sphere as an optimization
 * 2. Transform the ray to the torus's local coordinate system
 * 3. Compute the coefficients for the quartic equation that represents the intersection
 * 4. Solve the quartic equation using the solveQuartic method
 * 5. Find the smallest positive root, which represents the closest intersection point
 * 
 * The mathematical derivation of the quartic equation comes from the standard equation of a torus
 * and the parametric equation of a ray. When combined and simplified, they yield a 4th degree polynomial
 * whose roots represent the intersection points.
 * 
 * This method overrides the intersect method from the SceneObject base class.
 */
float Torus::intersect(glm::vec3 p0, glm::vec3 dir) {
    // First check bounding sphere
    float boundT;
    if (!intersectBoundingSphere(p0, dir, boundT)) {
        return -1.0f; // Early exit if no intersection with bounding sphere
    }

    // Transform ray to torus's local coordinate system using the transformation matrices
    glm::vec4 localP0_4 = invTransform * glm::vec4(p0, 1.0f);
    glm::vec4 localDir_4 = invTransform * glm::vec4(dir, 0.0f);
    
    glm::vec3 localP0 = glm::vec3(localP0_4) - center;
    glm::vec3 localDir = glm::vec3(localDir_4);

    // Precompute values
    double Ex = localDir.x, Ey = localDir.y, Ez = localDir.z;
    double Dx = localP0.x, Dy = localP0.y, Dz = localP0.z;
    double A = majorRadius, B = minorRadius;

    double Ex2 = Ex * Ex;
    double Ey2 = Ey * Ey;
    double Ez2 = Ez * Ez;
    double Dx2 = Dx * Dx;
    double Dy2 = Dy * Dy;
    double Dz2 = Dz * Dz;

    // Compute coefficients for quartic equation
    double G = 4 * A * A * (Ex2 + Ey2);
    double H = 8 * A * A * (Dx * Ex + Dy * Ey);
    double I = 4 * A * A * (Dx2 + Dy2);
    double J = Ex2 + Ey2 + Ez2;
    double K = 2 * (Dx * Ex + Dy * Ey + Dz * Ez);
    double L = Dx2 + Dy2 + Dz2 + (A * A - B * B);

    // Quartic coefficients: J^2*u^4 + 2*J*K*u^3 + (2*J*L + K^2 - G)*u^2 + (2*K*L - H)*u + (L^2 - I)
    vector<double> coeffs = {
        J * J,
        2 * J * K,
        2 * J * L + K * K - G,
        2 * K * L - H,
        L * L - I
    };

    // Solve quartic equation
    vector<float> roots = solveQuartic(coeffs);

    // Find smallest positive real root
    float t = -1;
    for (float root : roots) {
        if (root > 0 && (t < 0 || root < t)) {
            t = root;
        }
    }

    return t;
}

/**
 * Calculates the surface normal at a point on the torus
 * 
 * @param p The point on the torus surface in world space
 * @return The normalized surface normal vector at the point
 * 
 * This method computes the surface normal at a given point on the torus. The steps are:
 * 1. Transform the point to the torus's local coordinate system
 * 2. Find the center of the tube (Q) that is closest to the given point
 *    - This is done by projecting the point onto the xy-plane and scaling to the major radius
 * 3. The normal vector is the direction from Q to the point (P - Q)
 * 4. Transform the normal back to world space using the transpose of the inverse transform
 * 5. Normalize the resulting vector
 * 
 * The normal vector is essential for lighting calculations in the ray tracer.
 * 
 * This method overrides the normal method from the SceneObject base class.
 */
glm::vec3 Torus::normal(glm::vec3 p) {
    // Transform point to torus's local coordinate system using the transformation matrices
    glm::vec4 localP_4 = invTransform * glm::vec4(p, 1.0f);
    glm::vec3 localP = glm::vec3(localP_4) - center;

    // Compute Q (center of the tube nearest to P)
    double Px = localP.x, Py = localP.y, Pz = localP.z;
    double alpha = majorRadius / sqrt(Px * Px + Pz * Pz);
    glm::vec3 Q(alpha * Px, 0.0f, alpha * Pz);

    // Compute normal vector N = P - Q
    glm::vec3 N = localP - Q;

    // Transform normal back to world space and return normalized
    // For normals, we need to use the transpose of the inverse of the upper 3x3 part of the transform matrix
    glm::vec4 worldN_4 = glm::transpose(invTransform) * glm::vec4(N, 0.0f);
    return glm::normalize(glm::vec3(worldN_4));
}


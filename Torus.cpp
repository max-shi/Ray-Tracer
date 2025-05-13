#include "Torus.h"
#include <cmath>
#include <complex>
#include <algorithm>

using namespace std;


/*Considerations to make still: TODO
 *  The Durand-Kerner method is relatively expensive compared to analytical solutions for simpler shapes.
 The number of iterations (currently 100) can be adjusted based on quality needs.
 *
 */
#include "Torus.h"
#include <cmath>
#include <complex>
#include <algorithm>

Torus::Torus(glm::vec3 center, float majorRadius, float minorRadius)
    : center(center), majorRadius(majorRadius), minorRadius(minorRadius) {
    transform = glm::mat4(1.0f); // Identity matrix
    invTransform = glm::mat4(1.0f);
}

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

void Torus::translate(glm::vec3 translation) {
    glm::mat4 translationMat = glm::translate(glm::mat4(1.0f), translation);
    transform = translationMat * transform;
    invTransform = glm::inverse(transform);
}

void Torus::scale(glm::vec3 scale) {
    glm::mat4 scaleMat = glm::scale(glm::mat4(1.0f), scale);
    transform = scaleMat * transform;
    invTransform = glm::inverse(transform);
}
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
        return false; // No intersection
    }

    float sqrtDisc = sqrt(discriminant);
    float t1 = (-b - sqrtDisc) / (2.0f * a);
    float t2 = (-b + sqrtDisc) / (2.0f * a);

    t = (t1 > 0.0f) ? t1 : t2;
    return (t > 0.0f);
}

vector<float> Torus::solveQuartic(const vector<double>& coeffs) const {
    const int maxIterations = 100;
    const double tolerance = 1e-6;

    // Initial guesses - equally spaced around the unit circle
    complex<double> roots[4] = {
        complex<double>(0.4, 0.9),
        complex<double>(-0.9, 0.4),
        complex<double>(-0.4, -0.9),
        complex<double>(0.9, -0.4)
    };

    // Durand-Kerner iteration
    for (int iter = 0; iter < maxIterations; ++iter) {
        complex<double> newRoots[4];
        bool converged = true;

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


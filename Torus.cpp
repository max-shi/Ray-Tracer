#include "Torus.h"
#include <cmath>
#include <complex>
#include <algorithm>

using namespace std;


/*Considerations to make still: TODO
 *  The Durand-Kerner method is relatively expensive compared to analytical solutions for simpler shapes.
 You might want to add a bounding sphere check first to avoid the quartic solver when possible.
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
    glm::mat4 rotation = glm::rotate(glm::mat4(1.0f), glm::radians(angle), axis);
    transform = rotation * transform;
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
    // Transform ray to torus's local coordinate system
    glm::vec3 localP0 = p0 - center;
    glm::vec3 localDir = dir;

    // For simplicity, assume torus is axis-aligned (axis = (0,1,0))
    // In a complete implementation, we would transform the ray to object space

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
    // Transform point to torus's local coordinate system
    glm::vec3 localP = p - center;

    // For simplicity, assume torus is axis-aligned (axis = (0,1,0))
    // In a complete implementation, we would transform the point to object space

    // Compute Q (center of the tube nearest to P)
    double Px = localP.x, Py = localP.y, Pz = localP.z;
    double alpha = majorRadius / sqrt(Px * Px + Pz * Pz);
    glm::vec3 Q(alpha * Px, 0.0f, alpha * Pz);

    // Compute normal vector N = P - Q
    glm::vec3 N = localP - Q;

    // Return normalized normal
    return glm::normalize(N);
}
// C++ include
#include <iostream>
#include <string>
#include <vector>

// Utilities for the Assignment
#include "utils.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

void raytrace_sphere()
{
    std::cout << "Simple ray tracer, one sphere with orthographic projection" << std::endl;

    const std::string filename("sphere_orthographic.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask
 
    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is orthographic, pointing in the direction -z and covering the
    // unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = camera_view_direction;

            // Intersect with the sphere
            // NOTE: this is a special case of a sphere centered in the origin and for orthographic rays aligned with the z axis
            Vector2d ray_on_xy(ray_origin(0), ray_origin(1));
            const double sphere_radius = 0.9;

            if (ray_on_xy.norm() < sphere_radius)
            {
                // The ray hit the sphere, compute the exact intersection point
                Vector3d ray_intersection(
                    ray_on_xy(0), ray_on_xy(1),
                    sqrt(sphere_radius * sphere_radius - ray_on_xy.squaredNorm()));

                // Compute normal at the intersection point
                Vector3d ray_normal = ray_intersection.normalized();

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

bool intersect_ray_parallelogram(const Vector3d& ray_origin, const Vector3d& ray_direction, 
                                 const Vector3d& pgram_origin, const Vector3d& pgram_u, const Vector3d& pgram_v, 
                                 Vector3d& intersection)
{
    Matrix3d A;
    A << -pgram_u, -pgram_v, ray_direction;
    Vector3d b = pgram_origin - ray_origin;
    
    Vector3d x = A.inverse() * b;

    double u = x(0);
    double v = x(1);
    double t = x(2);

    if ((u >= 0 && u <= 1) && (v >= 0 && v <= 1) && t >= 0) 
    {
        intersection = ray_origin + t * ray_direction;
        return true;
    }
    return false;
}

void raytrace_parallelogram()
{
    std::cout << "Simple ray tracer, one parallelogram with orthographic projection" << std::endl;

    const std::string filename("plane_orthographic.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // Parameters of the parallelogram (position of the lower-left corner + two sides)
    const Vector3d pgram_origin(-0.5, -0.5, 0);
    const Vector3d pgram_u(0, 0.7, -10);
    const Vector3d pgram_v(1, 0.4, 0);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = camera_view_direction;

            // TODO: Check if the ray intersects with the parallelogram
            // TODO: The ray hit the parallelogram, compute the exact intersection point
            Vector3d ray_intersection(0, 0, 0);
            if (intersect_ray_parallelogram(ray_origin, ray_direction, pgram_origin, pgram_u, pgram_v, ray_intersection))
            {
            
                // TODO: Compute normal at the intersection point
                Vector3d ray_normal = -(pgram_u.cross(pgram_v).normalized());

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_perspective()
{
    std::cout << "Simple ray tracer, one parallelogram with perspective projection" << std::endl;

    const std::string filename("plane_perspective.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // TODO: Parameters of the parallelogram (position of the lower-left corner + two sides)
    const Vector3d pgram_origin(-0.5, -0.5, 0);
    const Vector3d pgram_u(0, 0.7, -10);
    const Vector3d pgram_v(1, 0.4, 0);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // TODO: Prepare the ray (origin point and direction)
            const Vector3d ray_origin = camera_origin;
            const Vector3d ray_direction = pixel_center - camera_origin;

            // TODO: Check if the ray intersects with the parallelogram
            // TODO: The ray hit the parallelogram, compute the exact intersection point
            Vector3d ray_intersection(0, 0, 0);
            if (intersect_ray_parallelogram(ray_origin, ray_direction, pgram_origin, pgram_u, pgram_v, ray_intersection))
            {
            
                // TODO: Compute normal at the intersection point
                Vector3d ray_normal = -(pgram_u.cross(pgram_v).normalized());

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_shading()
{
    std::cout << "Simple ray tracer, one sphere with different shading" << std::endl;

    const std::string filename("shading.png");
    MatrixXd R = MatrixXd::Zero(800, 800); // Store the red color
    MatrixXd G = MatrixXd::Zero(800, 800); // Store the green color
    MatrixXd B = MatrixXd::Zero(800, 800); // Store the blue color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / A.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / A.rows(), 0);

    // Sphere setup
    const Vector3d sphere_center(0, 0, 0);
    const double sphere_radius = 0.9;

    // Material params
    const Vector3d diffuse_color(1, 0, 1);
    const double specular_exponent = 100;
    const Vector3d specular_color(0., 0, 1);

    // Single light source
    const Vector3d light_position(-1, 1, 1);
    const double ambient = 0.1;

    for (unsigned i = 0; i < A.cols(); ++i)
    {
        for (unsigned j = 0; j < A.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            const Vector3d ray_origin = camera_origin;
            const Vector3d ray_direction = (pixel_center - camera_origin).normalized();

            // Intersect with the sphere
            // TODO: implement the generic ray sphere intersection
            Vector3d oc = ray_origin - sphere_center;
            double a = ray_direction.dot(ray_direction);
            double b = 2.0 * oc.dot(ray_direction);
            double c = oc.dot(oc) - sphere_radius * sphere_radius;

            double discriminant = b * b - 4 * a * c;

            if (discriminant > 0)
            {
                double t = (-b - sqrt(discriminant)) / (2.0 * a);
                if (t > 0)
                {
                    // TODO: The ray hit the sphere, compute the exact intersection point
                    Vector3d ray_intersection = ray_origin + t * ray_direction;

                    // TODO: Compute normal at the intersection point
                    Vector3d ray_normal = (ray_intersection - sphere_center).normalized();

                    // TODO: Add shading parameter here
                    Vector3d l = (light_position - ray_intersection).normalized();
                    Vector3d v = (camera_origin - ray_intersection).normalized();
                    Vector3d h = ((v + l) / ((v + l).norm())).normalized();

                    Vector3d diffuse = diffuse_color * fmax(0, (ray_normal.dot(l)));
                    Vector3d specular = specular_color * pow(fmax(0, (ray_normal.dot(h))), specular_exponent);

                    // Simple diffuse model
                    R(i, j) = ambient + diffuse(0) + specular(0);
                    G(i, j) = ambient + diffuse(1) + specular(1);
                    B(i, j) = ambient + diffuse(2) + specular(2);

                    // Clamp to zero
                    R(i, j) = std::max(R(i, j), 0.);
                    G(i, j) = std::max(G(i, j), 0.);
                    B(i, j) = std::max(B(i, j), 0.);

                    // Disable the alpha mask for this pixel
                    A(i, j) = 1;
            
        
                }
            }
        }
    }

    // Save to png
    write_matrix_to_png(R, G, B, A, filename);
}


int main()
{
    raytrace_sphere();
    raytrace_parallelogram();
    raytrace_perspective();
    raytrace_shading();

    return 0;
}

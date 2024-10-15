#pragma once

#include <Eigen/Core>

class VertexAttributes
{
public:
    VertexAttributes(double x = 0, double y = 0, double z = 0, double w = 1)
    {
        position << x, y, z, w;
        color << 1, 1, 1, 1;
    }

    // Interpolates the vertex attributes
    static VertexAttributes interpolate(
        const VertexAttributes &a,
        const VertexAttributes &b,
        const VertexAttributes &c,
        const double alpha,
        const double beta,
        const double gamma)
    {
        VertexAttributes r;
        r.position = alpha * a.position + beta * b.position + gamma * c.position;
        r.color = a.color * alpha + b.color * beta + c.color * gamma;
        return r;
    }

    Eigen::Vector4d position;
    Eigen::Vector4f color;
    Eigen::Vector3f normal;
};

class FragmentAttributes
{
public:
    FragmentAttributes(double r = 0, double g = 0, double b = 0, double a = 1)
    {
        color << r, g, b, a;
    }

    Eigen::Vector4d color;
    Eigen::Vector4f position;
};

class FrameBufferAttributes
{
public:
    FrameBufferAttributes(double r = 0, double g = 0, double b = 0, double a = 1)
    {
        color << r, g, b, a;
        depth = 0;
    }

    Eigen::Matrix<double, 4, 1> color;
    float depth;
    
};

class UniformAttributes
{
public:
    Eigen::Matrix4d view; // Add this line for the view transformation matrix
    Eigen::Vector4f color;
};
// C++ include
#include <iostream>
#include <string>
#include <vector>

// Utilities for the Assignment
#include "raster.h"

#include <gif.h>
#include <fstream>

#include <Eigen/Geometry>
// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

using namespace std;
using namespace Eigen;

//Image height
const int H = 480;

//Camera settings
const double near_plane = 1.5;       //AKA focal length
const double far_plane = near_plane * 100;
const double field_of_view = 0.7854; //45 degrees
const double aspect_ratio = 1.5;
const bool is_perspective = true;
const Vector3d camera_position(0, 0, 3);
const Vector3d camera_gaze(0, 0, -1);
const Vector3d camera_top(0, 1, 0);

//Object
const std::string data_dir = DATA_DIR;
const std::string mesh_filename(data_dir + "bunny.off");
MatrixXd vertices; // n x 3 matrix (n points)
MatrixXi facets;   // m x 3 matrix (m triangles)

//Material for the object
const Vector3d obj_diffuse_color(0.5, 0.5, 0.5);
const Vector3d obj_specular_color(0.2, 0.2, 0.2);
const double obj_specular_exponent = 256.0;

//Lights
std::vector<Vector3d> light_positions;
std::vector<Vector3d> light_colors;
//Ambient light
const Vector3d ambient_light(0.3, 0.3, 0.3);

//Fills the different arrays
void setup_scene()
{
    //Loads file
    std::ifstream in(mesh_filename);
    if (!in.good())
    {
        std::cerr << "Invalid file " << mesh_filename << std::endl;
        exit(1);
    }
    std::string token;
    in >> token;
    int nv, nf, ne;
    in >> nv >> nf >> ne;
    vertices.resize(nv, 3);
    facets.resize(nf, 3);
    for (int i = 0; i < nv; ++i)
    {
        in >> vertices(i, 0) >> vertices(i, 1) >> vertices(i, 2);
    }
    for (int i = 0; i < nf; ++i)
    {
        int s;
        in >> s >> facets(i, 0) >> facets(i, 1) >> facets(i, 2);
        assert(s == 3);
    }

    //Lights
    light_positions.emplace_back(8, 8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(6, -8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(4, 8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(2, -8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(0, 8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(-2, -8, 0);
    light_colors.emplace_back(16, 16, 16);

    light_positions.emplace_back(-4, 8, 0);
    light_colors.emplace_back(16, 16, 16);
}

void build_uniform(UniformAttributes &uniform)
{
    //TODO: setup uniform

    //TODO: setup camera, compute w, u, v
    Vector3d w = -camera_gaze.normalized();
    Vector3d u = camera_top.cross(w).normalized();
    Vector3d v = w.cross(u);

    //TODO: compute the camera transformation
    Matrix4d cam_transform;
    cam_transform << u(0), u(1), u(2), -u.dot(camera_position),
                     v(0), v(1), v(2), -v.dot(camera_position),
                     w(0), w(1), w(2), -w.dot(camera_position),
                     0, 0, 0, 1;

    //TODO: setup projection matrix
    double t = near_plane * tan(field_of_view / 2);
    double r = t * aspect_ratio;
    const float l = -r;
    const float b = -t;
    const float n = -near_plane;
    const float f = -far_plane;


    Matrix4d ortho;
    ortho << 2 / (r - l), 0, 0, -(r + l) / (r - l),
        0, 2 / (t - b), 0, -(t + b) / (t - b),
        0, 0, 2 / (n - f), -(n + f) / (n - f),
        0, 0, 0, 1;

    Matrix4d persp_to_ortho;
    persp_to_ortho << n, 0, 0, 0,
            0, n, 0, 0,
            0, 0, (n + f), (-f * n),
            0, 0, 1, 0;


    Matrix4d P;
    if (!is_perspective)
    {
        //TODO setup prespective camera
         P = ortho * persp_to_ortho;
    }
    else
    {
        P = ortho;
    }
    uniform.view = P * cam_transform;
}

void simple_render(Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
{
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;

    program.VertexShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        //TODO: fill the shader
        VertexAttributes out = va;
        out.position = uniform.view * va.position;
        out.position /= out.position(3); // Homogeneous divide to normalize coordinates
        return out;
    };

    program.FragmentShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        //TODO: fill the shader
        return FragmentAttributes(1, 0, 0);
    };

    program.BlendingShader = [](const FragmentAttributes &fa, const FrameBufferAttributes &previous) {
        //TODO: fill the shader
        return FrameBufferAttributes(fa.color[0], fa.color[1], fa.color[2], fa.color[3]);
    };

    std::vector<VertexAttributes> vertex_attributes;
    //TODO: build the vertex attributes from vertices and facets
    for (int i = 0; i < facets.rows(); ++i)
    {
        vertex_attributes.push_back(VertexAttributes(vertices(facets(i, 0), 0), vertices(facets(i, 0), 1), vertices(facets(i, 0), 2)));
        vertex_attributes.push_back(VertexAttributes(vertices(facets(i, 1), 0), vertices(facets(i, 1), 1), vertices(facets(i, 1), 2)));
        vertex_attributes.push_back(VertexAttributes(vertices(facets(i, 2), 0), vertices(facets(i, 2), 1), vertices(facets(i, 2), 2)));

    }

    rasterize_triangles(program, uniform, vertex_attributes, frameBuffer);
}

Matrix4d compute_rotation(const double alpha)
{
    //TODO: Compute the rotation matrix of angle alpha on the y axis around the object barycenter
    // Compute the barycenter of the object
    Matrix4d res;
    double c = cos(alpha);
    double s = sin(alpha);

    res << c, 0, s, 0,
        0, 1, 0, 0,
        -s, 0, c, 0,
        0, 0, 0, 1;

    Matrix4d I;
    I.setIdentity();
    res = res * I;
    return res;
}

void wireframe_render(const double alpha, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
{
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;

    Matrix4d trafo = compute_rotation(alpha);

    program.VertexShader = [&](const VertexAttributes &va, const UniformAttributes &uniform) {
        //TODO: fill the shader
        VertexAttributes out = va;
        out.position = uniform.view * trafo * va.position;
        out.position /= out.position(3); // Homogeneous divide to normalize coordinates
        return out;
    };

    program.FragmentShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        //TODO: fill the shader
        return FragmentAttributes(1, 0, 0);
    };

    program.BlendingShader = [](const FragmentAttributes &fa, const FrameBufferAttributes &previous) {
        //TODO: fill the shader
        return FrameBufferAttributes(fa.color[0], fa.color[1], fa.color[2], fa.color[3]);
    };

    std::vector<VertexAttributes> vertex_attributes;

    //TODO: generate the vertex attributes for the edges and rasterize the lines
    //TODO: use the transformation matrix
    for (int i = 0; i < facets.rows(); ++i)
    {
        // Extract the indices of the vertices of the triangle
        int v0 = facets(i, 0);
        int v1 = facets(i, 1);
        int v2 = facets(i, 2);

        // Create VertexAttributes for each edge of the triangle
        vertex_attributes.push_back(VertexAttributes(vertices(v0, 0), vertices(v0, 1), vertices(v0, 2)));
        vertex_attributes.push_back(VertexAttributes(vertices(v1, 0), vertices(v1, 1), vertices(v1, 2)));

        vertex_attributes.push_back(VertexAttributes(vertices(v1, 0), vertices(v1, 1), vertices(v1, 2)));
        vertex_attributes.push_back(VertexAttributes(vertices(v2, 0), vertices(v2, 1), vertices(v2, 2)));

        vertex_attributes.push_back(VertexAttributes(vertices(v2, 0), vertices(v2, 1), vertices(v2, 2)));
        vertex_attributes.push_back(VertexAttributes(vertices(v0, 0), vertices(v0, 1), vertices(v0, 2)));
    }

    rasterize_lines(program, uniform, vertex_attributes, 0.5, frameBuffer);
}

void get_shading_program(Program &program)
{
    program.VertexShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        //TODO: transform the position and the normal
        //TODO: compute the correct lighting
        VertexAttributes out;
        out.position = uniform.view * va.position;
        out.normal = va.normal;

        Vector3d accumulated_light_color = Vector3d::Zero();
        
        for (size_t i = 0; i < light_positions.size(); ++i)
        {
            const Vector3d light_position = light_positions[i];
            const Vector3d light_color = light_colors[i];

            Vector3d vertex_position = out.position.head<3>();
            Vector3d normal_vector = out.normal.cast<double>();

            Vector3d light_direction = (light_position - vertex_position).normalized();
            Vector3d diffuse = obj_diffuse_color * std::max(light_direction.dot(normal_vector), 0.0);

            Vector3d view_direction = (camera_position - vertex_position).normalized();
            Vector3d half_vector = (light_direction + view_direction).normalized();
            Vector3d specular = obj_specular_color * std::pow(std::max(normal_vector.dot(half_vector), 0.0), obj_specular_exponent);

            Vector3d light_distance = light_position - vertex_position;
            accumulated_light_color += (diffuse + specular).cwiseProduct(light_color) / light_distance.squaredNorm();
        }

        accumulated_light_color += ambient_light;
        out.color << accumulated_light_color[0], accumulated_light_color[1], accumulated_light_color[2], 1.0;
        
        return out;
    };
    

    program.FragmentShader = [](const VertexAttributes &va, const UniformAttributes &uniform) {
        //TODO: create the correct fragment
        return FragmentAttributes(va.color[0], va.color[1], va.color[2]);
    };

    program.BlendingShader = [](const FragmentAttributes &fa, const FrameBufferAttributes &previous) {
        //TODO: implement the depth check
         if (fa.position(2) < previous.depth)
        {
            FrameBufferAttributes out(fa.color[0], fa.color[1], fa.color[2], fa.color[3]);
            out.depth = fa.position(2);
            return out;
        }
        else
        {
            return previous;
        }

    };
}

void flat_shading(const double alpha, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
{
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;
    get_shading_program(program);
    Eigen::Matrix4d trafo = compute_rotation(alpha);

    std::vector<VertexAttributes> vertex_attributes;
    //TODO: compute the normals
    //TODO: set material colors
    for (int i = 0; i < facets.rows(); ++i)
    {
        Vector3d vertex_a, vertex_b, vertex_c;

        Vector3i vertex_indices = facets.row(i).transpose();

        vertex_a = vertices.row(vertex_indices(0)).transpose();
        vertex_b = vertices.row(vertex_indices(1)).transpose();
        vertex_c = vertices.row(vertex_indices(2)).transpose();

        // Apply rotation transformation
        vertex_a = (trafo * Vector4d(vertex_a[0], vertex_a[1], vertex_a[2], 1.0)).head<3>();
        vertex_b = (trafo * Vector4d(vertex_b[0], vertex_b[1], vertex_b[2], 1.0)).head<3>();
        vertex_c = (trafo * Vector4d(vertex_c[0], vertex_c[1], vertex_c[2], 1.0)).head<3>();

        VertexAttributes attribute_a(vertex_a[0], vertex_a[1], vertex_a[2]);
        VertexAttributes attribute_b(vertex_b[0], vertex_b[1], vertex_b[2]);
        VertexAttributes attribute_c(vertex_c[0], vertex_c[1], vertex_c[2]);

        Vector3d vector_u = vertex_b - vertex_a;
        Vector3d vector_v = vertex_c - vertex_a;

        Vector3f normal = vector_v.cross(vector_u).normalized().cast<float>();

        attribute_a.normal = normal;
        attribute_b.normal = normal;
        attribute_c.normal = normal;

        vertex_attributes.push_back(attribute_a);
        vertex_attributes.push_back(attribute_b);
        vertex_attributes.push_back(attribute_c);
    }

    rasterize_triangles(program, uniform, vertex_attributes, frameBuffer);
}

void pv_shading(const double alpha, Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> &frameBuffer)
{
    UniformAttributes uniform;
    build_uniform(uniform);
    Program program;
    get_shading_program(program);

    Eigen::Matrix4d trafo = compute_rotation(alpha);

    //TODO: compute the vertex normals as vertex normal average

    std::vector<VertexAttributes> vertex_attributes;
    //TODO: create vertex attributes
    //TODO: set material colors
    std::vector<Vector3d> vertex_normals(vertices.rows(), Vector3d(0, 0, 0));
    std::vector<int> normal_count(vertices.rows(), 0);

    // Calculate vertex normals by averaging face normals
    for (int i = 0; i < facets.rows(); ++i)
    {
        Vector3d vertex_a, vertex_b, vertex_c;

        Vector3i vertex_indices = facets.row(i).transpose();

        vertex_a = vertices.row(vertex_indices(0)).transpose();
        vertex_b = vertices.row(vertex_indices(1)).transpose();
        vertex_c = vertices.row(vertex_indices(2)).transpose();

        // Apply rotation transformation
        vertex_a = (trafo * Vector4d(vertex_a[0], vertex_a[1], vertex_a[2], 1.0)).head<3>();
        vertex_b = (trafo * Vector4d(vertex_b[0], vertex_b[1], vertex_b[2], 1.0)).head<3>();
        vertex_c = (trafo * Vector4d(vertex_c[0], vertex_c[1], vertex_c[2], 1.0)).head<3>();


        Vector3d vector_u = vertex_b - vertex_a;
        Vector3d vector_v = vertex_c - vertex_a;

        Vector3d normal = vector_v.cross(vector_u).normalized();

        for (int j = 0; j < 3; ++j)
        {
            vertex_normals[vertex_indices(j)] += normal;
            normal_count[vertex_indices(j)]++;
        }
    }

    // Normalize the accumulated vertex normals
    for (int i = 0; i < vertex_normals.size(); ++i)
    {
        vertex_normals[i] /= normal_count[i];
        vertex_normals[i].normalize();

        // Apply rotation to the normals
        vertex_normals[i] = (trafo.block<3, 3>(0, 0) * vertex_normals[i]).normalized();
    }

    // Create vertex attributes with vertex normals
    for (int i = 0; i < facets.rows(); ++i)
    {
        Vector3d vertex_a, vertex_b, vertex_c;

        Vector3i vertex_indices = facets.row(i).transpose();

        vertex_a = vertices.row(vertex_indices(0)).transpose();
        vertex_b = vertices.row(vertex_indices(1)).transpose();
        vertex_c = vertices.row(vertex_indices(2)).transpose();

        // Apply rotation transformation
        vertex_a = (trafo * Vector4d(vertex_a[0], vertex_a[1], vertex_a[2], 1.0)).head<3>();
        vertex_b = (trafo * Vector4d(vertex_b[0], vertex_b[1], vertex_b[2], 1.0)).head<3>();
        vertex_c = (trafo * Vector4d(vertex_c[0], vertex_c[1], vertex_c[2], 1.0)).head<3>();

        VertexAttributes attribute_a(vertex_a[0], vertex_a[1], vertex_a[2]);
        VertexAttributes attribute_b(vertex_b[0], vertex_b[1], vertex_b[2]);
        VertexAttributes attribute_c(vertex_c[0], vertex_c[1], vertex_c[2]);

        attribute_a.normal = vertex_normals[vertex_indices(0)].cast<float>();
        attribute_b.normal = vertex_normals[vertex_indices(1)].cast<float>();
        attribute_c.normal = vertex_normals[vertex_indices(2)].cast<float>();

        vertex_attributes.push_back(attribute_a);
        vertex_attributes.push_back(attribute_b);
        vertex_attributes.push_back(attribute_c);
    }

    rasterize_triangles(program, uniform, vertex_attributes, frameBuffer);
}

int main(int argc, char *argv[])
{
    setup_scene();

    int W = H * aspect_ratio;
    Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> frameBuffer(W, H);
    vector<uint8_t> image;

    frameBuffer.setConstant(FrameBufferAttributes()); 
    simple_render(frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("simple.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);

    frameBuffer.setConstant(FrameBufferAttributes()); 
    wireframe_render(0, frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("wireframe.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);

    frameBuffer.setConstant(FrameBufferAttributes()); 
    flat_shading(0, frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("flat_shading.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);

    frameBuffer.setConstant(FrameBufferAttributes()); 
    pv_shading(0, frameBuffer);
    framebuffer_to_uint8(frameBuffer, image);
    stbi_write_png("pv_shading.png", frameBuffer.rows(), frameBuffer.cols(), 4, image.data(), frameBuffer.rows() * 4);

    //TODO: add the animation

    
    int delay = 25;
    GifWriter g;
    int num_frames = 24;
    double step_size = EIGEN_PI / 10;

    // Wireframe render
    frameBuffer.setConstant(FrameBufferAttributes()); 
    GifBegin(&g, "wireframe.gif", frameBuffer.rows(), frameBuffer.cols(), delay);
    for (float i = 1; i < num_frames; i += step_size)
    {
        frameBuffer.setConstant(FrameBufferAttributes());
        wireframe_render(i, frameBuffer);
        framebuffer_to_uint8(frameBuffer, image);
        GifWriteFrame(&g, image.data(), frameBuffer.rows(), frameBuffer.cols(), delay);
    }
    GifEnd(&g);

    //flat_shading render
    frameBuffer.setConstant(FrameBufferAttributes()); 
    GifBegin(&g, "flat_shading.gif", frameBuffer.rows(), frameBuffer.cols(), delay);
    for (float i = 1; i < num_frames; i += step_size)
    {
        frameBuffer.setConstant(FrameBufferAttributes());
        flat_shading(i, frameBuffer);
        framebuffer_to_uint8(frameBuffer, image);
        GifWriteFrame(&g, image.data(), frameBuffer.rows(), frameBuffer.cols(), delay);
    }
    GifEnd(&g);

    //pv_shading render
    frameBuffer.setConstant(FrameBufferAttributes()); 
    GifBegin(&g, "pv_shading.gif", frameBuffer.rows(), frameBuffer.cols(), delay);
    for (float i = 1; i < num_frames; i += step_size)
    {
        frameBuffer.setConstant(FrameBufferAttributes());
        pv_shading(i, frameBuffer);
        framebuffer_to_uint8(frameBuffer, image);
        GifWriteFrame(&g, image.data(), frameBuffer.rows(), frameBuffer.cols(), delay);
    }
    GifEnd(&g);

    return 0;
}

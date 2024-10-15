////////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <complex>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

#include <Eigen/Dense>
// Shortcut to avoid  everywhere, DO NOT USE IN .h
using namespace Eigen;
////////////////////////////////////////////////////////////////////////////////

const std::string root_path = DATA_DIR;

// Computes the determinant of the matrix whose columns are the vector u and v
double inline det(const Vector2d &u, const Vector2d &v)
{
    // TODO
    return u.x() * v.y() - u.y() * v.x();
}

// Return true iff [a,b] intersects [c,d]
bool intersect_segment(const Vector2d &a, const Vector2d &b, const Vector2d &c, const Vector2d &d)
{
    // TODO
    Vector2d ab = b - a;
    Vector2d ac = c - a;
    Vector2d ad = d - a;
    Vector2d cd = d - c;
    Vector2d ca = a - c;
    Vector2d cb = b - c;

    double det1 = det(ab, ac);
    double det2 = det(ab, ad);
    double det3 = det(cd, ca);
    double det4 = det(cd, cb);

    if ((det1 * det2 < 0) && (det3 * det4 < 0))
        return true;
        
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool is_inside(const std::vector<Vector2d> &poly, const Vector2d &query)
{
    // 1. Compute bounding box and set coordinate of a point outside the polygon
    // TODO
    double min_x = poly[0].x();
    double max_x = poly[0].x();
    double min_y = poly[0].y();
    double max_y = poly[0].y();

    for (const auto& vertex : poly)
    {
        if (vertex.x() < min_x) {
            min_x = vertex.x();
        }
        if (vertex.x() > max_x) {
            max_x = vertex.x();
        }
        if (vertex.y() < min_y) {
            min_y = vertex.y();
        }
        if (vertex.y() > max_y) {
            max_y = vertex.y();
        }
    }

    Vector2d outside(max_x + 1, query.y());
    // 2. Cast a ray from the query point to the 'outside' point, count number of intersections
    // TODO
    
    int count = 0;
    for (size_t i = 0; i < poly.size(); ++i)
    {
        Vector2d a = poly[i];
        Vector2d b = poly[(i + 1) % poly.size()];
        if (intersect_segment(query, outside, a, b))
        {
            count++;
        }
    }

    return (count % 2) == 1;
   
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Vector2d> load_xyz(const std::string &filename)
{
    std::vector<Vector2d> points;
    std::ifstream in(filename);
    // TODO

    if (!in){
        std::cerr << "Error opening file: " << filename << std::endl;
        return points;
    }

    int num_points;
    in >> num_points;
    
    double x, y, z;
    while (in >> x >> y >> z)
    {
        points.push_back(Vector2d(x, y));
    }

    in.close();

    return points;
    
}

void save_xyz(const std::string &filename, const std::vector<Vector2d> &points)
{
    // TODO
    std::ofstream out(filename);

    if (!out)
    {
        std::cerr << "Error creating file: " << filename << std::endl;
        return;
    }
    
    out << points.size() << std::endl;   
    for (const auto &point : points)
    {
        out << point.x() << " " << point.y() << " 0" << std::endl;
    }
    
    out.close();
}

std::vector<Vector2d> load_obj(const std::string &filename)
{
    std::ifstream in(filename);
    std::vector<Vector2d> points;
    std::vector<Vector2d> poly;
    char key;
    while (in >> key)
    {
        if (key == 'v')
        {
            double x, y, z;
            in >> x >> y >> z;
            points.push_back(Vector2d(x, y));
        }
        else if (key == 'f')
        {
            std::string line;
            std::getline(in, line);
            std::istringstream ss(line);
            int id;
            while (ss >> id)
            {
                poly.push_back(points[id - 1]);
            }
        }
    }
    return poly;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    const std::string points_path = root_path + "/points.xyz";
    const std::string poly_path = root_path + "/polygon.obj";

    std::vector<Vector2d> points = load_xyz(points_path);

    ////////////////////////////////////////////////////////////////////////////////
    //Point in polygon
    std::vector<Vector2d> poly = load_obj(poly_path);
    std::vector<Vector2d> result;
    for (size_t i = 0; i < points.size(); ++i)
    {
        if (is_inside(poly, points[i]))
        {
            result.push_back(points[i]);
        }
    }
    save_xyz("output.xyz", result);


    return 0;
}

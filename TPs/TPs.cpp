#include <cmath>
#include <cstdlib>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

typedef unsigned long long size_t; // should already be defined in cstdib,
// but I redefine it in case size_t is not unsigned long long (even it's very unlikely),
// to ensure compatibility with std::stoull (see Mesh::readOFF)

double sqr(double a) { return a * a; }

std::vector<std::string> split(const std::string& line, const char& delim = ' ') {
    std::stringstream stream(line);
    std::string string;
    std::vector<std::string> result;
    while (std::getline(stream, string, delim)) {
        result.push_back(string);
    }
    return result;
}

class Vector3 {
public:
    explicit Vector3(double x = 0, double y = 0, double z = 0) {
        coord[0] = x;
        coord[1] = y;
        coord[2] = z;
    }
    double& operator[](int i) { return coord[i]; }
    double operator[](int i) const { return coord[i]; }

    Vector3& operator+=(const Vector3& v) {
        coord[0] += v[0];
        coord[1] += v[1];
        coord[2] += v[2];
        return *this;
    }

    double norm2() const {
        return sqr(coord[0]) + sqr(coord[1]) + sqr(coord[2]);
    }

    void normalize() {
        double norm{ sqrt(norm2()) };
        coord[0] /= norm;
        coord[1] /= norm;
        coord[2] /= norm;
    };

    double coord[3];
};
Vector3 operator+(const Vector3& a, const Vector3& b) {
    return Vector3(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector3 operator-(const Vector3& a, const Vector3& b) {
    return Vector3(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector3 operator*(const Vector3& a, double b) {
    return Vector3(a[0] * b, a[1] * b, a[2] * b);
}
 Vector3 operator*(double a, const Vector3& b) {
    return Vector3(a * b[0], a * b[1], a * b[2]);
}
Vector3 operator*(const Vector3& a, const Vector3& b) {
    return Vector3(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}
Vector3 operator/(const Vector3& a, double b) {
    return Vector3(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector3& a, const Vector3& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector3 cross(const Vector3& a, const Vector3& b) {
    return Vector3(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

class Vertex3 {
public:

    Vertex3(std::initializer_list<double> position, size_t index_adjacent_triangle, bool is_virtual = false) :
        index_adjacent_triangle(index_adjacent_triangle), is_virtual(is_virtual) { updatePosition(position); }

    void updatePosition(std::initializer_list<double> position) {

        if (position.size() != 3) throw std::invalid_argument("Invalid initializer list : 3 values are expected");
        for (size_t i{ 0 }; i < 3; i++) this->position[i] = position.begin()[i];
    }

    Vector3 position;
    size_t index_adjacent_triangle;
    bool is_virtual;
};

class Triangle {
public:

    Triangle(std::initializer_list<size_t> vertices_indices, std::initializer_list<size_t> neighbouring_triangles_indices) {
        updateVertices(vertices_indices);
        updateNT(neighbouring_triangles_indices);
    }

    void updateVertices(std::initializer_list<size_t> vertices_indices) {

        if (vertices_indices.size() != 3) throw std::invalid_argument("Invalid initializer list : 3 values are expected");
        for (size_t i{ 0 }; i < 3; i++) this->vertices_indices[i] = vertices_indices.begin()[i];
    }

    void updateNT(std::initializer_list<size_t> neighbouring_triangles_indices) {

        if (neighbouring_triangles_indices.size() != 3) throw std::invalid_argument("Invalid initializer list : 3 values are expected");
        for (size_t i{ 0 }; i < 3; i++) this->neighbouring_triangles_indices[i] = neighbouring_triangles_indices.begin()[i];
    }

    size_t vertices_indices[3];
    size_t neighbouring_triangles_indices[3];
    bool is_virtual{ false };
};

class Mesh {
public:

    void addVertex(const Vertex3& Vertex3) {
        vertices.push_back(Vertex3);
    }

    void addTriangle(const Triangle& triangle) {
        triangles.push_back(triangle);
    }

    void writeOFF(const char* filename) {

        check_virtual();

        std::ofstream ofs;
        ofs.open(filename);
        if (ofs.bad()) {
            std::cout << "Can not write file " << filename << std::endl;
            exit(1);
        }

        ofs << "OFF\n";
        ofs << vertices.size() - nb_virtual_vertices << " " << triangles.size() - nb_virtual_triangles << " 0\n";

        for (size_t i{ 0 }; i < vertices.size(); i++) {
            if (!vertices[i].is_virtual) {
                for (size_t j{ 0 }; j < 3; j++) {
                    ofs << vertices[i].position[j];
                    if (j < 2) ofs << " ";
                }
                ofs << "\n";
            }
        }

        for (size_t i{ 0 }; i < triangles.size(); i++) {
            if (!triangles[i].is_virtual) {
                ofs << "3";
                for (size_t j{ 0 }; j < 3; j++) {
                    ofs << " " << triangles[i].vertices_indices[j];
                }
                ofs << "\n";
            }
        }

        ofs.close();
        std::cout << filename << " successfully written" << std::endl;
    }

    void readOFF(const char* filename) {

        std::ifstream ifs;
        ifs.open(filename);
        if (ifs.bad()) {
            std::cout << "Can not read file " << filename << "\n";
            exit(1);
        }

        std::map<std::set<size_t>, std::pair<size_t, size_t>> map;

        size_t nb_vertices{ 0 }, nb_triangles{ 0 }, line_counter{ 0 }, bound1{ 2 }, bound2{ 0 }, bound3{ 0 };
        std::string line;
        while (std::getline(ifs, line)) { // parsing using order (comments not supported)
                
            std::vector<std::string> v_line{ split(line) };

            if (line_counter == 1) { // 2nd line (assuming there is a first line with OFF)

                nb_vertices = std::stoull(v_line[0]);
                bound2 = bound1 + nb_vertices;
                for (size_t i{ 0 }; i < nb_vertices; i++) addVertex(Vertex3({ 0, 0, 0 }, 0));

                nb_triangles = std::stoull(v_line[1]);
                bound3 = bound2 + nb_triangles;
                for (size_t i{ 0 }; i < nb_triangles; i++) addTriangle(Triangle({ 0, 0, 0 }, { 0, 0, 0 }));
                // default neighbouring triangle will be the first one in case the mesh is not closed
            }
            else if (bound1 <= line_counter && line_counter < bound2) {

                size_t i{ line_counter - bound1 };
                vertices[i].updatePosition({ std::stod(v_line[0]), std::stod(v_line[1]), std::stod(v_line[2]) });
            }
            else if (bound2 <= line_counter && line_counter < bound3) {

                size_t i{ line_counter - bound2 };

                // we trust the file that the vertices indices are in the right (trigonometric) order
                // only support triangles reading
                triangles[i].vertices_indices[0] = std::stoull(v_line[1]);
                triangles[i].vertices_indices[1] = std::stoull(v_line[2]);
                triangles[i].vertices_indices[2] = std::stoull(v_line[3]);
                // adjacent triangle for a Vertex3 will be the last adjacent triangle "seen"
                vertices[triangles[i].vertices_indices[0]].index_adjacent_triangle = i;
                vertices[triangles[i].vertices_indices[1]].index_adjacent_triangle = i;
                vertices[triangles[i].vertices_indices[2]].index_adjacent_triangle = i;

                std::vector<std::set<size_t>> edges;
                edges.reserve(3);
                edges.emplace_back(std::set<size_t>{ triangles[i].vertices_indices[1], triangles[i].vertices_indices[2] }); // relative index 0
                edges.emplace_back(std::set<size_t>{ triangles[i].vertices_indices[2], triangles[i].vertices_indices[0] }); // relative index 1
                edges.emplace_back(std::set<size_t>{ triangles[i].vertices_indices[0], triangles[i].vertices_indices[1] }); // relative index 2
                // relative index k = index (inside of the Triangle class) of a Vertex3, the edge opposite to the Vertex3 and the triangle connected to opposite edge

                for (size_t k{ 0 }; k < 3; k++) {
                    if (map.count(edges[k])) { // 1 ...

                        size_t i_bis{ map[edges[k]].first };
                        triangles[i].neighbouring_triangles_indices[k] = i_bis;
                        triangles[i_bis].neighbouring_triangles_indices[map[edges[k]].second] = i;
                    }
                    else { // ... or 0
                        map[edges[0]] = std::pair<size_t, size_t>(i, k);
                    }
                }
            }

            line_counter++;
        }

        ifs.close();
        std::cout << filename << " successfully loaded" << std::endl;
    }

    void check_virtual() {

        nb_virtual_vertices = 0;
        for (size_t i{ 0 }; i < vertices.size(); i++)
            if (vertices[i].is_virtual)
                nb_virtual_vertices += 1;

        nb_virtual_triangles = 0;
        for (size_t i{ 0 }; i < triangles.size(); i++) {

            bool triangle_is_virtual{ false };

            for (size_t j{ 0 }; j < 3; j++)
                if (vertices[triangles[i].vertices_indices[j]].is_virtual)
                    triangle_is_virtual = true;

            triangles[i].is_virtual = triangle_is_virtual;
            if (triangles[i].is_virtual) nb_virtual_triangles += 1;
        }
    }

    size_t nb_virtual_vertices{ 0 };
    size_t nb_virtual_triangles{ 0 };

    std::vector<Vertex3> vertices;
    std::vector<Triangle> triangles;
};


int main() {

    // Tetrahedron
    Mesh tetrahedron;
    // It is not regular but ... simpler to draw in a standard(x, y, z) reference frame
    tetrahedron.addVertex(Vertex3({ 0, 0, 0 }, 0));
    tetrahedron.addVertex(Vertex3({ 1, 0, 0 }, 0));
    tetrahedron.addVertex(Vertex3({ 0, 1, 0 }, 0));
    tetrahedron.addVertex(Vertex3({ 0, 0, 1 }, 1));
    tetrahedron.addTriangle(Triangle({ 0, 1, 2 }, { 3, 2, 1 }));
    tetrahedron.addTriangle(Triangle({ 0, 3, 1 }, { 3, 0, 2 }));
    tetrahedron.addTriangle(Triangle({ 0, 2, 3 }, { 3, 1, 0 }));
    tetrahedron.addTriangle(Triangle({ 1, 3, 2 }, { 2, 0, 1 }));
    tetrahedron.writeOFF("off_files/tetrahedron_correct.off"); // we can visualize it in 3dviewer.net

    // Checking that we get the good results for the loading
    Mesh new_tetrahedron;
    new_tetrahedron.readOFF("off_files/tetrahedron_correct.off");
    new_tetrahedron.writeOFF("off_files/tetrahedron.off"); // this file must be identical to the correct one (read as text file)
    
    // Square based pyramide
    Mesh pyramide;
    pyramide.addVertex(Vertex3({ 0, 0, 0 }, 0));
    pyramide.addVertex(Vertex3({ 1, 0, 0 }, 0));
    pyramide.addVertex(Vertex3({ 1, 1, 0 }, 1));
    pyramide.addVertex(Vertex3({ 0, 1, 0 }, 0));
    pyramide.addVertex(Vertex3({ 0.5, 0.5, 1 }, 2));
    pyramide.addTriangle(Triangle({ 0, 1, 3 }, { 1, 3, 4 }));
    pyramide.addTriangle(Triangle({ 1, 2, 3 }, { 2, 0, 5 }));
    pyramide.addTriangle(Triangle({ 2, 4, 3 }, { 3, 1, 5 }));
    pyramide.addTriangle(Triangle({ 0, 3, 4 }, { 2, 4, 0 }));
    pyramide.addTriangle(Triangle({ 0, 4, 1 }, { 5, 0, 3 }));
    pyramide.addTriangle(Triangle({ 1, 4, 2 }, { 2, 1, 4 }));
    pyramide.writeOFF("off_files/pyramide_correct.off"); // we can visualize it in 3dviewer.net

    // Checking that we get the good results for the loading
    Mesh new_pyramide;
    new_pyramide.readOFF("off_files/pyramide_correct.off");
    new_pyramide.writeOFF("off_files/pyramide.off"); // this file must be identical to the correct one (read as text file)

    // 2D bounding box
    Mesh bounding_box;
    bounding_box.addVertex(Vertex3({ 0, 0, 0 }, 0));
    bounding_box.addVertex(Vertex3({ 1, 0, 0 }, 0));
    bounding_box.addVertex(Vertex3({ 1, 2, 0 }, 1));
    bounding_box.addVertex(Vertex3({ 0, 2, 0 }, 0));
    bounding_box.addVertex(Vertex3({ 0, 0, 0 }, 2, true));
    bounding_box.addTriangle(Triangle({ 0, 1, 3 }, { 1, 3, 4 }));
    bounding_box.addTriangle(Triangle({ 1, 2, 3 }, { 0, 5, 2 }));
    bounding_box.addTriangle(Triangle({ 2, 4, 3 }, { 1, 5, 3 }));
    bounding_box.addTriangle(Triangle({ 0, 3, 4 }, { 0, 2, 4 }));
    bounding_box.addTriangle(Triangle({ 0, 4, 1 }, { 0, 3, 5 }));
    bounding_box.addTriangle(Triangle({ 1, 4, 2 }, { 1, 4, 2 }));
    bounding_box.writeOFF("off_files/bounding_box_correct.off");
    // we can visualize it in 3dviewer.net (virtual virtices and triangles are not written in the OFF file)

    // Checking that we get the good results for the loading
    Mesh new_bounding_box;
    new_bounding_box.readOFF("off_files/bounding_box_correct.off");
    new_bounding_box.writeOFF("off_files/bounding_box.off"); // this file must be identical to the correct one (read as text file)

    Mesh queen; // we can visualize both in 3dviewer.net
    queen.readOFF("off_files/queen.off");
    queen.writeOFF("off_files/queen_check.off");
    
    return 0;
}

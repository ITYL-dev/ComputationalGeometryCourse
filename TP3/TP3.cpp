#include <cmath>
#include <iostream>
#include <vector>
#include <initializer_list>
#include <set>

double sqr(double a) { return a * a; }

class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        coord[0] = x;
        coord[1] = y;
        coord[2] = z;
    }
    double& operator[](int i) { return coord[i]; }
    double operator[](int i) const { return coord[i]; }

    Vector& operator+=(const Vector& v) {
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
Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const Vector& a, double b) {
    return Vector(a[0] * b, a[1] * b, a[2] * b);
}
Vector operator*(double a, const Vector& b) {
    return Vector(a * b[0], a * b[1], a * b[2]);
}
Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}
Vector operator/(const Vector& a, double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

bool predicate_orientation(Vector p, Vector q, Vector r, Vector e_axis = Vector(0, 0, 1)) { // default normal axis : Oz

    return dot(cross(q - p, r - p), e_axis) > 0; // true: positive, false: negative
}

class Vertex {
public:

    Vertex(std::initializer_list<double> position, int index_adjacent_triangle) :
        index_adjacent_triangle(index_adjacent_triangle) {
        updatePosition(position);
    }

    void updatePosition(std::initializer_list<double> position) {
        for (int i{ 0 }; i < 3; i++) this->position[i] = position.begin()[i];
    }

    Vector position;
    int index_adjacent_triangle;
};

class Triangle {
public:

    Triangle(std::initializer_list<int> vertices_indices, std::initializer_list<int> neighbouring_triangles_indices) {
        updateVertices(vertices_indices);
        updateNT(neighbouring_triangles_indices);
    }

    void updateVertices(std::initializer_list<int> vertices_indices) {
        for (int i{ 0 }; i < 3; i++) this->vertices_indices[i] = vertices_indices.begin()[i];
    }

    void updateNT(std::initializer_list<int> neighbouring_triangles_indices) {
        for (int i{ 0 }; i < 3; i++) this->neighbouring_triangles_indices[i] = neighbouring_triangles_indices.begin()[i];
    }

    bool pointInTriangle(const Vector& P, const std::vector<Vertex>& vertices) const {
        // Get triangle vertices
        Vector A = vertices[vertices_indices[0]].position;
        Vector B = vertices[vertices_indices[1]].position;
        Vector C = vertices[vertices_indices[2]].position;

        // Compute cross products
        Vector AB = B - A, AP = P - A;
        Vector BC = C - B, BP = P - B;
        Vector CA = A - C, CP = P - C;

        // Normal of the triangle
        Vector normal = cross(AB, BC);

        // Check signs of cross products
        double d1 = dot(cross(AB, AP), normal);
        double d2 = dot(cross(BC, BP), normal);
        double d3 = dot(cross(CA, CP), normal);

        return (d1 >= 0 && d2 >= 0 && d3 >= 0) || (d1 <= 0 && d2 <= 0 && d3 <= 0);
    }

    int vertices_indices[3];
    int neighbouring_triangles_indices[3];
};

class Mesh {
public:

    void addVertex(const Vertex& Vertex) {
        vertices.push_back(Vertex);
    }

    void addTriangle(const Triangle& triangle) {
        triangles.push_back(triangle);
    }

    void splitTriangle(int triangle_index, const Vertex& new_vertex) {

        addVertex(new_vertex);
        int new_vertex_index = vertices.size() - 1;

        Triangle& triangle{ triangles[triangle_index] };
        int triangle_vertices_indices[3];
        int triangle_neighb_indices[3];
        for (int k{ 0 }; k < 3; k++) {
            triangle_vertices_indices[k] = triangle.vertices_indices[k];
            triangle_neighb_indices[k] = triangle.neighbouring_triangles_indices[k];
        }
        
        triangle.vertices_indices[2] = new_vertex_index;
        triangle.neighbouring_triangles_indices[0] = triangles.size();
        triangle.neighbouring_triangles_indices[1] = triangles.size() + 1;

        addTriangle(Triangle(
            { triangle_vertices_indices[1], triangle_vertices_indices[2], new_vertex_index},
            { triangle.neighbouring_triangles_indices[1], triangle_index, triangle_neighb_indices[0]}
        ));  

        addTriangle(Triangle(
            { triangle_vertices_indices[2], triangle_vertices_indices[0], new_vertex_index },
            { triangle_index, triangle.neighbouring_triangles_indices[0], triangle_neighb_indices[1] }
        ));
    }

    void handleOutsidePoint(int new_vertex_index) {
        std::vector<std::pair<int, int>> boundaryEdges;
        std::set<int> trianglesToRemove;

        Vector newPoint = vertices[new_vertex_index].position;

        // 1. Find the visible boundary edges
        for (int i = 0; i < triangles.size(); i++) {
            Triangle& t = triangles[i];

            bool isOutside = false;
            for (int j = 0; j < 3; j++) {
                Vector A = vertices[t.vertices_indices[j]].position;
                Vector B = vertices[t.vertices_indices[(j + 1) % 3]].position;

                // Check if the edge is "visible" from the new point
                if (predicate_orientation(A, B, newPoint)) {
                    boundaryEdges.push_back({ t.vertices_indices[j], t.vertices_indices[(j + 1) % 3] });
                    isOutside = true;
                }
            }

            if (isOutside) {
                trianglesToRemove.insert(i);
            }
        }

        // 2. Remove old triangles that are no longer valid
        std::vector<Triangle> newTriangles;
        for (int i = 0; i < triangles.size(); i++) {
            if (trianglesToRemove.find(i) == trianglesToRemove.end()) {
                newTriangles.push_back(triangles[i]);
            }
        }
        triangles = std::move(newTriangles);

        // 3. Create new triangles connecting the new point to the convex hull boundary
        for (auto& edge : boundaryEdges) {
            addTriangle(Triangle({ edge.first, edge.second, new_vertex_index }, { -1, -1, -1 }));
        }
    }


    void flipEdge(int triangle_index, int edge_index) {
        Triangle& t1 = triangles[triangle_index];

        // Get the two vertices forming the selected edge
        int v0 = t1.vertices_indices[edge_index];
        int v1 = t1.vertices_indices[(edge_index + 1) % 3];
        int v2 = t1.vertices_indices[(edge_index + 2) % 3];

        // Find the adjacent triangle sharing this edge
        int adjacent_index = t1.neighbouring_triangles_indices[edge_index];
        if (adjacent_index == -1) return; // No adjacent triangle to flip with

        Triangle& t2 = triangles[adjacent_index];

        // Find the opposite vertex in the adjacent triangle
        int opposite_index = -1;
        int v3 = -1;
        for (int i = 0; i < 3; i++) {
            if (t2.vertices_indices[i] != v0 && t2.vertices_indices[i] != v1) {
                opposite_index = i;
                v3 = t2.vertices_indices[i];
                break;
            }
        }

        if (v3 == -1) return; // Something went wrong

        // Identify edges in t2
        int edge_t2_index = -1;
        for (int i = 0; i < 3; i++) {
            if ((t2.vertices_indices[i] == v0 && t2.vertices_indices[(i + 1) % 3] == v1) ||
                (t2.vertices_indices[i] == v1 && t2.vertices_indices[(i + 1) % 3] == v0)) {
                edge_t2_index = i;
                break;
            }
        }

        if (edge_t2_index == -1) return; // Edge not found in adjacent triangle

        // Neighboring triangles before the flip
        int t1_n0 = t1.neighbouring_triangles_indices[(edge_index + 1) % 3];
        int t1_n1 = t1.neighbouring_triangles_indices[(edge_index + 2) % 3];
        int t2_n0 = t2.neighbouring_triangles_indices[(edge_t2_index + 1) % 3];
        int t2_n1 = t2.neighbouring_triangles_indices[(edge_t2_index + 2) % 3];

        // Update triangle vertex indices
        t1.vertices_indices[0] = v2;
        t1.vertices_indices[1] = v3;
        t1.vertices_indices[2] = v0;

        t2.vertices_indices[0] = v3;
        t2.vertices_indices[1] = v2;
        t2.vertices_indices[2] = v1;

        // Update neighbors
        t1.neighbouring_triangles_indices[0] = adjacent_index;
        t1.neighbouring_triangles_indices[1] = t1_n0;
        t1.neighbouring_triangles_indices[2] = t2_n0;

        t2.neighbouring_triangles_indices[0] = triangle_index;
        t2.neighbouring_triangles_indices[1] = t1_n1;
        t2.neighbouring_triangles_indices[2] = t2_n1;

        // Correct neighboring triangle references
        if (t1_n0 != -1) {
            for (int i = 0; i < 3; i++)
                if (triangles[t1_n0].neighbouring_triangles_indices[i] == triangle_index)
                    triangles[t1_n0].neighbouring_triangles_indices[i] = adjacent_index;
        }

        if (t1_n1 != -1) {
            for (int i = 0; i < 3; i++)
                if (triangles[t1_n1].neighbouring_triangles_indices[i] == adjacent_index)
                    triangles[t1_n1].neighbouring_triangles_indices[i] = triangle_index;
        }
    }

    void insertPoint(const Vector& point) {
        // 1. Find the triangle containing the point
        int containingTriangle = -1;
        for (int i = 0; i < triangles.size(); i++) {
            if (triangles[i].pointInTriangle(point, vertices)) {
                containingTriangle = i;
                break;
            }
        }

        int new_vertex_index = vertices.size();
        addVertex(Vertex({ point[0], point[1], point[2] }, containingTriangle));

        if (containingTriangle != -1) {
            // 2. Split the triangle
            splitTriangle(containingTriangle, vertices.back());
        }
        else {
            // 3. If the point is outside, handle convex hull expansion
            handleOutsidePoint(new_vertex_index);
        }
    }

    std::vector<Vertex> vertices;
    std::vector<Triangle> triangles;

    // num vertices et num triangles pour raylib gulliver
};


int main() {

    
    return 0;
}

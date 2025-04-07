#include <cmath>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <thread>
#include <tuple>

#ifdef _OPENMP
#include <omp.h>
#endif


#define NEZ 19606

// TODO : check intégrité maillage au chargement : vérifier que chaque paire d'arrête entre commune entre 2 triangles soient inversées (a, b) (b, a)

static int max_threads = std::thread::hardware_concurrency();

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

Vector blueRedColorScale(double val, double min_val = -1, double max_val = 1) {
    // linearly map a value between min_val and max_val to a RGB values
    // min_val = blue, 0 = white, max_val = red

    double R{ 0 }, B{ 0 };

    if (val > 0) {
        R = 255 * (1 - (val / max_val));
        return Vector(255, R, R);
    }
    else {
        B = 255 * (1 - (val / min_val));
        return Vector(B, B, 255);
    }
}

Vector greyColorScale(double val, double min_val = 0, double max_val = 1) {
    // linearly map a value between min_val and max_val to RGB values
    // min_val = black, max_val = white

    double grey{ 255 * (val - min_val) / (max_val - min_val) };

    return Vector(grey, grey, grey);
}

Vector redColorScale(double val, double min_val = 0, double max_val = 1) {
    // linearly map a value between min_val and max_val to RGB values
    // min_val = white, max_val = red

    double R{ 255 * (1 - (val - min_val) / (max_val - min_val)) };

    return Vector(255, R, R);
}

class Vertex {
public:

    Vertex(std::initializer_list<double> position, int index_adjacent_triangle, bool is_virtual = false) :
        index_adjacent_triangle(index_adjacent_triangle), is_virtual(is_virtual) { updatePosition(position); }

    void updatePosition(std::initializer_list<double> position) {
        for (int i{ 0 }; i < 3; i++) this->position[i] = position.begin()[i];
    }

    Vector position;
    int index_adjacent_triangle;
    bool is_virtual;
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

    void check_virtual(const std::vector<Vertex>& vertices) {

        bool triangle_is_virtual{ false };

        for (int i{ 0 }; i < 3; i++)
            if (vertices[vertices_indices[i]].is_virtual)
                triangle_is_virtual = true;

        is_virtual = triangle_is_virtual;
    }

    bool pointInTriangle(const Vector& point, const std::vector<Vertex>& vertices) {

        check_virtual(vertices);

        Vector p{ vertices[vertices_indices[0]].position };
        Vector q{ vertices[vertices_indices[1]].position };
        Vector r{ vertices[vertices_indices[2]].position };

        // Check if point is inside the triangle using orientation predicate
        bool orientation1{ predicate_orientation(p, q, point) };
        bool orientation2{ predicate_orientation(q, r, point) };
        bool orientation3{ predicate_orientation(r, p, point) };

        // point is inside the triangle if all three orientations match (and the triangle is not connected to the infinite virtual point)
        return (orientation1 == orientation2) && (orientation2 == orientation3) && !is_virtual;
    }


    int vertices_indices[3];
    int neighbouring_triangles_indices[3];
    bool is_virtual{ false };
};

class Mesh {
public:

    void addVertex(const Vertex& Vertex) {
        vertices.push_back(Vertex);
        u.push_back(0);
        Lu.push_back(0);
    }

    void addTriangle(const Triangle& triangle) {
        triangles.push_back(triangle);
    }

    void writeOFF(
        const std::string& filename,
        bool enable_vertex_color = false,
        std::function<Vector(double, double, double)> color_scale = [](double val, double minv_val, double max_val) { return Vector(127, 127, 127); },
        double min_scale = 0,
        double max_scale = 1
    ) {

        check_virtual();

        std::ofstream ofs;
        ofs.open(filename);
        if (ofs.bad()) {
            std::cout << "Can not write file " << filename << std::endl;
            exit(1);
        }

        ofs << "OFF\n";
        ofs << vertices.size() - nb_virtual_vertices << " " << triangles.size() - nb_virtual_triangles << " 0\n";

        for (int i{ 0 }; i < vertices.size(); i++) {

            Vector col{ color_scale(u[i], min_scale, max_scale) };

            if (!vertices[i].is_virtual) {
                for (int j{ 0 }; j < 3; j++)
                    ofs << vertices[i].position[j] / scaleUp << " ";
                if (enable_vertex_color)
                    for (int j{ 0 }; j < 3; j++)
                        ofs << static_cast<int>(col[j]) << " ";
                ofs << "\n";
            }
        }

        for (int i{ 0 }; i < triangles.size(); i++) {
            if (!triangles[i].is_virtual) {
                ofs << "3";
                for (int j{ 0 }; j < 3; j++) {
                    ofs << " " << triangles[i].vertices_indices[j] - nb_virtual_vertices; // virtual vertex should be at the beginning of the array for this to be true
                }
                ofs << "\n";
            }
        }

        ofs.close();
        std::cout << filename << " successfully written" << std::endl;
    }

    void readOFF(const std::string& filename, double scaleUp = 1) {

        this->scaleUp = scaleUp;

        std::ifstream ifs;
        ifs.open(filename);
        if (ifs.bad()) {
            std::cout << "Can not read file " << filename << "\n";
            exit(1);
        }

        std::map<std::set<int>, std::pair<int, int>> map;

        int nb_vertices{ 0 }, nb_triangles{ 0 }, line_counter{ 0 }, bound1{ 2 }, bound2{ 0 }, bound3{ 0 };
        std::string line;
        while (std::getline(ifs, line)) { // parsing using order (comments not supported)
                
            std::vector<std::string> v_line{ split(line) };

            if (line_counter == 1) { // 2nd line (assuming there is a first line with OFF)

                nb_vertices = std::stoi(v_line[0]);
                bound2 = bound1 + nb_vertices;
                for (int i{ 0 }; i < nb_vertices; i++) addVertex(Vertex({ 0, 0, 0 }, 0));

                nb_triangles = std::stoi(v_line[1]);
                bound3 = bound2 + nb_triangles;
                for (int i{ 0 }; i < nb_triangles; i++) addTriangle(Triangle({ 0, 0, 0 }, { 0, 0, 0 }));
                // default neighbouring triangle will be the first one in case the mesh is not closed
            }
            else if (bound1 <= line_counter && line_counter < bound2) {

                int i{ line_counter - bound1 };
                vertices[i].updatePosition({ std::stod(v_line[0]) * scaleUp, std::stod(v_line[1]) * scaleUp, std::stod(v_line[2]) * scaleUp });
            }
            else if (bound2 <= line_counter && line_counter < bound3) {

                int i{ line_counter - bound2 };

                // we trust the file that the vertices indices are in the right (trigonometric) order
                // only support triangles reading
                triangles[i].updateVertices({ std::stoi(v_line[1]), std::stoi(v_line[2]), std::stoi(v_line[3])});
                // adjacent triangle for a Vertex will be the last adjacent triangle "seen"
                vertices[triangles[i].vertices_indices[0]].index_adjacent_triangle = i;
                vertices[triangles[i].vertices_indices[1]].index_adjacent_triangle = i;
                vertices[triangles[i].vertices_indices[2]].index_adjacent_triangle = i;

                for (int k{ 0 }; k < 3; k++) {

                    std::set<int> edge = std::set<int>{ triangles[i].vertices_indices[(k + 1) % 3], triangles[i].vertices_indices[(k + 2) % 3] };

                    if (map.find(edge) != map.end()) {

                        int i_bis{ map[edge].first };
                        int k_bis{ map[edge].second };
                        triangles[i].neighbouring_triangles_indices[k] = i_bis;
                        triangles[i_bis].neighbouring_triangles_indices[k_bis] = i;
                        map.erase(edge);
                    }
                    else {
                        map[edge] = std::pair<int, int>(i, k);
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
        for (int i{ 0 }; i < vertices.size(); i++)
            if (vertices[i].is_virtual)
                nb_virtual_vertices += 1;

        nb_virtual_triangles = 0;
        for (int i{ 0 }; i < triangles.size(); i++) {

            triangles[i].check_virtual(vertices);
            if (triangles[i].is_virtual) nb_virtual_triangles += 1;
        }
    }

    void Laplacian(int i_global) {

        Vertex& center_vertex = vertices[i_global];
        Vector A{ center_vertex.position };

        Lu[i_global] = 0;
        double Ai{ 0 };

        // NT = neighboring triangle
        int indexFirstNT{ center_vertex.index_adjacent_triangle };
        int indexNT{ indexFirstNT };
        
        do {
            Triangle& NT = triangles[indexNT];

            int i_local{ 0 };
            for (int k{ 0 }; k < 3; k++)
                if (NT.vertices_indices[k] == i_global) i_local = k;

            int j_local{ (i_local + 1) % 3 }; // j = vertex following the vertex i in NT
            int k_local{ (j_local + 1) % 3 };  // j = vertex before vertex i in NT
            
            int j_global{ NT.vertices_indices[j_local] };
            int indexFollowingNT{ NT.neighbouring_triangles_indices[j_local] };
            
            Triangle& previousNT = triangles[NT.neighbouring_triangles_indices[k_local]];
            int i_local_previousNT{ 0 };
            for (int k{ 0 }; k < 3; k++)
                if (previousNT.vertices_indices[k] == i_global) i_local_previousNT = k;
            int j_local_previousNT{ (i_local_previousNT + 1) % 3 };

            Vector B{ vertices[previousNT.vertices_indices[j_local_previousNT]].position };
            Vector C{ vertices[j_global].position };
            Vector D{ vertices[NT.vertices_indices[k_local]].position};

            Vector AC{ C - A };
            Vector AD{ D - A };
            double A_NT{ sqrt(cross(AC, AD).norm2()) };
            A_NT /= 2;
            Ai += A_NT / 3;

            Vector BC{ C - B };
            Vector BA{ A - B };
            double cot_alpha_ij{ dot(BC, BA) / sqrt(cross(BC, BA).norm2()) };

            Vector DA{ A - D };
            Vector DC{ C - D };
            double cot_beta_ij{ dot(DA, DC) / sqrt(cross(DA, DC).norm2()) };

            Lu[i_global] += (cot_alpha_ij + cot_beta_ij) * (u[j_global] - u[i_global]);

            indexNT = indexFollowingNT;

        } while (indexFirstNT != indexNT);

        Lu[i_global] /= 2 * Ai;
    }

    void computeLaplacian() {

        #ifdef _OPENMP
        omp_set_num_threads(max_threads);
        #endif
        
        #pragma omp parallel for schedule(dynamic, 1000)
        for (int i{ 0 }; i < vertices.size(); i++)
            Laplacian(i);
    }

    void splitTriangle(int triangle_index, Vertex& new_vertex) {
        
        new_vertex.index_adjacent_triangle = triangle_index;

        int new_vertex_index = vertices.size();
        addVertex(new_vertex);

        Triangle& triangle{ triangles[triangle_index] };

        // copy of the original triangle original properties
        int triangle_vertices_indices[3];
        int triangle_neighb_indices[3];
        for (int k{ 0 }; k < 3; k++) {
            triangle_vertices_indices[k] = triangle.vertices_indices[k];
            triangle_neighb_indices[k] = triangle.neighbouring_triangles_indices[k];
        }

        triangle.vertices_indices[2] = new_vertex_index; // reduce original triangle to one of the sub-triangles

        // update original triangle neighbours with the soon to be created triangles
        int new_triangle1_index = triangles.size();
        int new_triangle2_index = triangles.size() + 1;
        triangle.neighbouring_triangles_indices[0] = new_triangle1_index;
        triangle.neighbouring_triangles_indices[1] = new_triangle2_index;

        // create 2 new sub-triangles
        addTriangle(Triangle(
            { triangle_vertices_indices[1], triangle_vertices_indices[2], new_vertex_index },
            { new_triangle2_index, triangle_index, triangle_neighb_indices[0] }
        ));

        addTriangle(Triangle(
            { triangle_vertices_indices[2], triangle_vertices_indices[0], new_vertex_index },
            { triangle_index, new_triangle1_index, triangle_neighb_indices[1] }
        ));

        // update now invalid adjacent triangle for the vertex replaced in the original original triangle
        vertices[triangle_vertices_indices[2]].index_adjacent_triangle = new_triangle1_index;

        // update neighouring triangles with new sub-triangles indices
        
        // finding the original triangle in the first neighbouring triangle
        Triangle& firstNT{ triangles[triangle_neighb_indices[0]] };
        int original_triangle_local_index{ -1 };
        for (int j{ 0 }; j < 3; j++)
            if (firstNT.neighbouring_triangles_indices[j] == triangle_index)
                original_triangle_local_index = j;
        if (original_triangle_local_index == -1) throw("Sewing problem: original triangle not found in neighbouring triangle");
        firstNT.neighbouring_triangles_indices[original_triangle_local_index] = new_triangle1_index; // update

        // finding the original triangle in the second neighbouring triangle
        Triangle& secondNT{ triangles[triangle_neighb_indices[1]] };
        original_triangle_local_index = -1 ;
        for (int j{ 0 }; j < 3; j++)
            if (secondNT.neighbouring_triangles_indices[j] == triangle_index)
                original_triangle_local_index = j;
        if (original_triangle_local_index == -1) throw("Sewing problem: original triangle not found in neighbouring triangle");
        secondNT.neighbouring_triangles_indices[original_triangle_local_index] = new_triangle2_index; // update
    }

    void handlePointOutsideConvexHull(Vertex& new_vertex) { // not const, we need to update the neighbouring triangle index

        check_virtual();

        if (nb_virtual_vertices == 1) { // only work if infinite virtual vertex is used

            std::vector<std::tuple<int, int, int>> visibleBoundaryEdges;
            int virtual_triangle_to_split_index{ -1 };

            Vector& P{ new_vertex.position };

            // Find the visible boundary edges
            for (int i{ 0 }; i < triangles.size(); i++) {

                const Triangle& T{ triangles[i] };

                if (T.is_virtual) {

                    int virtual_vertex_local_index{ -1 };
                    for (int j{ 0 }; j < 3; j++)
                        if (T.vertices_indices[j] == 0) // assuming that the infinite vertex is first in the vector
                            virtual_vertex_local_index = j;
                    if (virtual_vertex_local_index == -1) throw("Virtual vertex not found in a virtual triangle");

                    int edgeT_index{ T.neighbouring_triangles_indices[virtual_vertex_local_index] };

                    //std::cout << "Virtual triangle : " << i << ", connected real triangle : " << edgeT_index << std::endl;

                    const Triangle& edgeT{ triangles[edgeT_index] }; // this triangle is on the edge of the convex hull

                    int virtual_triangle_local_index{ -1 };
                    for (int j{ 0 }; j < 3; j++)
                        if (edgeT.vertices_indices[j] == i)
                            virtual_triangle_local_index = j;
                    if (virtual_vertex_local_index == -1) throw("Virtual triangle not found in supposedly connected triangle");

                    std::tuple<int, int, int> boundaryEdge(
                        edgeT.vertices_indices[virtual_triangle_local_index],
                        edgeT.vertices_indices[(virtual_triangle_local_index + 1) % 3],
                        i // index of the virtual triangle connected to the convex hull edge (we want to flip it to expand the convex hull)
                    );

                    Vector A{ vertices[std::get<0>(boundaryEdge)].position };
                    Vector B{ vertices[std::get<1>(boundaryEdge)].position };

                    /*
                    std::cout << "Edge vertices: " << boundaryEdge.first << ", " << boundaryEdge.second << " (triangle " << i << ")" << std::endl;
                    std::cout << "A: " << A[0] << ", " << A[1] << ", " << A[2] << std::endl;
                    std::cout << "B: " << B[0] << ", " << B[1] << ", " << B[2] << std::endl;
                    std::cout << "P: " << P[0] << ", " << P[1] << ", " << P[2] << std::endl;
                    */

                    if (!predicate_orientation(A, B, P)) { // visibility test of boundary edge (A,B) from P ( triangle ABP must be of negative orientation)
                        //std::cout << "Visible" << std::endl;
                        visibleBoundaryEdges.push_back(boundaryEdge);
                        virtual_triangle_to_split_index = i; 
                    }
                }   
            }
            if (virtual_triangle_to_split_index == -1) throw("Did not find a boundary edge");

            splitTriangle(virtual_triangle_to_split_index, new_vertex); // last valid (connected to visible boundary) virtual triangle is splitted
            visibleBoundaryEdges.pop_back();

            /*
                bool isOutside = false;
                for (int j{ 0 }; j < 3; j++) {
                    int ai{ T.vertices_indices[j] };
                    int bi{ T.vertices_indices[(j + 1) % 3] };
                    Vector A{ vertices[ai].position };
                    Vector B{ vertices[bi].position };

                    // Check if the edge is "visible" from the new point
                    if (predicate_orientation(A, B, P)) {
                        visibleBoundaryEdges.push_back({ ai, bi });
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
            */
        }
        else {
            throw("No unique infinite vertex, handling new point outside the convex hull is not supported");
        }
    }


    void flipEdge(int triangle_index, int local_edge_index) {

        int T1_i{ triangle_index };
        Triangle& T1{ triangles[T1_i] };

        int v1_i{ T1.vertices_indices[(local_edge_index + 1) % 3] }; // start of selected edge
        int v2_i{ T1.vertices_indices[(local_edge_index + 2) % 3] }; // end of selected edge

        int T2_i{ T1.neighbouring_triangles_indices[local_edge_index] };
        Triangle& T2{ triangles[T2_i] }; // connected triangle, through edge (v1, v2)

        // Find the opposite vertex in the connected triangle
        int opposite_local_index{ -1 };
        for (int i{ 0 }; i < 3; i++)
            if (T2.vertices_indices[(i + 1) % 3] == v2_i && T2.vertices_indices[(i + 2) % 3] == v1_i) // edge vertex should be in the reversed order
                opposite_local_index = i;
        if (opposite_local_index == -1) throw("Sewing problem: cannot find corresponding edge in connected triangle");
        
        int v1bis_i{ T2.vertices_indices[opposite_local_index] };
        int v2bis_i{ T1.vertices_indices[local_edge_index] };

        // Flip the edge
        T1.vertices_indices[(local_edge_index + 1) % 3] = v1bis_i;
        T2.vertices_indices[(opposite_local_index + 1) % 3] = v2bis_i;

        // Update vertex adjacent triangle where it has potentially been invalidated
        vertices[v2_i].index_adjacent_triangle = T1_i;
        vertices[v1_i].index_adjacent_triangle = T2_i;

        // Getting surrounding triangles to update
        int T1N_i{ T1.neighbouring_triangles_indices[(local_edge_index + 2) % 3] };
        Triangle& T1N{ triangles[T1N_i] }; // triange opposite of v2 initially
        int T2N_i{ T2.neighbouring_triangles_indices[(opposite_local_index + 2) % 3] };
        Triangle& T2N{ triangles[T2N_i] }; // triange opposite of v1 initially

        // Find the opposite vertex in T1N
        int opposite_local_index_T1N{ -1 };
        for (int i{ 0 }; i < 3; i++)
            if (T1N.vertices_indices[(i + 1) % 3] == v1_i  && T1N.vertices_indices[(i + 2) % 3] == v2bis_i) 
                opposite_local_index_T1N = i;
        if (opposite_local_index_T1N == -1) throw("Sewing problem: cannot find corresponding edge in connected triangle");

        // Update the corresponding neighbouring triangle in T1N (from T1 to T2)
        T1N.neighbouring_triangles_indices[opposite_local_index_T1N] = T2_i;

        // Find the opposite vertex in T2N
        int opposite_local_index_T2N{ -1 };
        for (int i{ 0 }; i < 3; i++)
            if (T2N.vertices_indices[(i + 1) % 3] == v2_i && T2N.vertices_indices[(i + 2) % 3] == v1bis_i)
                opposite_local_index_T2N = i;
            
        if (opposite_local_index_T2N == -1) throw("Sewing problem: cannot find corresponding edge in connected triangle");

        // Update the corresponding neighbouring triangle in T2N (from T2 to T1)
        T2N.neighbouring_triangles_indices[opposite_local_index_T2N] = T1_i;

        // Updating T1 and T2 neighbouring triangles
        T1.neighbouring_triangles_indices[local_edge_index] = T2N_i; // from T2 to T2N
        T2.neighbouring_triangles_indices[opposite_local_index] = T1N_i; // from T1 to T1N
        T1.neighbouring_triangles_indices[(local_edge_index + 2) % 3] = T2_i; // from T1N to T2
        T2.neighbouring_triangles_indices[(opposite_local_index + 2) % 3] = T1_i; // from T2N to T1
    }

    void insertPoint(const Vector& new_vertex_pos) {
        
        // Find the triangle containing the point
        // TODO: implement ANT tactic to optimize this process
        int containingTriangle = -1;
        for (int i{ 0 }; i < triangles.size(); i++)
            if (triangles[i].pointInTriangle(new_vertex_pos, vertices))
                containingTriangle = i;
                
        Vertex new_vertex({ new_vertex_pos[0], new_vertex_pos[1], new_vertex_pos[2] }, containingTriangle);

        if (containingTriangle == -1) {
            // If the point is outside the convex hull, expand it
            handlePointOutsideConvexHull(new_vertex);
        }
        else {
            // Otherwise, split the triangle where the point belongs
            splitTriangle(containingTriangle, new_vertex);
        }
    }

    void toSphereSpace() {
        for (int i{ 0 }; i < vertices.size(); i++)
            vertices[i].position[2] = vertices[i].position[0] * vertices[i].position[0] + vertices[i].position[1] * vertices[i].position[1];
    }

    void debug() {
        std::cout << "=================== DEBUG ===================\nVertices : " << std::endl;
        for (int i{ 0 }; i < vertices.size(); i++) {
            std::cout << "Index : " << i << ", is virtual ? " << vertices[i].is_virtual;
            std::cout << ", Position : " << vertices[i].position[0] << ", " << vertices[i].position[1] << ", " << vertices[i].position[2];
            std::cout << ", Connected triangle index: " << vertices[i].index_adjacent_triangle << std::endl;
        }
        std::cout << "Triangles : " << std::endl;
        for (int i{ 0 }; i < triangles.size(); i++) {
            std::cout << "Index : " << i << ", is virtual ? " << triangles[i].is_virtual;
            std::cout << ", Vertices indices : " << triangles[i].vertices_indices[0] << ", " << triangles[i].vertices_indices[1] << ", " << triangles[i].vertices_indices[2];
            std::cout << ", Connected triangles indices : " << triangles[i].neighbouring_triangles_indices[0] << ", " << triangles[i].neighbouring_triangles_indices[1] << ", " << triangles[i].neighbouring_triangles_indices[2] << std::endl;
        }
    }

    int nb_virtual_vertices{ 0 }; // must be 0 or 1 (handling point outside convex hull won't work if 0)
    int nb_virtual_triangles{ 0 };

    double scaleUp{ 1 };

    std::vector<Vertex> vertices; // virtual vertex should be at the beginning of the array, otherwise, it will mess up the OFF file output and the handling of points outside the convex hull
    std::vector<Triangle> triangles;
    std::vector<double> u;
    std::vector<double> Lu;
};


int main() {

    /*
    
    // ***************************
    // *           TP1           *
    // ***************************
    

    // Tetrahedron
    Mesh tetrahedron;
    // It is not regular but ... simpler to draw in a standard(x, y, z) reference frame
    tetrahedron.addVertex(Vertex({ 0, 0, 0 }, 0));
    tetrahedron.addVertex(Vertex({ 1, 0, 0 }, 0));
    tetrahedron.addVertex(Vertex({ 0, 1, 0 }, 0));
    tetrahedron.addVertex(Vertex({ 0, 0, 1 }, 1));
    tetrahedron.addTriangle(Triangle({ 0, 1, 2 }, { 3, 2, 1 }));
    tetrahedron.addTriangle(Triangle({ 0, 3, 1 }, { 3, 0, 2 }));
    tetrahedron.addTriangle(Triangle({ 0, 2, 3 }, { 3, 1, 0 }));
    tetrahedron.addTriangle(Triangle({ 1, 3, 2 }, { 2, 0, 1 }));
    tetrahedron.writeOFF("off_files/tetrahedron_correct.off"); // we can visualize it in 3dviewer.net

    // Checking that we get the good results for the loading
    Mesh new_tetrahedron;
    new_tetrahedron.readOFF("off_files/tetrahedron_correct.off");
    new_tetrahedron.writeOFF("off_files/tetrahedron.off"); // this file must be identical to the correct one (read as text file)
    
    // We can compare this output to the tetrahedron creation to check if sewing is correct
    for (int i{ 0 }; i < new_tetrahedron.triangles.size(); i++)
        std::cout << new_tetrahedron.triangles[i].neighbouring_triangles_indices[0] << "|" << new_tetrahedron.triangles[i].neighbouring_triangles_indices[1] << "|" << new_tetrahedron.triangles[i].neighbouring_triangles_indices[2] << std::endl;

    
    // Square based pyramide
    Mesh pyramide;
    pyramide.addVertex(Vertex({ 0, 0, 0 }, 0));
    pyramide.addVertex(Vertex({ 1, 0, 0 }, 0));
    pyramide.addVertex(Vertex({ 1, 1, 0 }, 1));
    pyramide.addVertex(Vertex({ 0, 1, 0 }, 0));
    pyramide.addVertex(Vertex({ 0.5, 0.5, 1 }, 2));
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

    // We can compare this output to the tetrahedron creation to check if sewing is correct
    for (int i{ 0 }; i < new_pyramide.triangles.size(); i++)
        std::cout << new_pyramide.triangles[i].neighbouring_triangles_indices[0] << "|" << new_pyramide.triangles[i].neighbouring_triangles_indices[1] << "|" << new_pyramide.triangles[i].neighbouring_triangles_indices[2] << std::endl;


    // 2D bounding box
    Mesh bounding_box;
    bounding_box.addVertex(Vertex({ 0, 0, 0 }, 0));
    bounding_box.addVertex(Vertex({ 1, 0, 0 }, 0));
    bounding_box.addVertex(Vertex({ 1, 2, 0 }, 1));
    bounding_box.addVertex(Vertex({ 0, 2, 0 }, 0));
    bounding_box.addVertex(Vertex({ 0, 0, 10 }, 2, true));
    bounding_box.addTriangle(Triangle({ 0, 1, 3 }, { 1, 3, 4 }));
    bounding_box.addTriangle(Triangle({ 1, 2, 3 }, { 0, 5, 2 }));
    bounding_box.addTriangle(Triangle({ 2, 4, 3 }, { 1, 5, 3 }));
    bounding_box.addTriangle(Triangle({ 0, 3, 4 }, { 0, 2, 4 }));
    bounding_box.addTriangle(Triangle({ 0, 4, 1 }, { 0, 3, 5 }));
    bounding_box.addTriangle(Triangle({ 1, 4, 2 }, { 1, 4, 2 }));
    bounding_box.writeOFF("off_files/bounding_box_correct.off");
    // we can visualize it in 3dviewer.net (virtual vertices and triangles are not written in the OFF file)

    // Checking that we get the good results for the loading
    Mesh new_bounding_box;
    new_bounding_box.readOFF("off_files/bounding_box_correct.off");
    new_bounding_box.writeOFF("off_files/bounding_box.off"); // this file must be identical to the correct one (read as text file)

    // We can compare this output to the tetrahedron creation to check if sewing is correct
    for (int i{ 0 }; i < new_bounding_box.triangles.size(); i++)
        std::cout << new_bounding_box.triangles[i].neighbouring_triangles_indices[0] << "|" << new_bounding_box.triangles[i].neighbouring_triangles_indices[1] << "|" << new_bounding_box.triangles[i].neighbouring_triangles_indices[2] << std::endl;


    system("pause");

    Mesh queen; // we can visualize both in 3dviewer.net
    queen.readOFF("off_files/queen.off", 1000); // x1000 sur les longueurs, pour que delta_t ne soit pas trop petit, ...
    queen.writeOFF("off_files/queen_check.off");


    // ***************************
    // *           TP2           *
    // ***************************

    // ... car le pas temporel doit être suffisament petit pour que le schéma d'approximation numérique soit stable,
    // le pas temporel dépend de h²/alpha avec h la plus petite longueur du maillage  


    // Test diffusion thermique sur un mesh simple
    double T{ 100 }; // °C, température de la source
    new_pyramide.u[0] = T;
    double delta_t{ 1e-1 }; // sec, pas temporel
    double max_t{ 10 }; // sec, temps de simulation
    double alpha{ 1 }; // m²/s, diffusivité thermique, suppose que les positions sont exprimés en mètres

    for (double t{ 0 }; t < max_t; t += delta_t) {

        new_pyramide.computeLaplacian();

        for (int i{ 0 };  i < new_pyramide.u.size(); i++) {

            if (i != 0) new_pyramide.u[i] += delta_t * alpha * new_pyramide.Lu[i];
            std::cout << t+delta_t << "|" << new_pyramide.u[i] << "|" << new_pyramide.Lu[i] << std::endl;
        }
        std::cout << "---------------------------------" << std::endl;
    }

    system("pause");

    std::vector<Vector> curvature(queen.u.size());

    // Laplacien de x
    for (int i{ 0 }; i < queen.u.size(); i++) queen.u[i] = queen.vertices[i].position[0];
    queen.computeLaplacian();
    for (int i{ 0 }; i < queen.Lu.size(); i++) curvature[i][0] = queen.Lu[i];

    // Laplacien de y
    for (int i{ 0 }; i < queen.u.size(); i++) queen.u[i] = queen.vertices[i].position[1];
    queen.computeLaplacian();
    for (int i{ 0 }; i < queen.Lu.size(); i++) curvature[i][1] = queen.Lu[i];

    // Laplacien de z
    for (int i{ 0 }; i < queen.u.size(); i++) queen.u[i] = queen.vertices[i].position[2];
    queen.computeLaplacian();
    for (int i{ 0 }; i < queen.Lu.size(); i++) curvature[i][2] = queen.Lu[i];

    // Enregistrement de la valeur absolue de la courbure
    for (int i{ 0 }; i < queen.u.size(); i++) queen.u[i] = sqrt(curvature[i].norm2()) / 2;
    double min_curve{ *std::min_element(queen.u.begin(), queen.u.end()) };
    double max_curve{ *std::max_element(queen.u.begin(), queen.u.end()) };

    std::cout << "Minimum of absolute curvature: " << min_curve << std::endl;
    std::cout << "Maximum of absolute curvature: " << max_curve << std::endl;

    queen.writeOFF("off_files/queen_curvature_greyScale.off", true, greyColorScale, min_curve, max_curve);
    queen.writeOFF("off_files/queen_curvature_redScale.off", true,redColorScale, min_curve, max_curve);
    queen.writeOFF("off_files/queen_curvature_greyScale_saturated.off", true, greyColorScale, min_curve, 0.5 * max_curve);
    queen.writeOFF("off_files/queen_curvature_redScale_saturated.off", true, redColorScale, min_curve, 0.5 * max_curve);

    // Enregistrement de la courbure
    for (int i_global{ 0 }; i_global < queen.u.size(); i_global++) {

        int nb_positive{ 0 }, nb_negative{ 0 };

        Vertex& vertex{ queen.vertices[i_global] };
        Vector curv{ curvature[i_global] };

        // Parcours autour de vertex pour récupérer les normales des triangles adjacents
        // pour faire un vote majoritaire

        // NT = neighboring triangle
        int indexFirstNT{ vertex.index_adjacent_triangle };
        int indexNT{ indexFirstNT };

        Vector A{ vertex.position };

        do {
            Triangle& NT = queen.triangles[indexNT];

            int i_local{ 0 };
            for (int k{ 0 }; k < 3; k++)
                if (NT.vertices_indices[k] == i_global)
                    i_local = k;

            int j_local{ (i_local + 1) % 3 }; // j = vertex following the vertex i in NT
            int k_local{ (j_local + 1) % 3 };  // j = vertex before vertex i in NT

            int indexFollowingNT{ NT.neighbouring_triangles_indices[j_local] };

            int j_global{ NT.vertices_indices[j_local] };
            int k_global{ NT.vertices_indices[k_local] };

            Vector B{ queen.vertices[j_global].position };
            Vector C{ queen.vertices[k_global].position };

            Vector AB{ B - A };
            Vector AC{ C - A };

            Vector N{ cross(AB, AC) };
            
            if (dot(N, curv) > 0) nb_positive++;
            else nb_negative++;

            indexNT = indexFollowingNT;

        } while (indexFirstNT != indexNT);

        int sign{ nb_positive >= nb_negative ? 1 : -1 };

        queen.u[i_global] = sign * sqrt(curv.norm2()) / 2;
    }

    min_curve = *std::min_element(queen.u.begin(), queen.u.end());
    max_curve = *std::max_element(queen.u.begin(), queen.u.end());

    std::cout << "Minimum of signed curvature: " << min_curve << std::endl;
    std::cout << "Maximum of signed curvature: " << max_curve << std::endl;

    queen.writeOFF("off_files/queen_signed_curvature.off", true, blueRedColorScale, min_curve, max_curve);
    queen.writeOFF("off_files/queen_signed_curvature_saturated.off", true, blueRedColorScale, 0.5 * min_curve, 0.5 * max_curve);
    
    system("pause");
    
    // Test diffusion thermique sur queen.off
    delta_t = 75e-2; // si + grand (ex: 0.77), le schéma diverge
    int max_iter = 15e4; // nombre d'itérations
    max_t = delta_t * max_iter;
    for (int i{ 0 }; i < queen.u.size(); i++) queen.u[i] = 0; // reset car des valeurs de courbures étaient présentes dans u
    queen.u[NEZ] = T;

    for (double t{ 0 }; t < max_t; t += delta_t) {

        queen.computeLaplacian();

        for (int i{ 0 }; i < queen.u.size(); i++) if(i != NEZ) queen.u[i] += delta_t * alpha * queen.Lu[i];

        int iter{ static_cast<int>(t / delta_t) }; // itération actuelle

        if (iter % (max_iter / 1500) == 0) {

            double sum{ 0 };
            std::for_each(queen.u.begin(), queen.u.end(), [&sum](int el) { sum += el; });
            std::cout << iter << "/" << max_iter << " iterations | Check Laplacian tip nose: " << queen.Lu[NEZ] << " | Mean temperature: " << sum / queen.u.size() << std::endl;
        }

        if ((iter % (max_iter / 30) == 0) || (iter < 5000 && iter % 1000 == 0) || (iter < 1000 && iter % 200 == 0)) {
            // on enregistre toutes les 200it jusqu'à 1000it, toutes les 1000it jusqu'à 5000, puis toutes les 5000
            
            std::string str_t{ std::to_string(iter) };
            queen.writeOFF("off_files/temperature/queen_greyScale_" + str_t + ".off", true, greyColorScale, 0, T);
            queen.writeOFF("off_files/temperature/queen_redScale_" + str_t + ".off", true, redColorScale, 0, T);
        }
    }
    
    queen.writeOFF("off_files/temperature/queen_greyScale_final.off", true, greyColorScale, 0, T);
    queen.writeOFF("off_files/temperature/queen_redScale_final.off", true, redColorScale, 0, T);

    */


    // ***************************
    // *           TP3           *
    // ***************************

    Mesh small_triangulation;
    
    small_triangulation.addVertex(Vertex({0, 0, 10}, 1, true)); // "infinite" virtual vertex

    // Creating a first triangle and 3 virtual triangle to connect it to the infinite vertex
    small_triangulation.addVertex(Vertex({ 0, 0, 0 }, 0));
    small_triangulation.addVertex(Vertex({ 1, 0, 0 }, 0));
    small_triangulation.addVertex(Vertex({ 0, 1, 0 }, 0));
    small_triangulation.addTriangle(Triangle({ 1, 2, 3 }, { 2, 3, 1 }));
    small_triangulation.addTriangle(Triangle({ 1, 0, 2 }, { 2, 0, 3 }));
    small_triangulation.addTriangle(Triangle({ 2, 0, 3 }, { 3, 0, 1 }));
    small_triangulation.addTriangle(Triangle({ 0, 1, 3 }, { 0, 2, 1 }));

    // Insert point inside the convex hull    
    small_triangulation.insertPoint(Vector(0.25, 0.25, 0));

    // Insert point outside the convex hull to expand it
    small_triangulation.insertPoint(Vector(1, 1, 0));

    // Insert point inside the newly expanded convex hull
    small_triangulation.insertPoint(Vector(0.75, 0.75, 0));

    small_triangulation.debug(); // checked manually for sewing 

    // Flip an edge
    small_triangulation.flipEdge(4, 2);
    small_triangulation.debug(); // checked manually for sewing

    // Flip the resulting edge
    small_triangulation.flipEdge(4, 1);
    small_triangulation.debug(); // checked manually for sewing

    // Flip the resulting edge
    small_triangulation.flipEdge(4, 0);
    small_triangulation.debug();
    // Flip the resulting edge
    small_triangulation.flipEdge(4, 2);
    small_triangulation.debug(); // here we get back the same structure there was before any flip :
    // a flip is like rotating an edge and "pushing" triangles along !
    
    // Insert another point outside the convex hull, expansion need flips in this case
    //small_triangulation.insertPoint(Vector(-1, -1, 0));
    
    // Write the off file with the virtual infinite vertex
    small_triangulation.vertices[0].is_virtual = false;
    small_triangulation.writeOFF("small_triangulation.off");

    // triangulation.toSphereSpace(); Lift the triangulation on the paraboloid

    
    return 0;
}

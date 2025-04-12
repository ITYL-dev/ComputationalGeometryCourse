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
#include <filesystem>

#ifdef _OPENMP
#include <omp.h>
#endif


#define NEZ 19606 // vertices on the nose of the queen.off

// TODO : check mesh integrity at loading : check that all pair of "half-edges" are in the opposite sense (a,b) (b,a)

static int max_threads = std::thread::hardware_concurrency(); // number of threads for speeding up laplacian computing using openMP

double sqr(double a) { return a * a; }

std::vector<std::string> split(const std::string& line, const char& delim = ' ') {
    // split string into substring using separator delim

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
    // check if triangle pqr is direct in the plan normal to e_axis

    return dot(cross(q - p, r - p), e_axis) > 0; // true: positive, false: negative
}

Vector liftingOperator(const Vector& point) {

    Vector output{ point };
    output[2] = output[1] * output[1] + output[0] * output[0];

    return output;
}

bool predicate_in_circumscribed_cirle(const Vector& p, const Vector& q, const Vector& r, const Vector& s) { 
    // check if s is in the circumscribed circle to triangle pqr

    Vector phi_p{ liftingOperator(p) }, phi_q{ liftingOperator(q) }, phi_r{ liftingOperator(r) }, phi_s{ liftingOperator(s) };

    return dot(cross(phi_q - phi_p, phi_r - phi_p), phi_s - phi_p) < 0;
}

bool predicate_quad_convex(const Vector& A, const Vector& B, const Vector& C, const Vector& D, const Vector& N = Vector(0, 0, 1)) {
    // check if projection of quadrilateral ABCD on plan normal to N is convex
    // In the edge case where 3 points are aligned (a cross product is 0), return false

    Vector AB{ B - A };
    Vector BC{ C - B };
    Vector CD{ D - C };
    Vector DA{ A - D };

    double cross1{ dot(cross(AB, BC), N) };
    double cross2{ dot(cross(BC, CD), N) };
    double cross3{ dot(cross(CD, DA), N) };
    double cross4{ dot(cross(DA, AB), N) };

    // All cross products should have the same sign
    return (cross1 > 0 && cross2 > 0 && cross3 > 0 && cross4 > 0) ||
        (cross1 < 0 && cross2 < 0 && cross3 < 0 && cross4 < 0);
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
        // check if the triangle one of the vertices of the triangle is virtual
        // and update the triangle is_virtual property accordingly

        bool triangle_is_virtual{ false };

        for (int i{ 0 }; i < 3; i++)
            if (vertices[vertices_indices[i]].is_virtual)
                triangle_is_virtual = true;

        is_virtual = triangle_is_virtual;
    }

    bool pointInTriangle(const Vector& point, const std::vector<Vertex>& vertices) {
        // check if a point is in the triangle

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
        // write the mesh to an OFF file
        // coloring of the vertex can be enabled : mapping function u to a color using provided color scale

        check_virtual();

        std::filesystem::path file_path(filename);
        std::filesystem::path parent_dir{ file_path.parent_path() };
        std::filesystem::create_directories(parent_dir);

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
        // load a mesh from an OFF file (including sewing)
        // scaleUp can be used to multiply the coordinates of all vertices

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

                    std::set<int> edge{ triangles[i].vertices_indices[(k + 1) % 3], triangles[i].vertices_indices[(k + 2) % 3] };

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
        // check every triangle to see if they are virtual

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
        // compute the Laplacian of u[i_global] meaning Laplacian of u at vertices[i_global]

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
        // compute Laplacian of u at all vertices

        #ifdef _OPENMP
        omp_set_num_threads(max_threads);
        #endif
        
        #pragma omp parallel for schedule(dynamic, 1000)
        for (int i{ 0 }; i < vertices.size(); i++)
            Laplacian(i);
    }

    bool isEdgeDelaunay(int triangle_index, int edge_local_index) {
        // check if an edge can be flipped or not in Lawson algorithm (checking Delaunay property of edge
        // and convexity of quadrilateral made by the 2 triangle connected through the edge)

        const Triangle& T1{ triangles[triangle_index] };

        int T2_i{ T1.neighbouring_triangles_indices[edge_local_index] };
        const Triangle& T2{ triangles[T2_i] };

        int T1_local_index{ -1 }; // in T2
        for (int j{ 0 }; j < 3; j++)
            if (T2.neighbouring_triangles_indices[j] == triangle_index) T1_local_index = j;
        if (T1_local_index == -1) throw("Sewing problem : a triangle points to another triangle that does not point back to it");

        // if quadrilateral made by fusing T1 and T2 is not convex, we skip the edge (return true if edge is not delaunay)
        const Vector& V1{ vertices[T1.vertices_indices[edge_local_index]].position };
        const Vector& V2{ vertices[T1.vertices_indices[(edge_local_index + 1) % 3]].position };
        const Vector& V3{ vertices[T2.vertices_indices[T1_local_index]].position };
        const Vector& V4{ vertices[T2.vertices_indices[(T1_local_index + 1) % 3]].position };

        if (!predicate_quad_convex(V1, V2, V3, V4)) return true;

        bool edge_partially_delaunay_1{ !predicate_in_circumscribed_cirle(
            vertices[T1.vertices_indices[0]].position, vertices[T1.vertices_indices[1]].position, vertices[T1.vertices_indices[2]].position, // T1
            vertices[T2.vertices_indices[T1_local_index]].position // vertex opposite T1 in T2
        ) };

        bool edge_partially_delaunay_2{ !predicate_in_circumscribed_cirle(
            vertices[T2.vertices_indices[0]].position, vertices[T2.vertices_indices[1]].position, vertices[T2.vertices_indices[2]].position, // T2
            vertices[T1.vertices_indices[edge_local_index]].position // vertex opposite T2 in T1
        ) };

        return edge_partially_delaunay_1 && edge_partially_delaunay_2;
    }

    void splitTriangle(int triangle_index, Vertex& new_vertex) {
        // insert a new vertex in a triangle by splitting it into 3
        
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

    void flipEdge(int triangle_index, int local_edge_index, bool update_lawson_queue = false) {
        // flip an edge between 2 triangles, and optionnaly update lawson queue with the edges that are 
        // not locally Delaunay anymore because of the flip

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
            if (T1N.vertices_indices[(i + 1) % 3] == v1_i && T1N.vertices_indices[(i + 2) % 3] == v2bis_i)
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

        T1.check_virtual(vertices);
        T1N.check_virtual(vertices);
        T2.check_virtual(vertices);
        T2N.check_virtual(vertices);

        if (update_lawson_queue) {

            // Edge (v1bis, v2) connecting T1 and T2N after flip if edge is not virtual
            if (!T1.is_virtual && !T2N.is_virtual && !isEdgeDelaunay(T1_i, local_edge_index)) {
                std::set<int> edge{ T1.vertices_indices[(local_edge_index + 1) % 3], T1.vertices_indices[(local_edge_index + 2) % 3] };
                non_delaunay_edges.insert(edge);
            }

            // Edge (v1, v2bis) connecting T2 and T1N after flip if edge is not virtual
            if (!T2.is_virtual && !T1N.is_virtual && !isEdgeDelaunay(T2_i, opposite_local_index)) {
                std::set<int> edge{ T2.vertices_indices[(opposite_local_index + 1) % 3], T2.vertices_indices[(opposite_local_index + 2) % 3] };
                non_delaunay_edges.insert(edge);
            }

            // "constant" = did not change with flip

            // Edge (v2, v2bis) connecting T1 and its constant neighbour if edge is not virtual
            const Triangle& T1_constantN{ triangles[T1.neighbouring_triangles_indices[(local_edge_index + 1) % 3]] };
            if (!T1.is_virtual && !T1_constantN.is_virtual && !isEdgeDelaunay(T1_i, (local_edge_index + 1) % 3)) {
                std::set<int> edge{ T1.vertices_indices[(local_edge_index + 2) % 3], T1.vertices_indices[local_edge_index] };
                non_delaunay_edges.insert(edge);
            }

            // Edge (v1, v1bis) connecting T2 and its constant neighbour if edge is not virtual
            const Triangle& T2_constantN{ triangles[T2.neighbouring_triangles_indices[(opposite_local_index + 1) % 3]] };
            if (!T1.is_virtual && !T2_constantN.is_virtual && !isEdgeDelaunay(T2_i, (opposite_local_index + 1) % 3)) {
                std::set<int> edge{ T2.vertices_indices[(opposite_local_index + 2) % 3], T2.vertices_indices[opposite_local_index] };
                non_delaunay_edges.insert(edge);
            }
        }
    }

    void handlePointOutsideConvexHull(Vertex& new_vertex) { // not const, we need to update the neighbouring triangle index
        // insert a point that is outside the convex hull inserting it in a virtual triangle
        // and by expanding the convex hull using flips

        check_virtual();

        if (nb_virtual_vertices == 1) { // only work if infinite virtual vertex is used

            std::vector<std::tuple<int, int, int>> visibleBoundaryEdges;

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
                        if (edgeT.neighbouring_triangles_indices[j] == i)
                            virtual_triangle_local_index = j;
                    if (virtual_triangle_local_index == -1) throw("Virtual triangle not found in supposedly connected triangle");

                    std::tuple<int, int, int> boundaryEdge(
                        edgeT.vertices_indices[(virtual_triangle_local_index + 1) % 3],
                        edgeT.vertices_indices[(virtual_triangle_local_index + 2) % 3],
                        i // index of the virtual triangle connected to the convex hull edge (we want to flip it to expand the convex hull)
                    );

                    //std::cout << "Local index of virtual triangle in real triangle : " << virtual_vertex_local_index << " Edge : " << std::get<0>(boundaryEdge) << " " << std::get<1>(boundaryEdge) << std::endl;

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
                    }
                }   
            }
            if (visibleBoundaryEdges.empty()) throw("Did not find a boundary edge");

            std::tuple<int, int, int> starterEdge{ visibleBoundaryEdges.back() };
            visibleBoundaryEdges.pop_back();
            splitTriangle(std::get<2>(starterEdge), new_vertex); // last valid (connected to visible boundary) virtual triangle is splitted

            std::tuple<int, int, int> backwardEdge{ starterEdge };
            std::tuple<int, int, int> forwardEdge{ starterEdge };

            int nb_processed_edge{ 0 };
            
            // expanding convex hull by flipping specific edges
            while (nb_processed_edge < visibleBoundaryEdges.size()) {

                int previous_edge_i{ -1 };
                int next_edge_i{ -1 };

                for (int i{ 0 }; i < visibleBoundaryEdges.size(); i++) {

                    if (std::get<1>(visibleBoundaryEdges[i]) == std::get<0>(backwardEdge)) previous_edge_i = i;
                    if (std::get<0>(visibleBoundaryEdges[i]) == std::get<1>(forwardEdge)) next_edge_i = i;
                }

                if (previous_edge_i != -1) {

                    std::tuple<int, int, int> currentEdge{ visibleBoundaryEdges[previous_edge_i] };

                    int edge_to_flip_local_index{ -1 };
                    for (int j{ 0 }; j < 3; j++)
                        if (triangles[std::get<2>(currentEdge)].vertices_indices[j] == std::get<0>(currentEdge))
                            edge_to_flip_local_index = j;
                    if (edge_to_flip_local_index == -1) throw("Mismatch between bondary edge and corresponding virtual triangle");

                    flipEdge(std::get<2>(currentEdge), edge_to_flip_local_index);

                    backwardEdge = currentEdge;
                    nb_processed_edge++;
                }

                if (next_edge_i != -1) {

                    std::tuple<int, int, int> currentEdge{ visibleBoundaryEdges[next_edge_i] };

                    int edge_to_flip_local_index{ -1 };
                    for (int j{ 0 }; j < 3; j++)
                        if (triangles[std::get<2>(currentEdge)].vertices_indices[j] == std::get<1>(currentEdge))
                            edge_to_flip_local_index = j;
                    if (edge_to_flip_local_index == -1) throw("Mismatch between bondary edge and corresponding virtual triangle");

                    flipEdge(std::get<2>(currentEdge), edge_to_flip_local_index);

                    forwardEdge = currentEdge;
                    nb_processed_edge++;
                }
            }

            check_virtual();
        }
        else {
            throw("No unique infinite vertex, handling new point outside the convex hull is not supported");
        }
    }

    void insertPoint(const Vector& new_vertex_pos, bool in_delaunay = false) {
        // TODO: insert a point in a delaunay triangulation
        // --> update the lawson queue after inserting a point 
        // (without checking all the triangles if possible)
        // --> launch lawson from there (without buildLawsonQueue)
        
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

    std::pair<std::vector<double>, std::vector<int>> naiveTriangulationFromFile(const std::string& filename, double scaleUp = 1) {
        // read a file with x, y, z points in it and add them to a 2D (x, y, o) triangulation
        // returns z values and a mapping from the point indices to vertices indices, i.e. z[i] corresponds to vertices[index_map[i]]
        // coordinates are multiplied by scaleUp (used for paraboloid visualisation on a reasonable scale)

        std::cout << "Building naive triangulation from points in the .txt file ..." << std::endl;

        this->scaleUp = scaleUp;

        std::ifstream ifs;
        ifs.open(filename);
        if (ifs.bad()) {
            std::cout << "Can not read file " << filename << "\n";
            exit(1);
        }

        int nb_points{ 0 };
        std::vector<Vector> points;
        std::vector<int> index_map; // map point index to vertx index, i.e. points[i] corresponds to vertices[index_map[i]]
        std::vector<double> elevation;

        // read everything and put it in points
        std::string line;
        int line_counter{ 0 };
        while (std::getline(ifs, line)) {

            std::vector<std::string> v_line{ split(line) };

            if (line_counter == 0) {
                nb_points = std::stoi(v_line[0]);
                points.reserve(nb_points);
                index_map.reserve(nb_points);
                elevation.reserve(nb_points);
            }
            else if (line_counter <= nb_points) {
                Vector new_point(std::stod(v_line[0]), std::stod(v_line[1]), std::stod(v_line[2]));
                new_point = new_point * scaleUp;
                points.push_back(new_point);
            } 

            line_counter++;
        }

        addVertex(Vertex({ 0, 0, 10 }, 1, true)); // add a first vertex : "infinite" virtual vertex

        // Creating a first triangle from the points (in the plane), along with 3 virtual triangle to connect it to the infinite vertex

        Vector A{ points[0] };
        addVertex(Vertex({ A[0], A[1], 0 }, 0)); // we operate in the plane for now

        Vector B{ points[1] };
        addVertex(Vertex({ B[0], B[1], 0 }, 0)); // we operate in the plane for now

        // looking for a third point C so that triangle ABC faces upward (we would be very unlucky if it does not exist)
        int c_i{ 2 };
        Vector C{ points[c_i] };
        C[2] = 0;
        while (!predicate_orientation(A, B, C)) {
            c_i++;
            C = points[c_i]; // if C does not exist, will go out of bounds here
            C[2] = 0;
        }

        addVertex(Vertex({ C[0], C[1], 0 }, 0));
        addTriangle(Triangle({ 1, 2, 3 }, { 2, 3, 1 }));
        addTriangle(Triangle({ 1, 0, 2 }, { 2, 0, 3 }));
        addTriangle(Triangle({ 2, 0, 3 }, { 3, 0, 1 }));
        addTriangle(Triangle({ 0, 1, 3 }, { 0, 2, 1 }));

        // insert the points in the mesh through naive triangulation
        for (int i{ 0 }; i < nb_points; i++) {

            Vector P{ points[i] };
            elevation.push_back(P[2]);

            if (i == 0) index_map.push_back(1);
            else if (i == 1) index_map.push_back(2);
            else if (i == c_i) index_map.push_back(3);
            else {

                P[2] = 0;
                insertPoint(P);

                if (i < c_i) index_map.push_back(i + 2);
                else index_map.push_back(i + 1);
            }         
        }

        return std::pair<std::vector<double>, std::vector<int>>(elevation, index_map);
    }

    std::pair<int, int> translateEdgeRepresentation(const std::pair<int, int> edge) {
        // Translate an edge representation from (v1 index, v2 index) to (triangle index, local edge index)

        int v1_i{ edge.first };
        int v2_i{ edge.second };
        const Vertex& v1{ vertices[v1_i] };
        const Vertex& v2{ vertices[v2_i] };

        // find a triangle with both v1 and v2

        int NT_i{ v1.index_adjacent_triangle }; // NT = neighbouring triangle
        int first_NT_i{ NT_i };
        Triangle NT{ triangles[NT_i] };
        bool NT_with_edge{ false };
        int v1_local_edge{ -1 };
        int v2_local_edge{ -1 };

        do {

            //std::cout << "Current NT : " << NT_i << " Looking for edge : (" << v1_i << ", " << v2_i << ")" << std::endl;
            //std::cout << "NT vertices : ";

            v1_local_edge = -1;
            v2_local_edge = -1;
            for (int j{ 0 }; j < 3; j++) {
                //std::cout << "Local index : " << j << " Global Index : " << NT.vertices_indices[j] << " | ";
                if (NT.vertices_indices[j] == v1_i) v1_local_edge = j;
                if (NT.vertices_indices[j] == v2_i) v2_local_edge = j;
            }
            //std::cout << std::endl;
            if (v1_local_edge == -1) throw("Sewing problem : a vertex points to a triangle it does not belong to");

            if (v2_local_edge == -1) { // try next triangle
                NT_i = NT.neighbouring_triangles_indices[(v1_local_edge + 1) % 3];
                NT = triangles[NT_i];
            }
            else {
                NT_with_edge = true;
            }
        }
        while (!NT_with_edge && NT_i != first_NT_i);
        if (NT_i == first_NT_i && !NT_with_edge) throw("(v1, v2) is not an edge : no triangle in common");

        int opposite_NT_local_index{ 3 - v1_local_edge - v2_local_edge }; // the leftover local index that is not v1 or v2

        int opposite_NT_i{ NT.neighbouring_triangles_indices[opposite_NT_local_index] };
        const Triangle& opposite_NT{ triangles[opposite_NT_i] }; // triangle connected to NT through edge (v1, v2)

        // not necessary but allow to catch error
        int NT_local_edge{ -1 };
        for (int j{ 0 }; j < 3; j++)
            if (opposite_NT.neighbouring_triangles_indices[j] == NT_i) NT_local_edge = j;
        if (NT_local_edge == -1) throw("Sewing problem : a triangle points to another triangle that does not point back to it");

        return std::pair<int, int>(NT_i, opposite_NT_local_index);
    }

    void buildLawsonQueue() {
        // build the initial lawson queue with all edges that must be flipped for the lawson algorithm

        check_virtual();

        for (int i{ 0 }; i < triangles.size(); i++) {
            for (int j{ 0 }; j < 3; j++) {

                const Triangle& T{ triangles[i] };

                if (!T.is_virtual && !triangles[T.neighbouring_triangles_indices[j]].is_virtual) { // check if edge is not virtual
                    // edge in front of vertex j (local index) in triangle i
                    std::set<int> edge{ T.vertices_indices[(j + 1) % 3] , T.vertices_indices[(j + 2) % 3] };
                    if (!isEdgeDelaunay(i, j)) non_delaunay_edges.insert(edge); // using sets enforce unicity
                }
            }
        }
    }

    void delaunay(int max_iter = 1000000) {
        // lawson algorithm : flip all non locally delaunay edges, diagonal of convex quadrilateral
        // until there are none left. At each flip we check if edges should be added to the queue

        std::cout << "Improving the triangulation (Delaunay) using Lawson algorithm ..." << std::endl;

        int c{ 0 };

        buildLawsonQueue();

        /*std::cout << "Initial Lawson Queue :" << std::endl;
        for (std::set<int> edge : non_delaunay_edges) {
            std::set<int>::iterator it{ edge.begin() };
            int v1_i{ *it };
            int v2_i{ *std::next(it) };
            std::cout << v1_i << " " << v2_i << std::endl;
        }
        debug(true);*/

        while (!non_delaunay_edges.empty() && c < max_iter) { // hard limit on the number of loop, just in case (very big file, or bug)
            if (c % 1000 == 0) std::cout << "Loop : " << c << " Edges to flip : " << non_delaunay_edges.size() << std::endl;
            c++;

            // extract the first edge in the set
            std::set<int> edge{ *non_delaunay_edges.begin() };
            non_delaunay_edges.erase(non_delaunay_edges.begin());

            std::set<int>::iterator it{ edge.begin() };
            int v1_i{ *it };
            int v2_i{ *std::next(it) };

            // We use the edge (v1, v2) representation for storing edge in the queue,
            // instead of the (triangle index, local edge index) representation (which is more pratical),
            // because it may be invalidated by a flip.
            // Consequently, we need to translate from vertex to triangle representation of the edge :

            std::pair<int, int> edge_translated{ translateEdgeRepresentation(std::pair<int, int>(v1_i, v2_i)) };
            int triangle_index{ edge_translated.first };
            int edge_local_index{ edge_translated.second };

            flipEdge(triangle_index, edge_local_index, true); // flip the edge and update the queue

            /*std::cout << "Next Lawson Queue :" << std::endl;
            for (std::set<int> edge : non_delaunay_edges) {
                std::set<int>::iterator it{ edge.begin() };
                int v1_i{ *it };
                int v2_i{ *std::next(it) };
                std::cout << v1_i << " " << v2_i << std::endl;
            }
            debug(true);*/
        }
    }

    void toSphereSpace() {
        // lift the triangulation on a paraboloid (sphere space)

        for (int i{ 0 }; i < vertices.size(); i++)
            vertices[i].position = liftingOperator(vertices[i].position);
    }

    void debug(bool vertices_only = false) {
        // print the current triangulation
        std::cout << "====================================== DEBUG ======================================\nVertices : " << std::endl;
        for (int i{ 0 }; i < vertices.size(); i++) {
            std::cout << "Index : " << i << ", is virtual ? " << vertices[i].is_virtual;
            std::cout << ", Position : " << vertices[i].position[0] << ", " << vertices[i].position[1] << ", " << vertices[i].position[2];
            std::cout << ", Connected triangle index: " << vertices[i].index_adjacent_triangle << std::endl;
        }
        if (!vertices_only) {
            std::cout << "Triangles : " << std::endl;
            for (int i{ 0 }; i < triangles.size(); i++) {
                std::cout << "Index : " << i << ", is virtual ? " << triangles[i].is_virtual;
                std::cout << ", Vertices indices : " << triangles[i].vertices_indices[0] << ", " << triangles[i].vertices_indices[1] << ", " << triangles[i].vertices_indices[2];
                std::cout << ", Connected triangles indices : " << triangles[i].neighbouring_triangles_indices[0] << ", " << triangles[i].neighbouring_triangles_indices[1] << ", " << triangles[i].neighbouring_triangles_indices[2] << std::endl;
            }
        }
    }

    int nb_virtual_vertices{ 0 }; // must be 0 or 1 (handling point outside convex hull won't work if 0)
    int nb_virtual_triangles{ 0 };

    double scaleUp{ 1 };

    std::vector<Vertex> vertices; // virtual vertex should be at the beginning of the array, otherwise, it will mess up the OFF file output and the handling of points outside the convex hull
    std::vector<Triangle> triangles;
    std::vector<double> u; // u[i] is function u evaluated in vertices[i]
    std::vector<double> Lu; // Lu[i] is Laplacian of u in vertices[i]
    std::set<std::set<int>> non_delaunay_edges; // lawson algorithm queue
};


int main() {
    
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
    tetrahedron.writeOFF("files/TP1_out/tetrahedron_correct.off"); // we can visualize it in 3dviewer.net

    // Checking that we get the good results for the loading
    Mesh new_tetrahedron;
    new_tetrahedron.readOFF("files/TP1_out/tetrahedron_correct.off");
    new_tetrahedron.writeOFF("files/TP1_out/tetrahedron.off"); // this file must be identical to the correct one (read as text file)
    new_tetrahedron.debug(); // to check sewing by comparing to tetrahedron creation
    
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
    pyramide.writeOFF("files/TP1_out/pyramide_correct.off"); // we can visualize it in 3dviewer.net

    // Checking that we get the good results for the loading
    Mesh new_pyramide;
    new_pyramide.readOFF("files/TP1_out/pyramide_correct.off");
    new_pyramide.writeOFF("files/TP1_out/pyramide.off"); // this file must be identical to the correct one (read as text file)
    new_pyramide.debug(); // to check sewing by comparing to pyramide creation
    
    // 2D bounding box
    Mesh bounding_box;

    bounding_box.addVertex(Vertex({ 0, 0, 10 }, 2, true)); // virtual infini vertex
    bounding_box.addVertex(Vertex({ 0, 0, 0 }, 0));
    bounding_box.addVertex(Vertex({ 1, 0, 0 }, 0));
    bounding_box.addVertex(Vertex({ 1, 2, 0 }, 1));
    bounding_box.addVertex(Vertex({ 0, 2, 0 }, 0));
    bounding_box.addTriangle(Triangle({ 1, 2, 4 }, { 1, 5, 2 }));
    bounding_box.addTriangle(Triangle({ 2, 3, 4 }, { 4, 0, 3 }));
    bounding_box.addTriangle(Triangle({ 1, 0, 2 }, { 3, 0, 5 }));
    bounding_box.addTriangle(Triangle({ 3, 2, 0 }, { 2, 4, 3 }));
    bounding_box.addTriangle(Triangle({ 4, 3, 0 }, { 3, 5, 4 }));
    bounding_box.addTriangle(Triangle({ 4, 0, 1 }, { 2, 0, 4 }));
    bounding_box.writeOFF("files/TP1_out/bounding_box_correct.off");
    // we can visualize it in 3dviewer.net (virtual vertices and triangles are not written in the OFF file)
    
    // Checking that we get the good results for the loading
    Mesh new_bounding_box;
    new_bounding_box.readOFF("files/TP1_out/bounding_box_correct.off");
    new_bounding_box.writeOFF("files/TP1_out/bounding_box.off"); // this file must be identical to the correct one (read as text file)
    new_bounding_box.debug(); // to check sewing by comparing to bounding_box creation
    
    system("pause");

    Mesh queen; // we can visualize both in 3dviewer.net
    queen.readOFF("files/queen.off", 1000); // x1000 sur les longueurs, pour que delta_t ne soit pas trop petit, ...
    queen.writeOFF("files/TP1_out/queen_check.off");
    

    // ***************************
    // *           TP2           *
    // ***************************
   
    system("pause");
    
    // ... car le pas temporel doit être suffisament petit pour que le schéma d'approximation numérique soit stable,
    // le pas temporel dépend de h²/alpha avec h la plus petite longueur du maillage  


    // Test diffusion thermique sur un mesh simple
    double T{ 100 }; // °C, température de la source
    new_pyramide.u[0] = T;
    double delta_t{ 1e-1 }; // sec, pas temporel
    double max_t{ 10 }; // sec, temps de simulation
    double alpha{ 1 }; // m²/s, diffusivité thermique, suppose que les positions sont exprimés en mètres

    std::cout << "Timestep| Vertex temperature | Laplacian of temperature at vertex" << std::endl;

    for (double t{ 0 }; t < max_t; t += delta_t) {

        new_pyramide.computeLaplacian();

        for (int i{ 0 };  i < new_pyramide.u.size(); i++) {

            if (i != 0) new_pyramide.u[i] += delta_t * alpha * new_pyramide.Lu[i];
            std::cout << t+delta_t << "| " << new_pyramide.u[i] << " | " << new_pyramide.Lu[i] << std::endl;
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

    //queen.writeOFF("files/TP2_out/queen_curvature_greyScale.off", true, greyColorScale, min_curve, max_curve);
    queen.writeOFF("files/TP2_out/queen_curvature_redScale.off", true,redColorScale, min_curve, max_curve);
    //queen.writeOFF("files/TP2_out/queen_curvature_greyScale_saturated.off", true, greyColorScale, min_curve, 0.5 * max_curve);
    queen.writeOFF("files/TP2_out/queen_curvature_redScale_saturated.off", true, redColorScale, min_curve, 0.5 * max_curve);
    
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

    queen.writeOFF("files/TP2_out/queen_signed_curvature.off", true, blueRedColorScale, min_curve, max_curve);
    queen.writeOFF("files/TP2_out/queen_signed_curvature_saturated.off", true, blueRedColorScale, 0.5 * min_curve, 0.5 * max_curve);
    
    system("pause");
    
    // Test diffusion thermique sur queen.off
    delta_t = 75e-2; // si + grand (ex: 0.77), le schéma diverge
    int max_iter = 15e2; // nombre d'itérations (ne pas le diminuer sous 15e2, 15e4 montre une diffusion significative)
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
            // on enregistre toutes les 200it jusqu'à 1000it, toutes les 1000it jusqu'à 5000, puis toutes les 5000 (si max_iter = 15e4)
            
            std::string str_t{ std::to_string(iter) };
            queen.writeOFF("files/TP2_out/temperature/queen_greyScale_" + str_t + ".off", true, greyColorScale, 0, T);
            //queen.writeOFF("files/TP2_out/temperature/queen_redScale_" + str_t + ".off", true, redColorScale, 0, T);
        }   
    }
    
    queen.writeOFF("files/TP2_out/temperature/queen_greyScale_" + std::to_string(max_iter) + ".off", true, greyColorScale, 0, T);
    //queen.writeOFF("files/TP2_out/temperature/queen_redScale_" + std::to_string(max_iter) + ".off", true, redColorScale, 0, T);

    
    // ***************************
    // *           TP3           *
    // ***************************

    system("pause");

    // Test on a small handcrafted mesh to check our features

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

    //small_triangulation.debug(); // checked manually for sewing 

    // Flip an edge
    small_triangulation.flipEdge(4, 2);
    //small_triangulation.debug(); // checked manually for sewing

    // Flip the resulting edge
    small_triangulation.flipEdge(4, 1);
    //small_triangulation.debug(); // checked manually for sewing

    // Flip the resulting edge
    small_triangulation.flipEdge(4, 0);
    //small_triangulation.debug();

    // Flip the resulting edge
    small_triangulation.flipEdge(4, 2);
    //small_triangulation.debug(); // here we get back the same structure there was before any flip :
    // a flip is like rotating an edge and "pushing" triangles along !
    
    // Insert others points outside the convex hull, whose expansion needs flips in this case
    small_triangulation.insertPoint(Vector(-1, -1, 0));
    small_triangulation.insertPoint(Vector(-0.8, 0, 0));
    small_triangulation.insertPoint(Vector(1.5, 1.2, 0));
    small_triangulation.insertPoint(Vector(-1.1, 1.1, 0));
    small_triangulation.insertPoint(Vector(0.5, -0.5, 0));
    small_triangulation.insertPoint(Vector(1, -1, 0));
    small_triangulation.insertPoint(Vector(0, 1.5, 0));
    
    // Write the off file with optionnaly the virtual infinite vertex
    small_triangulation.vertices[0].is_virtual = true; // set to false to visualize infinite virtual vertex and virtual triangles
    small_triangulation.writeOFF("files/TP3_out/small_triangulation.off");

    // Use lawson algorithm to make the triangulation delaunay    
    small_triangulation.delaunay();
    small_triangulation.writeOFF("files/TP3_out/small_triangulation_delaunay.off");

    // Lift the triangulation on the paraboloid 
    small_triangulation.toSphereSpace(); // to lift the triangulation on the paraboloid
    small_triangulation.writeOFF("files/TP3_out/small_triangulation_lifted.off");

    // Test on one of the provided point clouds
    std::string filename("alpes_poisson");

    Mesh triangulation;
    std::pair<std::vector<double>, std::vector<int>> res{ triangulation.naiveTriangulationFromFile("files/" + filename + ".txt", 0.00015) };
    triangulation.writeOFF("files/TP3_out/" + filename + "_naive.off");

    // Apply Lawson algorithm
    triangulation.delaunay();
    triangulation.writeOFF("files/TP3_out/" + filename + "_delaunay.off");

    // Lift on the paraboloid
    triangulation.toSphereSpace();
    triangulation.writeOFF("files/TP3_out/" + filename + "_lifted.off");

    // Reapply elevation
    std::vector<double> elevation{ res.first };
    std::vector<int> index_map{ res.second };
    for (int i{ 0 }; i < elevation.size(); i++) triangulation.vertices[index_map[i]].position[2] = elevation[i];
    triangulation.writeOFF("files/TP3_out/" + filename + "_elevation.off");
    
    return 0;
}

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
                for (int j{ 0 }; j < 3; j++)
                    ofs << static_cast<int>(col[j]) << " ";
                ofs << "\n";
            }
        }

        for (int i{ 0 }; i < triangles.size(); i++) {
            if (!triangles[i].is_virtual) {
                ofs << "3";
                for (int j{ 0 }; j < 3; j++) {
                    ofs << " " << triangles[i].vertices_indices[j];
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

            bool triangle_is_virtual{ false };

            for (int j{ 0 }; j < 3; j++)
                if (vertices[triangles[i].vertices_indices[j]].is_virtual)
                    triangle_is_virtual = true;

            triangles[i].is_virtual = triangle_is_virtual;
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

    int nb_virtual_vertices{ 0 };
    int nb_virtual_triangles{ 0 };

    double scaleUp{ 1 };

    std::vector<Vertex> vertices;
    std::vector<Triangle> triangles;
    std::vector<double> u;
    std::vector<double> Lu;
};


int main() {

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

    queen.writeOFF("off_files/queen_curvature_greyScale.off", greyColorScale, min_curve, max_curve);
    queen.writeOFF("off_files/queen_curvature_redScale.off", redColorScale, min_curve, max_curve);
    queen.writeOFF("off_files/queen_curvature_greyScale_saturated.off", greyColorScale, min_curve, 0.5 * max_curve);
    queen.writeOFF("off_files/queen_curvature_redScale_saturated.off", redColorScale, min_curve, 0.5 * max_curve);

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

    queen.writeOFF("off_files/queen_signed_curvature.off", blueRedColorScale, min_curve, max_curve);
    queen.writeOFF("off_files/queen_signed_curvature_saturated.off", blueRedColorScale, 0.5 * min_curve, 0.5 * max_curve);
    
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
            queen.writeOFF("off_files/temperature/queen_greyScale_" + str_t + ".off", greyColorScale, 0, T);
            queen.writeOFF("off_files/temperature/queen_redScale_" + str_t + ".off", redColorScale, 0, T);
        }
    }
    
    queen.writeOFF("off_files/temperature/queen_greyScale_final.off", greyColorScale, 0, T);
    queen.writeOFF("off_files/temperature/queen_redScale_final.off", redColorScale, 0, T);
    
    return 0;
}

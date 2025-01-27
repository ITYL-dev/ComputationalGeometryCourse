#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <string>
#include <initializer_list>

#define DTYPE double


template <size_t dim>
class Vertex {
public:

    Vertex(std::initializer_list<DTYPE> position, size_t index_adjacent_triangle): index_adjacent_triangle(index_adjacent_triangle) {
        for (size_t i{ 0 }; i < dim; i++) this->position[i] = position.begin()[i];
    }

    DTYPE position[dim];
    size_t index_adjacent_triangle;
};

class Triangle {
public:

    Triangle(std::initializer_list<size_t> vertices_indices, std::initializer_list<size_t> neighbouring_triangles_indices) {
        for (size_t i{ 0 }; i < 3; i++) {
            this->vertices_indices[i] = vertices_indices.begin()[i];
            this->neighbouring_triangles_indices[i] = neighbouring_triangles_indices.begin()[i];
        }
    }

    size_t vertices_indices[3];
    size_t neighbouring_triangles_indices[3];
};

template <size_t dim>
class Mesh {
public:

    void addVertex(const Vertex<dim> & vertex) {
        vertices.push_back(vertex);
    }

    void addTriangle(const Triangle& triangle) {
        triangles.push_back(triangle);
    }

    void writeOFF(const char* filename) {

        if (dim == 3) {

            std::ofstream ofs;
            ofs.open(filename);
            if (ofs.bad()) {
                std::cout << "Can not write file " << filename << std::endl;
                exit(1);
            }

            ofs << "OFF\n";
            ofs << vertices.size() << " " << triangles.size() << " 0\n";

            for (size_t i{ 0 }; i < vertices.size(); i++) {
                for (size_t j{ 0 }; j < dim; j++) {

                    ofs << vertices[i].position[j];
                    if (j < dim - 1) ofs << " ";
                }
                ofs << "\n";
            }

            for (size_t i{ 0 }; i < triangles.size(); i++) {
                ofs << "3";
                for (size_t j{ 0 }; j < 3; j++) {
                    ofs << " " << triangles[i].vertices_indices[j];
                }
                ofs << "\n";
            }

            ofs.close();
            std::cout << filename << " successfully written" << std::endl;
        }

        else std::cout << "Can not write " << dim << "D mesh as OFF file" << std::endl;
    }

    void readOFF(const char* filename) {

        if (dim == 3) {

            std::ifstream ifs;
            ifs.open(filename);
            if (ifs.bad()) {
                std::cout << "Can not read file " << filename << "\n";
                exit(1);
            }

            
            size_t nb_vertices{ 0 }, nb_triangles{ 0 }, counter{ 0 }, bound1{ 3 }, bound2{ 0 }, bound3{ 0 };
            std::string string;
            while (ifs >> string) { // parsing using order (variable "counter"), therefore, comments (starting with "#") are not handled
                if (counter == 1) {
                    nb_vertices = std::stoull(string);
                    bound2 = bound1 + dim * nb_vertices;
                    for (size_t i{ 0 }; i < nb_vertices; i++) addVertex(Vertex<dim>({ 0, 0, 0 }, 0));
                }
                else if (counter == 2) {
                    nb_triangles = std::stoull(string);
                    bound3 = bound2 + 4 * nb_triangles; // parsing only handle triangles = faces of 3 vertices / edges (4 = 3 + 1)
                    for (size_t i{ 0 }; i < nb_triangles; i++) addTriangle(Triangle({0, 0, 0}, {0, 0, 0}));
                }
                else if (bound1 < counter && counter <= bound2) {
                    std::cout << string << "|vertices|";
                    size_t i{ (counter - bound1) / dim };
                    size_t j{ (counter - bound1) % dim };
                    vertices[i].position[j] = std::stod(string);
                }
                else if (bound2 < counter && counter <= bound3) {
                    std::cout << string << "|triangles|";
                    size_t j{ (counter - bound2) % 4 };
                    if (j > 0) {
                        size_t i{ (counter - bound2) / 4 };
                        triangles[i].vertices_indices[j] = std::stod(string);
                    }
                }

                counter += 1;
            }

            ifs.close();
            std::cout << filename << " successfully loaded" << std::endl;
        }

        else std::cout << "Can not load OFF file for a " << dim << "D mesh" << std::endl;
    }

    std::vector<Vertex<dim>> vertices;
    std::vector<Triangle> triangles;
};


int main() {

    // Tetrahedron
    Mesh<3> tetrahedron;
    // It is not regular but ... simpler to draw in a standard(x, y, z) reference frame
    tetrahedron.addVertex(Vertex<3>({ 0, 0, 0 }, 0));
    tetrahedron.addVertex(Vertex<3>({ 1, 0, 0 }, 0));
    tetrahedron.addVertex(Vertex<3>({ 0, 1, 0 }, 0));
    tetrahedron.addVertex(Vertex<3>({ 0, 0, 1 }, 1));
    tetrahedron.addTriangle(Triangle({ 0, 1, 2 }, { 3, 2, 1 }));
    tetrahedron.addTriangle(Triangle({ 0, 3, 1 }, { 3, 0, 2 }));
    tetrahedron.addTriangle(Triangle({ 0, 2, 3 }, { 3, 1, 0 }));
    tetrahedron.addTriangle(Triangle({ 1, 3, 2 }, { 2, 0, 1 }));
    tetrahedron.writeOFF("tetrahedron.off");
    tetrahedron.readOFF("tetrahedron.off");

    return 0;
}

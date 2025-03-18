#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <string>
#include <initializer_list>
#include <set>
#include <map>
#include <utility>

#define DTYPE double // if modifying datatype, remember to swap std::stod for the correct function in Mesh::readOFF

typedef unsigned long long size_t ; // should already be defined in cstdib,
// but I redefine it in case size_t is not unsigned long long (even it's very unlikely),
// to ensure compatibility with std::stoull (see Mesh::readOFF)


template <size_t dim>
class Vertex {
public:

    Vertex(std::initializer_list<DTYPE> position, size_t index_adjacent_triangle, bool is_virtual = false ):
        index_adjacent_triangle(index_adjacent_triangle), is_virtual(is_virtual) {
        for (size_t i{ 0 }; i < dim; i++) this->position[i] = position.begin()[i];
    }

    DTYPE position[dim];
    size_t index_adjacent_triangle;
    bool is_virtual;
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
    bool is_virtual{ false };
};

template <size_t dim>
class Mesh {
public:

    void addVertex(const Vertex<dim> & vertex) {
        vertices.push_back(vertex);
        // check_virtual() ?
    }

    void addTriangle(const Triangle& triangle) {
        triangles.push_back(triangle);
        // check_virtual() ?
    }

    void writeOFF(const char* filename) {

        check_virtual();

        if (dim == 3) { // TODO : handle 2D meshes

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
                    for (size_t j{ 0 }; j < dim; j++) {
                        ofs << vertices[i].position[j];
                        if (j < dim - 1) ofs << " ";
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

        else std::cout << "Can not write " << dim << "D mesh as OFF file" << std::endl;
    }

    void readOFF(const char* filename) {

        if (dim == 3) { // TODO : handle 2D meshes

            std::ifstream ifs;
            ifs.open(filename);
            if (ifs.bad()) {
                std::cout << "Can not read file " << filename << "\n";
                exit(1);
            }
            
            std::map<std::set<size_t>, std::pair<size_t, size_t>> map;
   
            size_t nb_vertices{ 0 }, nb_triangles{ 0 }, counter{ 0 }, bound1{ 3 }, bound2{ 0 }, bound3{ 0 };
            std::string string;
            while (ifs >> string) { // parsing using order (using variable "counter")
                                    // therefore, some things like comments (starting with "#") are not handled

                if (counter == 1) { // 2nd line

                    nb_vertices = std::stoull(string);
                    bound2 = bound1 + dim * nb_vertices;
                    for (size_t i{ 0 }; i < nb_vertices; i++) addVertex(Vertex<dim>({ 0, 0, 0 }, 0));
                }
                else if (counter == 2) { // 3rd line

                    nb_triangles = std::stoull(string);
                    bound3 = bound2 + 4 * nb_triangles; // parsing only handle triangles = faces of 3 vertices / edges (4 = 3 + 1)
                    for (size_t i{ 0 }; i < nb_triangles; i++) addTriangle(Triangle({0, 0, 0}, {0, 0, 0}));
                }
                else if (bound1 < counter && counter <= bound2) {

                    // std::cout << string << "|vertices|";
                    size_t i{ (counter - bound1 - 1) / dim };
                    size_t j{ (counter - bound1 - 1) % dim };
                    vertices[i].position[j] = std::stod(string);
                }
                else if (bound2 < counter && counter <= bound3) {

                    // std::cout << string << "|triangles|";
                    size_t j{ (counter - bound2 - 1) % 4 };

                    if (j > 0) {

                        j = j - 1;
                        size_t i{ (counter - bound2 - 1) / 4 };
                        triangles[i].vertices_indices[j] = std::stoull(string); // we trust the file that the vertices indices are in the right (trigonometric) order
                        vertices[triangles[i].vertices_indices[j]].index_adjacent_triangle = i;

                        if (j == 2) {

                            std::vector<std::set<size_t>> edges;
                            edges.reserve(3);
                            edges.emplace_back(std::set<size_t>{ triangles[i].vertices_indices[1], triangles[i].vertices_indices[2] }); // relative index 0
                            edges.emplace_back(std::set<size_t>{ triangles[i].vertices_indices[2], triangles[i].vertices_indices[0] }); // relative index 1
                            edges.emplace_back(std::set<size_t>{ triangles[i].vertices_indices[0], triangles[i].vertices_indices[1] }); // relative index 2
                            // relative index k = index (inside of the Triangle class) of a vertex, the edge opposite to the vertex and the triangle connected to opposite edge

                            for (size_t k{ 0 }; k < 3; k++) {
                                if (map.count(edges[k])) { // 1 ...

                                    size_t i_bis{ map[edges[k]].first };
                                    triangles[i].neighbouring_triangles_indices[k] = i_bis;
                                    triangles[i_bis].neighbouring_triangles_indices[map[edges[k]].second] = i;
                                }
                                else { // ... or 0
                                    map[edges[k]] = std::pair<size_t, size_t>(i, k);
                                }
                            }
                        }
                    }
                }

                counter += 1;
            }

            ifs.close();
            std::cout << filename << " successfully loaded" << std::endl;
        }

        else std::cout << "Can not load OFF file for a " << dim << "D mesh" << std::endl;
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
    tetrahedron.writeOFF("tetrahedron_correct.off"); // we can visualize it in 3dviewer.net

    // Checking that we get the good results for the loading
    Mesh<3> new_tetrahedron;
    new_tetrahedron.readOFF("tetrahedron_correct.off");
    new_tetrahedron.writeOFF("tetrahedron.off"); // this file must be identical to the correct one (read as text file)

    // We can compare this output to the tetrahedron creation to check if sewing is correct
    for (int i{ 0 }; i < new_tetrahedron.triangles.size(); i++)
        std::cout << new_tetrahedron.triangles[i].neighbouring_triangles_indices[0] << "|" << new_tetrahedron.triangles[i].neighbouring_triangles_indices[1] << "|" << new_tetrahedron.triangles[i].neighbouring_triangles_indices[2] << std::endl;


    // Square based pyramide
    Mesh<3> pyramide;
    pyramide.addVertex(Vertex<3>({ 0, 0, 0 }, 0));
    pyramide.addVertex(Vertex<3>({ 1, 0, 0 }, 0));
    pyramide.addVertex(Vertex<3>({ 1, 1, 0 }, 1));
    pyramide.addVertex(Vertex<3>({ 0, 1, 0 }, 0));
    pyramide.addVertex(Vertex<3>({ 0.5, 0.5, 1 }, 2));
    pyramide.addTriangle(Triangle({ 0, 1, 3 }, { 1, 3, 4 }));
    pyramide.addTriangle(Triangle({ 1, 2, 3 }, { 2, 0, 5 }));
    pyramide.addTriangle(Triangle({ 2, 4, 3 }, { 3, 1, 5 }));
    pyramide.addTriangle(Triangle({ 0, 3, 4 }, { 2, 4, 0 }));
    pyramide.addTriangle(Triangle({ 0, 4, 1 }, { 5, 0, 3 }));
    pyramide.addTriangle(Triangle({ 1, 4, 2 }, { 2, 1, 4 }));
    pyramide.writeOFF("pyramide_correct.off"); // we can visualize it in 3dviewer.net

    // Checking that we get the good results for the loading
    Mesh<3> new_pyramide;
    new_pyramide.readOFF("pyramide_correct.off");
    new_pyramide.writeOFF("pyramide.off"); // this file must be identical to the correct one (read as text file)

    // We can compare this output to the tetrahedron creation to check if sewing is correct
    for (int i{ 0 }; i < new_pyramide.triangles.size(); i++)
        std::cout << new_pyramide.triangles[i].neighbouring_triangles_indices[0] << "|" << new_pyramide.triangles[i].neighbouring_triangles_indices[1] << "|" << new_pyramide.triangles[i].neighbouring_triangles_indices[2] << std::endl;


    // 2D bounding box
    Mesh<3> bounding_box;
    bounding_box.addVertex(Vertex<3>({ 0, 0, 0 }, 0));
    bounding_box.addVertex(Vertex<3>({ 1, 0, 0 }, 0));
    bounding_box.addVertex(Vertex<3>({ 1, 2, 0 }, 1));
    bounding_box.addVertex(Vertex<3>({ 0, 2, 0 }, 0));
    bounding_box.addVertex(Vertex<3>({ 0, 0, 0 }, 2, true));
    bounding_box.addTriangle(Triangle({ 0, 1, 3 }, { 1, 3, 4 }));
    bounding_box.addTriangle(Triangle({ 1, 2, 3 }, { 0, 5, 2 }));
    bounding_box.addTriangle(Triangle({ 2, 4, 3 }, { 1, 5, 3 }));
    bounding_box.addTriangle(Triangle({ 0, 3, 4 }, { 0, 2, 4 }));
    bounding_box.addTriangle(Triangle({ 0, 4, 1 }, { 0, 3, 5 }));
    bounding_box.addTriangle(Triangle({ 1, 4, 2 }, { 1, 4, 2 }));
    bounding_box.writeOFF("bounding_box_correct.off");
    // we can visualize it in 3dviewer.net (virtual virtices and triangles are not written in the OFF file)

    // Checking that we get the good results for the loading
    Mesh<3> new_bounding_box;
    new_bounding_box.readOFF("bounding_box_correct.off");
    new_bounding_box.writeOFF("bounding_box.off"); // this file must be identical to the correct one (read as text file)
    
    // We can compare this output to the tetrahedron creation to check if sewing is correct
    for (int i{ 0 }; i < new_bounding_box.triangles.size(); i++)
        std::cout << new_bounding_box.triangles[i].neighbouring_triangles_indices[0] << "|" << new_bounding_box.triangles[i].neighbouring_triangles_indices[1] << "|" << new_bounding_box.triangles[i].neighbouring_triangles_indices[2] << std::endl;

    
    Mesh<3> queen; // we can visualize both in 3dviewer.net
    queen.readOFF("queen.off");
    queen.writeOFF("queen_check.off");

    return 0;
}

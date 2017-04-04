#ifndef FLUIDQUADTREE
#define FLUIDQUADTREE

// #define _USE_MATH_DEFINES

#include <cstdlib>
#include <unordered_map>
#include <cmath>
#include <memory>
#include <tuple>
#include <vector>

#include "array2_utils.h"
#include "vec.h"

enum FacePosition {
	RIGHT, LEFT, TOP, BOTTOM
};

// Four directions are corresponding to lower left, lower right, upper left, upper right.
enum DiagonalDirection {
	SW, SE, NW, NE
};

// Most member fields should be public so struct is used. 
struct Cell {
    
	int depth, i, j;
    
	// Cell constructor.
    Cell(){}
    Cell(int depth_, int i_, int j_) :depth(depth_), i(i_), j(j_) {}

	friend std::ostream& operator<<(std::ostream& os, const Cell& c);
	bool operator==(const Cell& rhs) {
		return depth == rhs.depth && i == rhs.i && j == rhs.j;
	}
    
};

// The default postion of a face respect to a cell is Right for u and Top for v
// except for T junctions.
struct Face {
    
	// Location of the owner cell.
	int depth, i, j;

	// The face is on the right, left, top, or bottom of that cell
	FacePosition position;

	// Face constructor.
	Face(int depth_, int i_, int j_, FacePosition pos_) : i(i_), j(j_), depth(depth_), position(pos_) {
	}

	bool operator==(const Face& y) {
		return i == y.i && j == y.j && depth == y.depth && position == y.position;
	}
    
};


// The default position of a node respect to a cell is NE.
struct Node {
    
	// Location of the owner cell.
	int depth, i, j;

	// The node is on the SE, SW, NE, NW of its owner cell.
	DiagonalDirection position;

	// Node constructor.
	Node(int depth_, int i_, int j_, DiagonalDirection pos_) : i(i_), j(j_), depth(depth_), position(pos_) {
	}

	bool operator==(const Node& y) {
		return i == y.i && j == y.j && depth == y.depth && position == y.position;
	}
    
	bool operator< (const Node& y) const {
		return std::tie(i, j, depth, position) < std::tie(y.i, y.j, y.depth, y.position);
	}
    
};

// Define a quadtree structure and its functions. 
struct FluidQuadTree {
    
    // FluidQuadTree constructor.
    FluidQuadTree(float width, Array2f liquid_phi_);
    
    //tree overall manipulation
    void reindex();
    
    //boolean tests on tree
    bool is_cell_active(const Cell& c);
    bool is_cell_in_bounds(const Cell& c);
    bool is_leaf_cell(const Cell& c);
    
    //Accessors by neighbour relationships
    Cell above(const Cell& c);
    Cell below(const Cell& c);
    Cell right(const Cell& c);
    Cell left(const Cell& c);
    Cell child(const Cell& c); // The lower left child.
    
    Cell get_active_parent(Cell c);
    
    Cell get_child(const Cell& c, const DiagonalDirection& dir);
    Cell get_leaf_cell(const Cell& c, const DiagonalDirection& dir);
    Cell get_parent(const Cell& c);
    
    // Some geometry lookups
    Vec2f get_cell_centre(const Cell& c);
    Vec2f get_face_centre(const Face& f);
    Vec2f get_node_position(const Node& n);
    
    // Cell indexing info
    int get_cell_index(const Cell& c);
    int get_level_dims(int depth);
    float get_cell_width(int depth);
    
    std::vector<Array2i> cell_markers;
	std::vector<Array2i> cellInds;
    
    // All cells, faces and nodes on the quadtree. 
	std::vector<Cell> leaf_cells;
	std::vector<Face> qt_faces;
	std::vector<Node> nodes;
    
    //// The maps below linking faces, nodes and cells together.
    
	// The order of faces is (right, left, top, bottom). The dimension
    // of the map is n_cell * 4;
    Array2i cell_to_face_map;
    // The order of nodes is (SW, SE, NW, NE). The dimension
    // of the map is n_cell * 4;
    Array2i cell_to_node_map;
	// The order of faces is (right, left, top, bottom). 
	// If a face doesn't exist(T junction), then it is stored as -1.
	std::vector<std::vector<int>> node_to_face_map;
	// The order of nodes is (right, left) or (top, bottom).
	std::unordered_map<int, Vec2i> face_to_node_map;
	// The order of cells is (right, left) or (top, bottom).
    // For a T junction face, the cells on each side are also
    // in the preceding order.
	std::unordered_map<int, std::pair<std::vector<int>, std::vector<int>>> face_to_cell_map;

    
    int leaf_cell_count;
    int active_cell_count;
	int face_count;
	int node_count;
    int max_depth;
    
private:
    
    Array2f liquid_phi;
    float domain_width;
    float dx;
    int ni, nj;
    
};

#endif 

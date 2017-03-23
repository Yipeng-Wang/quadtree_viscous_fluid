#ifndef FLUIDQUADTREE
#define FLUIDQUADTREE

// #define _USE_MATH_DEFINES

#include <cstdlib>
#include <iostream>
#include <map>
#include <cmath>
#include <memory>
#include <tuple>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/SparseCore>

#include "array2_utils.h"
#include "vec.h"

typedef Eigen::SparseMatrix<float> SpMat;

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
	Cell() {}
	Cell(int depth_, int i_, int j_) :depth(depth_), i(i_), j(j_) {
	}

	friend std::ostream& operator<<(std::ostream& os, const Cell& c);
	bool operator==(const Cell& rhs) {
		return depth == rhs.depth && i == rhs.i && j == rhs.j;
	}
};

struct FluidFace {
	// Location of the owner cell.
	int depth, i, j;

	// The face is on the right, left, top, or bottom of that cell
	FacePosition position;

	// FluidFace constructor.
	FluidFace(int depth_, int i_, int j_, FacePosition pos_) : i(i_), j(j_), depth(depth_), position(pos_) {
	}

	bool operator==(const FluidFace& y) {
		return i == y.i && j == y.j && depth == y.depth && position == y.position;
	}
};

struct FluidNode {
	// Location of the owner cell.
	int depth, i, j;

	// The node is on the SE, SW, NE, NW of its owner cell.
	DiagonalDirection position;

	// FluidNode constructor.
	FluidNode(int depth_, int i_, int j_, DiagonalDirection pos_) : i(i_), j(j_), depth(depth_), position(pos_) {
	}

	bool operator==(const FluidNode& y) {
		return i == y.i && j == y.j && depth == y.depth && position == y.position;
	}
	bool operator< (const FluidNode& y) const {
		return std::tie(i, j, depth, position) < std::tie(y.i, y.j, y.depth, y.position);
	}
};

// Define a quadtree structure and its functions. 
struct FluidQuadTree {
    Array2f liquid_phi;
    
    // H matrix transform the liquid velocities on regular grid
    // to that on a quadtree grid.
    SpMat H;
    
	std::vector<Array2i> cells;
	std::vector<Array2i> cellInds;
    
    // signed distance at each level.
    std::vector<Array2f> all_phi;
    
	std::vector<Cell> leaf_cells;
	std::vector<FluidFace> velocity_faces;
	
	// The code is changed to include nodes at T junctions.  
	std::vector<FluidNode> nodes;
		
	// The order of faces is (right, left, top, bottom). 
    Array2i cell_to_face_map;

	// The order of faces is (right, left, top, bottom). 
	// If a face doesn't exist(T junction), then it is stored as -1.
	std::vector<std::vector<int>> node_to_face_map;

	// The order of nodes is (SW, SE, NW, NE).
	//std::vector<std::vector<int>>
    Array2i cell_to_node_map;

	// The order of nodes is (right, left) or (top, bottom).
	std::map<int, Vec2i> face_to_node_map;

	// The order of cells is (right, left) or (top, bottom). For a T junction face,
	// the cells on each side are also in the preceding order.
	std::map<int, std::pair<std::vector<int>, std::vector<int>>> face_to_cell_map;

	int leaf_cell_count, active_cell_count;
	int face_count;
	int node_count;
    
	// FluidQuadTree constructor.
	FluidQuadTree(float width, Array2f liquid_phi_);
    
	//tree overall manipulation
	void reindex();

	//boolean tests on tree
	bool is_cell_active(const Cell& c);
	bool is_cell_in_bounds(const Cell& c);
	bool is_leaf_cell(const Cell& c);

	//Accessors by neighbour relationships
	Cell above(Cell c);
	Cell below(Cell c);
	Cell left(Cell c);
	Cell right(Cell c);
	Cell child(Cell c); // The lower left child.

	Cell get_active_parent(Cell c);
	Cell get_child(const Cell &c, const DiagonalDirection &dir);
	Cell get_leaf_cell(const Cell &c, const DiagonalDirection &dir);
	Cell get_parent(const Cell& c);

	// Some geometry lookups
	Vec2f get_cell_centre(const Cell& c);
	Vec2f get_face_centre(const FluidFace& f);
	Vec2f get_node_position(const FluidNode& n);

	// Cell indexing info
	int get_cell_index(const Cell& c);
	int get_level_dims(int depth);
	float get_cell_width(int depth);
    
    void get_trans_matrix();
    
    // Data for viscosity solve, c_vol is cell volumes and n_vol is node volumes.
    Array2f u_vol, v_vol, c_vol, n_vol;
    Array2f viscosity;
    
    void compute_reg_grid_weights();

private:
    float domain_width;
    int max_depth;
    float dx;
    int ni, nj;

};

#endif 

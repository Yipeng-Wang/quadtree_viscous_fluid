//
//  viscositysolver.hpp
//  SurfaceFluid
//
//  Created by Yipeng Wang on 2017-03-30.
//  Copyright © 2017 Yipeng Wang. All rights reserved.
//

#ifndef viscositysolver_hpp
#define viscositysolver_hpp

#include <Eigen/Sparse>
#include <Eigen/SparseCore>

#include "fluidquadtree.h"

typedef Eigen::SparseMatrix<float> SpMat;
typedef Eigen::VectorXf Vecf;


class VisSolver {
    
    float EPSILON;
    
    // Functions for viscosity solve.
    void get_velocities();
    
    void compute_reg_grid_weights();
    void compute_trans_matrices();
    void compute_deformation_operator();
    void compute_factor_matrix();
    void compute_volume_matrices();
    void compute_tau_vis();
    void compute_rhs();
    
    void solve_viscosity(float dt);
    void get_grid_taus();
    void get_node_taus();
    
    
    float dx;
    int ni, nj;
    int n_u, n_v, n_t11, n_t22, n_t12;
    
    // c_vol is cell volumes and n_vol is node volumes.
    Array2f u_vol, v_vol, c_vol, n_vol;
    
    // H matrix transform the liquid velocities on regular grid
    // to that on a quadtree grid.
    SpMat H_u, H_v;
    SpMat H;
    
    // Data used in the linear system solver.
    SpMat D;
    SpMat factor;
    SpMat M_uv, M_t;
    SpMat M_vis;
    Vecf rhs;
    SpMat vis_operator;
    Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper> cg_solver;
    
    std::vector<Face> u_reg_faces;  // u faces on the regular grids.
    std::vector<Face> v_reg_faces;  // v faces on the regular grids.
    
    // Map i + ni *j to idx in u_reg_faces.
    std::unordered_map<int, int> u_reg_idx_to_face;
    
    // Map i + ni * j to idx in v_reg_faces.
    std::unordered_map<int, int> v_reg_idx_to_face;
    
    Vecf reg_u;
    Vecf reg_v;
    Vecf reg_uv;
    Vecf M_reg_uv;                  // Control volumes of regular grids.
    
    std::vector<Face> u_tree_faces;  // u faces on the tree grids.
    std::vector<Face> v_tree_faces;  // v faces on the tree grids.
    Vecf tree_u;
    Vecf tree_v;
    
    int n_u_tree, n_v_tree;
    int n_u_reg, n_v_reg;
    
    // Link u's with faces.
    std::vector<int> u_to_face;
    std::unordered_map<int, int> face_to_u;
    
    // Link v's with faces.
    std::vector<int> v_to_face;
    std::unordered_map<int, int> face_to_v;
    
    // Link tau11's with grids.
    std::vector<int> tau11_to_cell;
    std::unordered_map<int, int> grid_to_tau11;
    
    // Link tau22's with grids.
    std::vector<int> tau22_to_cell;
    std::unordered_map<int, int> grid_to_tau22;
    
    // Link tau12's with nodes.
    std::vector<int> tau12_to_node;
    std::unordered_map<int, int> node_to_tau12;
    
public:
    // VisSolver constructors
    VisSolver(Array2f u_, Array2f v_, Array2f viscosity_, float width,
              Array2f liquid_phi_, Array2f solid_phi_, float dt);

    FluidQuadTree tree;
    
    Array2f u; // (ni+1) * nj 2d array.
    Array2f v; // ni * (nj+1) 2d array.
    Array2f viscosity; // ni * nj 2d array. It's stored at cell centre.
    Array2f liquid_phi; // Phi is at the cell centre.
    Array2f solid_phi; // Phi is at the nodal position.
};



#endif /* viscositysolver_hpp */

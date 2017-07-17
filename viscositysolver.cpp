//
//  viscositysolver.cpp
//  SurfaceFluid
//
//  Created by Yipeng Wang on 2017-03-30.
//  Copyright © 2017 Yipeng Wang. All rights reserved.
//

#include "util.h"
#include "pcgsolver/blas_wrapper.h"
#include "pcgsolver/pcg_solver.h"
#include "pcgsolver/sparse_matrix.h"
#include "viscositysolver.h"


using namespace std;

void compute_volume_fractions(const Array2f& levelset, Array2f& fractions,
                              Vec2f fraction_origin, int subdivision);

// Build one large matrix by adding m1 to its upper left corner and m2 to its
// lower right corner.
SpMat build_large_diagnol(SpMat m1, SpMat m2);

// Copy an Eigen sparse matrix to Robert's SparseMatrixd.
static SparseMatrixd copy_sparse_matrix(SpMat m);


VisSolver::VisSolver(Array2f u_, Array2f v_, Array2f viscosity_, float width,
                     Array2f liquid_phi_, Array2f solid_phi_) :
u(u_), v(v_), viscosity(viscosity_), liquid_phi(liquid_phi_),
solid_phi(solid_phi_), tree(FluidQuadTree(width, liquid_phi_)) {
  
  ni = liquid_phi.ni;
  nj = liquid_phi.nj;
  dx = width / (float) ni;
  
  
  EPSILON = 1E-2 * sqr(dx);
  
  u_vol.resize(ni+1,nj);
  v_vol.resize(ni,nj+1);
  c_vol.resize(ni,nj);
  n_vol.resize(ni+1,nj+1);
  
}

void VisSolver::solve_viscosity(float dt) {
  get_velocities();
  get_grid_taus();
  get_node_taus();
  
  n_u = (int)u_to_face.size();
  n_v = (int)v_to_face.size();
  n_t11 = (int)tau11_to_cell.size();
  n_t22 = (int)tau22_to_cell.size();
  n_t12 = (int)tau12_to_node.size();
  
  int n_uv = n_u + n_v;
  int n_t = n_t11 + n_t22 + n_t12;
  
  // Generate the deformation operator.
  D.resize(n_t, n_uv);
  compute_deformation_operator();
  D = 0.5 * D;
  
  // Compute factoring matrix.
  factor.resize(n_t, n_t);
  factor.reserve(n_t);
  compute_factor_matrix();
  
  // Control volume matrices for uv's and
  M_uv.resize(n_uv, n_uv);
  M_t.resize(n_t, n_t);
  M_uv.reserve(n_uv);
  M_t.reserve(n_t);
  
  compute_volume_matrices();
  
  M_vis.resize(n_t, n_t);
  M_vis.reserve(n_t);
  
  compute_tau_vis();
  
  vis_operator = -2 * D.transpose() * (factor * M_vis * M_t) * D;
  
  SpMat sym_mat = M_uv - dt * vis_operator;
  
  rhs.resize(n_u + n_v);
  
  compute_rhs();
  
  // Solve the linear system using a conjugate gradient solver.
  cg_solver.compute(sym_mat);
  Vecf tree_uv = cg_solver.solve(rhs);
  cout << "#iterations:     " << cg_solver.iterations() << endl;
  cout << "estimated error: " << cg_solver.error()      << endl;
  
  /***********************************************
   // Solve the linear system with Robert's solver.
   SparseMatrixd r_matrix= copy_sparse_matrix(sym_mat);
   vector<double> r_rhs;
   for (int i = 0; i < rhs.size(); ++i) {
   r_rhs.push_back(rhs[i]);
   }
   double res_out;
   int iter_out;
   PCGSolver<double> solver;
   vector<double> r_velocities(rhs.size());
   solver.solve(r_matrix, r_rhs, r_velocities, res_out, iter_out);
   Vecf tree_uv(r_velocities.size());
   for (int i = 0; i < r_velocities.size(); ++i) {
   tree_uv[i] = r_velocities[i];
   }
   ************************************************/
  
  // M multiply u reg.
  Vecf Mu_reg = reg_uv.cwiseProduct(M_reg_uv) + (dt * H.transpose()
                                                 * vis_operator * tree_uv);
  
  // This data only used to update uv's with zero control volumes.
  Vecf fs_uv = H.transpose() * tree_uv;
  
  for (int i = 0; i < n_u_reg; ++i) {
    Face& u_f = u_reg_faces[i];
    
    if (M_reg_uv[i] < EPSILON) {
      u(u_f.i + 1, u_f.j) = fs_uv[i];
    } else {
      u(u_f.i + 1, u_f.j) = Mu_reg[i] / M_reg_uv[i];
    }
  }
  
  for (int i = 0; i < n_v_reg; ++i) {
    Face& v_f = v_reg_faces[i];
    int idx = n_u_reg + i;
    
    if (M_reg_uv[idx] < EPSILON) {
      v(v_f.i, v_f.j + 1) = fs_uv[idx];
    } else {
      v(v_f.i, v_f.j + 1) = Mu_reg[idx] / M_reg_uv[idx];
    }
  }
  
}

// Compute the cell, node, uv face weights on the regular grid.
void VisSolver::compute_reg_grid_weights() {
  compute_volume_fractions(liquid_phi, c_vol, Vec2f(-0.5, -0.5), 2);
  compute_volume_fractions(liquid_phi, n_vol, Vec2f(-1, -1), 2);
  compute_volume_fractions(liquid_phi, u_vol, Vec2f(-1, -0.5), 2);
  compute_volume_fractions(liquid_phi, v_vol, Vec2f(-0.5, -1), 2);
}

// Get actual u and v faces which the solver runs on, and store them in
// u_tree_faces and v_tree_faces.
void VisSolver::get_velocities() {
  compute_reg_grid_weights();
  
  Array2c u_state(ni + 1, nj, (const char&)0);
  Array2c v_state(ni, nj + 1, (const char&)0);
  const int SOLID = 1;
  const int FLUID = 0;
  
  // Just determine if the face position is inside the wall.
  for (int j = 0; j < nj; ++j) {
    for (int i = 0; i < ni + 1; ++i) {
      if (i - 1 < 0 || i >= ni
          || (solid_phi(i, j + 1) + solid_phi(i, j)) / 2 <= 0)
        u_state(i, j) = SOLID;
      else
        u_state(i, j) = FLUID;
    }
  }
  
  for (int j = 0; j < nj + 1; ++j) {
    for (int i = 0; i < ni; ++i) {
      if (j - 1 < 0 || j >= nj
          || (solid_phi(i + 1, j) + solid_phi(i, j)) / 2 <= 0)
        v_state(i, j) = SOLID;
      else
        v_state(i, j) = FLUID;
    }
  }
  
  // Start checking which u, v faces should be used in the solver.
  for (int k = 0; k < tree.face_count; ++k) {
    
    Face& f = tree.qt_faces[k];
    if (f.depth != tree.max_depth - 1) {
      if (f.position == LEFT || f.position == RIGHT) {
        if (f.position == LEFT) {
          f = Face(f.depth, f.i-1, f.j, RIGHT);
        }
        u_tree_faces.push_back(f);
        face_to_u[k] = (int)u_to_face.size();
        u_to_face.push_back(k);
      } else {
        if (f.position == BOTTOM) {
          f = Face(f.depth, f.i, f.j-1, TOP);
        }
        v_tree_faces.push_back(f);
        face_to_v[k] = (int)v_to_face.size();
        v_to_face.push_back(k);
      }
      continue;
    }
    
    int i = f.i, j = f.j;
    
    if (f.position == RIGHT)  {
      i++;
      
      if (u_state(i, j) == SOLID) {
        continue;
      }
      if (u_vol(i, j) > EPSILON || c_vol(i, j) > EPSILON
          || c_vol(i-1, j) > EPSILON || n_vol(i, j+1) > EPSILON
          || n_vol(i, j) > EPSILON) {
        u_tree_faces.push_back(f);
        face_to_u[k] = (int)u_to_face.size();
        u_to_face.push_back(k);
      }
      continue;
    }
    
    if (f.position == TOP) {
      j++;
      
      if (v_state(i, j) == SOLID) {
        continue;
      }
      if (v_vol(i, j) > EPSILON || c_vol(i, j) > EPSILON
          || c_vol(i, j-1) > EPSILON || n_vol(i+1, j) > EPSILON
          || n_vol(i, j) > EPSILON) {
        v_tree_faces.push_back(f);
        face_to_v[k] = (int)v_to_face.size();
        v_to_face.push_back(k);
      }
    }
  }
  compute_trans_matrices();
}


// Calculate two transforming matrices H_u and H_v, and the actual u's and
// v's used in the solver.
void VisSolver::compute_trans_matrices() {
  
  // Set up the matrix H_u transforming u's on regular grids to tree grids.
  
  n_u_tree = (int)u_tree_faces.size();
  H_u.resize(n_u_tree, n_u_tree);
  H_u.setIdentity();
  
  u_reg_faces = u_tree_faces;
  
  int n_cur = n_u_tree;
  vector<Face>* cur_level_fu = &u_tree_faces;
  vector<Face>* next_level_fu;
  
  for (int k = 0; k < tree.max_depth-1; ++k) {
    next_level_fu = new vector<Face>();
    unordered_map<int, int> next_idx_to_u;
    SpMat H_level(n_cur, ni * nj);
    
    for (int i = 0; i < n_cur; ++i) {
      Face& f = (*cur_level_fu)[i];
      
      int idx;
      if (f.depth != k) {
        idx = tree.get_face_idx(f);
        if (!next_idx_to_u.count(idx)) {
          next_idx_to_u[idx] = (int)next_level_fu->size();
          next_level_fu->push_back(f);
        }
        H_level.insert(i, next_idx_to_u[idx]) = 1;
        continue;
      }
      
      Cell c(f.depth, f.i, f.j);
      
      Cell upper_c = tree.get_child(c, NE);
      Cell lower_c = tree.get_child(c, SE);
      
      
      Face ul_f(upper_c.depth, upper_c.i-1, upper_c.j, RIGHT);
      Face ll_f(lower_c.depth, lower_c.i-1, lower_c.j, RIGHT);
      
      Face um_f(upper_c.depth, upper_c.i, upper_c.j, RIGHT);
      Face lm_f(lower_c.depth, lower_c.i, lower_c.j, RIGHT);
      
      Face ur_f(upper_c.depth, upper_c.i+1, upper_c.j, RIGHT);
      Face lr_f(lower_c.depth, lower_c.i+1, lower_c.j, RIGHT);
      
      
      idx = tree.get_face_idx(ul_f);
      if (!next_idx_to_u.count(idx)) {
        next_idx_to_u[idx] = (int)next_level_fu->size();
        next_level_fu->push_back(ul_f);
      }
      H_level.insert(i, next_idx_to_u[idx]) = 0.125;
      
      
      idx = tree.get_face_idx(ll_f);
      if (!next_idx_to_u.count(idx)) {
        next_idx_to_u[idx] = (int)next_level_fu->size();
        next_level_fu->push_back(ll_f);
      }
      H_level.insert(i, next_idx_to_u[idx]) = 0.125;
      
      
      idx = tree.get_face_idx(um_f);
      if (!next_idx_to_u.count(idx)) {
        next_idx_to_u[idx] = (int)next_level_fu->size();
        next_level_fu->push_back(um_f);
      }
      H_level.insert(i, next_idx_to_u[idx]) = 0.25;
      
      
      idx = tree.get_face_idx(lm_f);
      if (!next_idx_to_u.count(idx)) {
        next_idx_to_u[idx] = (int)next_level_fu->size();
        next_level_fu->push_back(lm_f);
      }
      H_level.insert(i, next_idx_to_u[idx]) = 0.25;
      
      
      idx = tree.get_face_idx(ur_f);
      if (!next_idx_to_u.count(idx)) {
        next_idx_to_u[idx] = (int)next_level_fu->size();
        next_level_fu->push_back(ur_f);
      }
      H_level.insert(i, next_idx_to_u[idx]) = 0.125;
      
      
      idx = tree.get_face_idx(lr_f);
      if (!next_idx_to_u.count(idx)) {
        next_idx_to_u[idx] = (int)next_level_fu->size();
        next_level_fu->push_back(lr_f);
      }
      H_level.insert(i, next_idx_to_u[idx]) = 0.125;
    }
    
    cur_level_fu = next_level_fu;
    int n_next = (int)next_level_fu->size();
    H_level.conservativeResize(n_cur, n_next);
    n_cur = n_next;
    
    H_u = H_u * H_level;
    
    if (k == tree.max_depth - 2) {
      u_reg_faces = *next_level_fu;
    }
  }
  
  n_u_reg = (int)u_reg_faces.size();
  reg_u.resize(n_u_reg);
  
  for (int i = 0; i < n_u_reg; ++i) {
    Face& f = u_reg_faces[i];
    reg_u[i] = u(f.i + 1, f.j);
  }
  
  tree_u = H_u * reg_u;
  
  
  // Set up the matrix H_v transforming v's on regular grids to tree grids.
  n_v_tree = (int)v_tree_faces.size();
  H_v.resize(n_v_tree, n_v_tree);
  H_v.setIdentity();
  
  v_reg_faces = v_tree_faces;
  
  n_cur = n_v_tree;
  vector<Face>* cur_level_fv = &v_tree_faces;
  vector<Face>* next_level_fv;
  
  for (int k = 0; k < tree.max_depth-1; ++k) {
    next_level_fv = new vector<Face>();
    unordered_map<int, int> next_idx_to_v;
    SpMat H_level(n_cur, ni * nj);
    
    for (int i = 0; i < n_cur; ++i) {
      Face& f = (*cur_level_fv)[i];
      
      int idx;
      if (f.depth != k) {
        idx = tree.get_face_idx(f);
        if (!next_idx_to_v.count(idx)) {
          next_idx_to_v[idx] = (int)next_level_fv->size();
          next_level_fv->push_back(f);
        }
        H_level.insert(i, next_idx_to_v[idx]) = 1;
        continue;
      }
      
      Cell c(f.depth, f.i, f.j);
      
      Cell left_c = tree.get_child(c, NW);
      Cell right_c = tree.get_child(c, NE);
      
      Face ul_f(left_c.depth, left_c.i, left_c.j+1, TOP);
      Face ur_f(right_c.depth, right_c.i, right_c.j+1, TOP);
      
      Face ml_f(left_c.depth, left_c.i, left_c.j, TOP);
      Face mr_f(right_c.depth, right_c.i, right_c.j, TOP);
      
      Face ll_f(left_c.depth, left_c.i, left_c.j-1, TOP);
      Face lr_f(right_c.depth, right_c.i, right_c.j-1, TOP);
      
      idx = tree.get_face_idx(ul_f);
      if (!next_idx_to_v.count(idx)) {
        next_idx_to_v[idx] = (int)next_level_fv->size();
        next_level_fv->push_back(ul_f);
      }
      H_level.insert(i, next_idx_to_v[idx]) = 0.125;
      
      
      idx = tree.get_face_idx(ur_f);
      if (!next_idx_to_v.count(idx)) {
        next_idx_to_v[idx] = (int)next_level_fv->size();
        next_level_fv->push_back(ur_f);
      }
      H_level.insert(i, next_idx_to_v[idx]) = 0.125;
      
      
      idx = tree.get_face_idx(ml_f);
      if (!next_idx_to_v.count(idx)) {
        next_idx_to_v[idx] = (int)next_level_fv->size();
        next_level_fv->push_back(ml_f);
      }
      H_level.insert(i, next_idx_to_v[idx]) = 0.25;
      
      
      idx = tree.get_face_idx(mr_f);
      if (!next_idx_to_v.count(idx)) {
        next_idx_to_v[idx] = (int)next_level_fv->size();
        next_level_fv->push_back(mr_f);
      }
      H_level.insert(i, next_idx_to_v[idx]) = 0.25;
      
      
      idx = tree.get_face_idx(ll_f);
      if (!next_idx_to_v.count(idx)) {
        next_idx_to_v[idx] = (int)next_level_fv->size();
        next_level_fv->push_back(ll_f);
      }
      H_level.insert(i, next_idx_to_v[idx]) = 0.125;
      
      
      idx = tree.get_face_idx(lr_f);
      if (!next_idx_to_v.count(idx)) {
        next_idx_to_v[idx] = (int)next_level_fv->size();
        next_level_fv->push_back(lr_f);
      }
      H_level.insert(i, next_idx_to_v[idx]) = 0.125;
    }
    
    cur_level_fv = next_level_fv;
    int n_next = (int)next_level_fv->size();
    H_level.conservativeResize(n_cur, n_next);
    n_cur = n_next;
    
    H_v = H_v * H_level;
    
    if (k == tree.max_depth - 2) {
      v_reg_faces = *next_level_fv;
    }
  }
  
  n_v_reg = (int)v_reg_faces.size();
  reg_v.resize(n_v_reg);
  
  for (int i = 0; i < n_v_reg; ++i) {
    Face& f = v_reg_faces[i];
    reg_v[i] = v(f.i, f.j + 1);
  }
  tree_v = H_v * reg_v;
  
  
  reg_uv.resize(n_u_reg + n_v_reg);
  reg_uv.head(n_u_reg) = reg_u;
  reg_uv.tail(n_v_reg) = reg_v;
  
  M_reg_uv.resize(n_u_reg + n_v_reg);
  
  for (int i = 0; i < n_u_reg; ++i) {
    Face& u_f = u_reg_faces[i];
    M_reg_uv[i] = u_vol(u_f.i + 1, u_f.j) * sqr(dx);
  }
  
  for (int i = 0; i < n_v_reg; ++i) {
    Face& v_f = v_reg_faces[i];
    M_reg_uv[i + n_u_reg] = v_vol(v_f.i, v_f.j + 1) * sqr(dx);
  }
  
  H = build_large_diagnol(H_u, H_v);
}


// Get tau's at cell centres, e.g. tau11 and tau22.
void VisSolver::get_grid_taus() {
  for (int i = 0; i < tree.leaf_cell_count; ++i) {
    Cell& c = tree.leaf_cells[i];
    if (c.depth != tree.max_depth - 1) {
      grid_to_tau11[i] = (int)tau11_to_cell.size();
      tau11_to_cell.push_back(i);
      
      grid_to_tau22[i] = (int) tau22_to_cell.size();
      tau22_to_cell.push_back(i);
      continue;
    }
    
    // The order of surrounding faces is (right, left, top, bottom).
    int right_f_idx = tree.cell_to_face_map(i, 0);
    int left_f_idx = tree.cell_to_face_map(i, 1);
    
    if (c_vol(c.i, c.j) > EPSILON && (face_to_u.count(right_f_idx)
                                      || face_to_u.count(left_f_idx))) {
      grid_to_tau11[i] = (int)tau11_to_cell.size();
      tau11_to_cell.push_back(i);
    }
    
    int top_f_idx = tree.cell_to_face_map(i, 2);
    int bottom_f_idx = tree.cell_to_face_map(i, 3);
    
    if (c_vol(c.i, c.j) > EPSILON && (face_to_v.count(top_f_idx)
                                      || face_to_v.count(bottom_f_idx))) {
      grid_to_tau22[i] = (int)tau22_to_cell.size();
      tau22_to_cell.push_back(i);
    }
  }
}


// Get tau's at grid nodes, e.g. tau12.
void VisSolver::get_node_taus() {
  for (int i = 0; i < tree.node_count; ++i) {
    Node& n = tree.nodes[i];
    
    int n_d = tree.get_level_dims(n.depth);
    if ((n.i == 0 && n.position == SW)
        || (n.i == n_d - 1 && n.position == NE)
        || (n.j == 0 && n.position == SW)
        || (n.j == n_d - 1 && n.position == NE)) {
      continue;
    }
    
    if (n.depth != tree.max_depth - 1 || n.position != NE) {
      node_to_tau12[i] = (int)tau12_to_node.size();
      tau12_to_node.push_back(i);
      continue;
    }
    
    std::vector<int>& neighbor_faces = tree.node_to_face_map[i];
    
    int right_f_idx = neighbor_faces[0];
    int left_f_idx = neighbor_faces[1];
    int top_f_idx = neighbor_faces[2];
    int bottom_f_idx = neighbor_faces[3];
    
    if (n_vol(n.i + 1, n.j + 1) > EPSILON && (face_to_u.count(top_f_idx)
                                              || face_to_u.count(bottom_f_idx) || face_to_v.count(right_f_idx)
                                              || face_to_v.count(left_f_idx))) {
      node_to_tau12[i] = (int)tau12_to_node.size();
      tau12_to_node.push_back(i);
    }
  }
}

void VisSolver::compute_deformation_operator() {
  // Set matrix for the velocity of t11.
  for (int i = 0; i < n_t11; ++i) {
    
    int idx = i;
    int c_idx = tau11_to_cell[i];
    Cell& c = tree.leaf_cells[c_idx];
    
    // The order of belonging faces are (right, left, top, bottom).
    int right_f_idx = tree.cell_to_face_map(c_idx, 0);
    int left_f_idx = tree.cell_to_face_map(c_idx, 1);
    
    float dx = tree.get_cell_width(c.depth);
    
    if (face_to_u.count(right_f_idx)) {
      D.insert(idx, face_to_u[right_f_idx]) = -2 / dx;
    }
    if (face_to_u.count(left_f_idx)) {
      D.insert(idx, face_to_u[left_f_idx]) = 2 / dx;
    }
  }
  
  // Set matrix for the velocity of t22.
  for (int i = 0; i < n_t22; ++i) {
    
    int idx = n_t11 + i;
    
    int c_idx = tau22_to_cell[i];
    Cell& c = tree.leaf_cells[c_idx];
    
    // The order of belonging faces are (right, left, top, bottom).
    int top_f_idx = tree.cell_to_face_map(c_idx, 2);
    int bottom_f_idx = tree.cell_to_face_map(c_idx, 3);
    
    float dx = tree.get_cell_width(c.depth);
    
    if (face_to_v.count(top_f_idx)) {
      D.insert(idx, n_u + face_to_v[top_f_idx]) = -2 / dx;
    }
    if (face_to_v.count(bottom_f_idx)) {
      D.insert(idx, n_u + face_to_v[bottom_f_idx]) = 2 / dx;
    }
  }
  
  // Set matrix for the veloc ity of t12.
  for (int i = 0; i < n_t12; ++i) {
    int idx = n_t11 + n_t22 + i;
    
    int n_idx = tau12_to_node[i];
    Node& n = tree.nodes[n_idx];
    
    // The order of faces are (right, left, top, bottom) respective to the
    // node.
    std::vector<int>& neighbor_faces = tree.node_to_face_map[n_idx];
    
    int right_f_idx = neighbor_faces[0];
    int left_f_idx = neighbor_faces[1];
    int top_f_idx = neighbor_faces[2];
    int bottom_f_idx = neighbor_faces[3];
    
    // If there is a T junction and above cell is larger.
    if (top_f_idx == -1) {
      // Right face to the node, which is also the large face.
      int large_f_idx = right_f_idx;
      // Large cell at the T junction.
      int large_c_idx = tree.face_to_cell_map[large_f_idx].first[0];
      Cell& large_c = tree.leaf_cells[large_c_idx];
      
      // The node's depth is equal to the smallest cell's depth.
      float dx_u = 0.5 * (tree.get_cell_width(n.depth) +
                          tree.get_cell_width(large_c.depth));
      int top_right_u_f_idx = tree.cell_to_face_map(large_c_idx, 0);
      int top_left_u_f_idx = tree.cell_to_face_map(large_c_idx, 1);
      
      D.insert(idx, face_to_u[top_right_u_f_idx]) = 0.5 * (-1 / dx_u);
      
      D.insert(idx, face_to_u[top_left_u_f_idx]) = 0.5 * (-1 / dx_u);
      
      if (face_to_u.count(bottom_f_idx)) {
        D.insert(idx, face_to_u[bottom_f_idx]) = 1 / dx_u;
      }
      
      Cell left_to_large_c = tree.left(large_c);
      Cell right_to_large_c = tree.right(large_c);
      
      float dx_left = 0.5 * tree.get_cell_width(large_c.depth);
      float dx_right = dx_left;
      
      int left_v_f_idx = -1, right_v_f_idx = -1;
      
      Cell left_leaf_c = tree.get_leaf_cell(left_to_large_c, SE);
      int left_leaf_c_idx = tree.get_cell_idx(left_leaf_c);
      dx_left += 0.5 * tree.get_cell_width(left_leaf_c.depth);
      left_v_f_idx = tree.cell_to_face_map(left_leaf_c_idx, 3);
      
      Cell right_leaf_c = tree.get_leaf_cell(right_to_large_c, SW);
      int right_leaf_c_idx = tree.get_cell_idx(right_leaf_c);
      dx_right += 0.5 * tree.get_cell_width(right_leaf_c.depth);
      right_v_f_idx = tree.cell_to_face_map(right_leaf_c_idx, 3);
      
      float dx_v = dx_left + dx_right;
      
      if (face_to_v.count(left_f_idx)) {
        D.insert(idx, n_u + face_to_v[left_v_f_idx]) = 1 / dx_v;
      }
      
      if (face_to_v.count(right_v_f_idx)) {
        D.insert(idx, n_u + face_to_v[right_v_f_idx]) = -1 / dx_v;
      }
      // If there is a T junction and below cell is larger.
    } else if (bottom_f_idx == -1) {
      // Right face to the node.
      int large_f_idx = right_f_idx;
      int large_c_idx = tree.face_to_cell_map[large_f_idx].second[0];
      Cell& large_c = tree.leaf_cells[large_c_idx];
      
      float dx_u = 0.5*(tree.get_cell_width(n.depth) +
                        tree.get_cell_width(large_c.depth));
      int bottom_right_u_f_idx = tree.cell_to_face_map(large_c_idx, 0);
      int bottom_left_u_f_idx = tree.cell_to_face_map(large_c_idx, 1);
      
      D.insert(idx, face_to_u[bottom_right_u_f_idx]) = 0.5 * (1 / dx_u);
      
      D.insert(idx, face_to_u[bottom_left_u_f_idx]) = 0.5 * (1 / dx_u);
      
      if (face_to_u.count(top_f_idx)) {
        D.insert(idx, face_to_u[top_f_idx]) = -1 / dx_u;
      }
      
      Cell left_to_large_c = tree.left(large_c);
      Cell right_to_large_c = tree.right(large_c);
      float dx_left = 0.5 * tree.get_cell_width(large_c.depth);
      float dx_right = dx_left;
      int left_v_f_idx = -1, right_v_f_idx = -1;
      
      Cell left_leaf_c = tree.get_leaf_cell(left_to_large_c, NE);
      int left_leaf_c_idx = tree.get_cell_idx(left_leaf_c);
      dx_left += 0.5 * tree.get_cell_width(left_leaf_c.depth);
      left_v_f_idx = tree.cell_to_face_map(left_leaf_c_idx, 2);
      
      Cell right_leaf_c = tree.get_leaf_cell(right_to_large_c, NW);
      int right_leaf_c_idx = tree.get_cell_idx(right_leaf_c);
      dx_right += 0.5 * tree.get_cell_width(right_leaf_c.depth);
      right_v_f_idx = tree.cell_to_face_map(right_leaf_c_idx, 2);
      
      float dx_v = dx_left + dx_right;
      
      if (face_to_v.count(right_v_f_idx)) {
        D.insert(idx, n_u + face_to_v[right_v_f_idx]) = -1 / dx_v;
      }
      
      if (face_to_v.count(left_v_f_idx)) {
        D.insert(idx, n_u + face_to_v[left_v_f_idx]) = 1 / dx_v;
      }
      // If there is a T junction and left cell is larger.
    } else if (left_f_idx == -1) {
      // Top face to the node, which is also the large face.
      int large_f_idx = top_f_idx;
      int large_c_idx = tree.face_to_cell_map[large_f_idx].second[0];
      Cell& large_c = tree.leaf_cells[large_c_idx];
      
      // The node's depth is equal to the smallest cell's depth.
      float dx_v = 0.5*(tree.get_cell_width(n.depth) +
                        tree.get_cell_width(large_c.depth));
      int left_top_v_f_idx = tree.cell_to_face_map(large_c_idx, 2);
      int left_bottom_v_f_idx = tree.cell_to_face_map(large_c_idx, 3);
      
      D.insert(idx, n_u + face_to_v[left_top_v_f_idx]) = 0.5 * (1 / dx_v);
      
      D.insert(idx, n_u + face_to_v[left_bottom_v_f_idx]) = 0.5 * (1 / dx_v);
      
      if (face_to_v.count(right_f_idx)) {
        D.insert(idx, n_u + face_to_v[right_f_idx]) = -1 / dx_v;
      }
      
      Cell above_to_large_c = tree.above(large_c);
      Cell below_to_large_c = tree.below(large_c);
      float dx_above = 0.5 * tree.get_cell_width(large_c.depth);
      float dx_below = dx_above;
      int bottom_u_f_idx = -1, top_u_f_idx = -1;
      
      Cell above_leaf_c = tree.get_leaf_cell(above_to_large_c, SE);
      int above_leaf_c_idx = tree.get_cell_idx(above_leaf_c);
      dx_above += 0.5 * tree.get_cell_width(above_leaf_c.depth);
      
      top_u_f_idx = tree.cell_to_face_map(above_leaf_c_idx, 0);
      
      Cell below_leaf_c = tree.get_leaf_cell(below_to_large_c, NE);
      int below_leaf_cell_idx = tree.get_cell_idx(below_leaf_c);
      dx_below += 0.5 * tree.get_cell_width(below_leaf_c.depth);
      bottom_u_f_idx = tree.cell_to_face_map(below_leaf_cell_idx, 0);
      
      float dx_u = dx_above + dx_below;
      
      if (face_to_u.count(bottom_f_idx)) {
        D.insert(idx, face_to_u[bottom_u_f_idx]) = 1 / dx_u;
      }
      
      if (face_to_u.count(top_u_f_idx)) {
        D.insert(idx, face_to_u[top_u_f_idx]) = -1 / dx_u;
      }
      // If there is a T junction and right cell is larger.
    } else if (right_f_idx == -1) {
      // Top face to the node, which is also the large face.
      int large_f_idx = neighbor_faces[2];
      int large_c_idx = tree.face_to_cell_map[large_f_idx].first[0];
      Cell& large_c = tree.leaf_cells[large_c_idx];
      float dx_v = 0.5*(tree.get_cell_width(n.depth) +
                        tree.get_cell_width(large_c.depth));
      int right_top_v_f_idx = tree.cell_to_face_map(large_c_idx, 2);
      int right_bottom_v_f_idx = tree.cell_to_face_map(large_c_idx, 3);
      
      D.insert(idx, n_u + face_to_v[right_top_v_f_idx]) = 0.5 * (-1 / dx_v);
      
      D.insert(idx, n_u + face_to_v[right_bottom_v_f_idx]) = 0.5 * (-1 / dx_v);
      
      if (face_to_v.count(left_f_idx)) {
        D.insert(idx, n_u + face_to_v[left_f_idx]) = 1 / dx_v;
      }
      
      Cell above_to_large_c = tree.above(large_c);
      Cell below_to_large_c = tree.below(large_c);
      float dx_above = 0.5 * tree.get_cell_width(large_c.depth);
      float dx_below = dx_above;
      int bottom_u_f_idx = -1, top_u_f_idx = -1;
      
      Cell above_leaf_c = tree.get_leaf_cell(above_to_large_c, SW);
      int above_leaf_c_idx = tree.get_cell_idx(above_leaf_c);
      dx_above += 0.5 * tree.get_cell_width(above_leaf_c.depth);
      top_u_f_idx = tree.cell_to_face_map(above_leaf_c_idx, 1);
      
      
      Cell below_leaf_c = tree.get_leaf_cell(below_to_large_c, NW);
      int below_leaf_c_idx = tree.get_cell_idx(below_leaf_c);
      dx_below += 0.5 * tree.get_cell_width(below_leaf_c.depth);
      bottom_u_f_idx = tree.cell_to_face_map(below_leaf_c_idx, 1);
      
      float dx_u = dx_above + dx_below;
      
      if (face_to_u.count(bottom_u_f_idx)) {
        D.insert(idx, face_to_u[bottom_u_f_idx]) = 1 / dx_u;
      }
      
      if (face_to_u.count(top_u_f_idx)) {
        D.insert(idx, face_to_u[top_u_f_idx]) = -1 / dx_u;
      }
    } else {
      
      Face& right_f = tree.qt_faces[right_f_idx];
      Face& left_f = tree.qt_faces[left_f_idx];
      
      // dx is equal to the distance between two face centers.
      float dx_v = 0.5*(tree.get_cell_width(right_f.depth) +
                        tree.get_cell_width(left_f.depth));
      
      if (face_to_v.count(right_f_idx)) {
        D.insert(idx, n_u + face_to_v[right_f_idx]) = -1 / dx_v;
      }
      
      if (face_to_v.count(left_f_idx)) {
        D.insert(idx, n_u + face_to_v[left_f_idx]) = 1 / dx_v;
      }
      
      Face& top_f = tree.qt_faces[top_f_idx];
      Face& bottom_f = tree.qt_faces[bottom_f_idx];
      
      // dx is equal to the distance between two face centers.
      float dx_u = 0.5*(tree.get_cell_width(top_f.depth) +
                        tree.get_cell_width(bottom_f.depth));
      
      if (face_to_u.count(top_f_idx)) {
        D.insert(idx, face_to_u[top_f_idx]) = -1 / dx_u;
      }
      
      if (face_to_u.count(bottom_f_idx)) {
        D.insert(idx, face_to_u[bottom_f_idx]) = 1 / dx_u;
      }
      
    }
  }
}


void VisSolver::compute_factor_matrix() {
  
  for (int i = 0; i < n_t11 + n_t22; ++i) {
    factor.insert(i, i) = 1;
  }
  
  for (int i = 0; i < n_t12; ++i) {
    int idx = n_t11 + n_t22 + i;
    factor.insert(idx, idx) = 2;
  }
}


void VisSolver::compute_volume_matrices() {
  
  // Compute control volumes of u's and v's.
  for (int i = 0; i < n_u; ++i) {
    Face& f = tree.qt_faces[u_to_face[i]];
    if (f.depth == tree.max_depth - 1) {
      M_uv.insert(i, i) = u_vol(f.i + 1, f.j) * sqr(dx);
    } else {
      int f_idx = u_to_face[i];
      float area = sqr(tree.get_cell_width(f.depth));
      if (tree.face_to_cell_map[f_idx].first.size() != 1
          || tree.face_to_cell_map[f_idx].second.size() != 1) {
        // it's only valid for a graded tree.
        area *= 0.75;
      }
      M_uv.insert(i, i) = area;
    }
  }
  
  for (int i = 0; i < n_v; ++i) {
    Face& f = tree.qt_faces[v_to_face[i]];
    if (f.depth == tree.max_depth - 1) {
      M_uv.insert(n_u + i, n_u + i) = v_vol(f.i, f.j + 1) * sqr(dx);
    } else {
      int f_idx = v_to_face[i];
      float area = sqr(tree.get_cell_width(f.depth));
      if (tree.face_to_cell_map[f_idx].first.size() != 1
          || tree.face_to_cell_map[f_idx].second.size() != 1) {
        // it's only valid for a graded tree.
        area *= 0.75;
      }
      M_uv.insert(n_u + i, n_u + i) = area;
    }
  }
  
  // Compute control volumes of taus.
  for (int i = 0; i < n_t11; ++i) {
    Cell& c = tree.leaf_cells[tau11_to_cell[i]];
    if (c.depth == tree.max_depth - 1) {
      M_t.insert(i, i) = c_vol(c.i, c.j) * sqr(dx);
    } else {
      M_t.insert(i, i) = sqr(tree.get_cell_width(c.depth));
    }
  }
  
  for (int i = 0; i < n_t22; ++i) {
    Cell& c = tree.leaf_cells[tau22_to_cell[i]];
    if (c.depth == tree.max_depth - 1) {
      M_t.insert(n_t11 + i, n_t11 + i) = c_vol(c.i, c.j) * sqr(dx);
    } else {
      M_t.insert(n_t11 + i, n_t11 + i) = sqr(tree.get_cell_width(c.depth));
    }
  }
  
  for (int i = 0; i < n_t12; ++i) {
    int node_idx = tau12_to_node[i];
    Node& n = tree.nodes[node_idx];
    float area = 0;
    std::vector<int> neighbor_faces = tree.node_to_face_map[node_idx];
    
    int right_f_idx = neighbor_faces[0];
    int left_f_idx = neighbor_faces[1];
    int top_f_idx = neighbor_faces[2];
    int bottom_f_idx = neighbor_faces[3];
    
    float dx1 = 0;
    float dx2 = 0;
    
    if (n.depth == tree.max_depth - 1 && right_f_idx != -1
        && left_f_idx != -1 && top_f_idx != -1 && bottom_f_idx != -1) {
      M_t.insert(n_t11 + n_t22 + i, n_t11 + n_t22 + i)
      = n_vol(n.i + 1, n.j + 1) *  sqr(dx);
      continue;
    }
    
    if (right_f_idx == -1 || left_f_idx == -1) {
      int large_f_idx = neighbor_faces[2];
      int large_c_idx = tree.face_to_cell_map[large_f_idx].first[0];
      if (left_f_idx == -1) {
        large_c_idx = tree.face_to_cell_map[large_f_idx].second[0];
      }
      Cell& large_c = tree.leaf_cells[large_c_idx];
      dx1 = tree.get_cell_width(n.depth);
      dx2 = tree.get_cell_width(large_c.depth);
      area = min(dx1, dx2) * (dx1 + dx2) * 0.5;
      
    } else if (top_f_idx == -1 || bottom_f_idx == -1) {
      int large_f_idx = neighbor_faces[0];
      int large_c_idx = tree.face_to_cell_map[large_f_idx].first[0];
      if (bottom_f_idx == -1) {
        large_c_idx = tree.face_to_cell_map[large_f_idx].second[0];
      }
      Cell& large_c = tree.leaf_cells[large_c_idx];
      dx1 = tree.get_cell_width(n.depth);
      dx2 = tree.get_cell_width(large_c.depth);
      area = min(dx1, dx2) * (dx1 + dx2) * 0.5;
      
    } else {
      Cell top_right_c = tree.leaf_cells[tree.face_to_cell_map[top_f_idx].first[0]];
      if (tree.face_to_cell_map[top_f_idx].first.size() == 2) {
        top_right_c = tree.leaf_cells[tree.face_to_cell_map[top_f_idx].first[1]];
      }
      float tr_dx = tree.get_cell_width(top_right_c.depth);
      
      Cell top_left_c = tree.leaf_cells[tree.face_to_cell_map[top_f_idx].second[0]];
      if (tree.face_to_cell_map[top_f_idx].second.size() == 2) {
        top_left_c = tree.leaf_cells[tree.face_to_cell_map[top_f_idx].second[1]];
      }
      float tl_dx = tree.get_cell_width(top_left_c.depth);
      
      Cell& bottom_right_c = tree.leaf_cells[tree.face_to_cell_map[bottom_f_idx].first[0]];
      float br_dx = tree.get_cell_width(bottom_right_c.depth);
      
      Cell& bottom_left_c = tree.leaf_cells[tree.face_to_cell_map[bottom_f_idx].second[0]];
      float bl_dx = tree.get_cell_width(bottom_left_c.depth);
      
      float min_dx = min(tr_dx, tl_dx, br_dx, bl_dx);
      int num_small_cells = 0;
      if (tr_dx == min_dx) {
        num_small_cells++;
      }
      if (tl_dx == min_dx) {
        num_small_cells++;
      }
      if (br_dx == min_dx) {
        num_small_cells++;
      }
      if (bl_dx == min_dx) {
        num_small_cells++;
      }
      if (num_small_cells == 1) {
        area = sqr(min_dx) * 9 / 4;
      } else if (num_small_cells == 2) {
        if (tl_dx == br_dx || tr_dx == bl_dx) {
          area = sqr(min_dx);
        } else {
          area = sqr(min_dx) * 6 / 4;
        }
      } else {
        area = sqr(min_dx);
      }
    }
    M_t.insert(n_t11 + n_t22 + i, n_t11 + n_t22 + i) = area;
  }
}

void VisSolver::compute_tau_vis() {
  for (int i = 0; i < n_t11; ++i) {
    Cell& c = tree.leaf_cells[tau11_to_cell[i]];
    M_vis.insert(i, i) = viscosity(c.i, c.j);
  }
  
  for (int i = 0; i < n_t22; ++i) {
    Cell& c = tree.leaf_cells[tau22_to_cell[i]];
    M_vis.insert(n_t11 + i, n_t11 + i) = viscosity(c.i, c.j);
  }
  
  for (int i = 0; i < n_t12; ++i) {
    Node& n = tree.nodes[tau12_to_node[i]];
    // Suppose the node is at the NE corner of the cell.
    M_vis.insert(n_t11 + n_t22 + i, n_t11 + n_t22 + i) = 0.25
    * (viscosity(n.i, n.j) + viscosity(n.i + 1, n.j)
       + viscosity(n.i, n.j + 1) + viscosity(n.i + 1, n.j + 1));
  }
}

void VisSolver::compute_rhs() {
  for (int i = 0; i < n_u; ++i) {
    rhs[i] = M_uv.coeff(i, i) * tree_u[i];
  }
  for (int i = 0; i < n_v; ++i) {
    int idx = n_u + i;
    rhs[idx] = M_uv.coeff(idx, idx) * tree_v[i];
  }
}


// Getter functions for u and v faces on the quadtree.
vector<Face> VisSolver::get_u_tree_faces() {
  return u_tree_faces;
}


vector<Face> VisSolver::get_v_tree_faces() {
  return v_tree_faces;
}


// Divide each control volume into several subdivision and count the actual
// volume on the free surface.
void compute_volume_fractions(const Array2f& levelset, Array2f& fractions,
                              Vec2f fraction_origin, int subdivision) {
  
  //Assumes levelset and fractions have the same dx
  float sub_dx = 1.0 / subdivision;
  int sample_max = subdivision*subdivision;
  for (int j = 0; j < fractions.nj; ++j) {
    for (int i = 0; i < fractions.ni; ++i) {
      float start_x = fraction_origin[0] + (float)i;
      float start_y = fraction_origin[1] + (float)j;
      int incount = 0;
      
      for (int sub_j = 0; sub_j < subdivision; ++sub_j) {
        for (int sub_i = 0; sub_i < subdivision; ++sub_i) {
          float x_pos = start_x + (sub_i + 0.5)*sub_dx;
          float y_pos = start_y + (sub_j + 0.5)*sub_dx;
          float phi_val = interpolate_value(Vec2f(x_pos, y_pos), levelset);
          if (phi_val < 0)
            ++incount;
        }
      }
      fractions(i, j) = (float)incount / (float)sample_max;
    }
  }
}



SpMat build_large_diagnol(SpMat m1, SpMat m2) {
  SpMat m(m1.rows() + m2.rows(), m1.cols() + m2.cols());
  m.reserve(m1.nonZeros() + m2.nonZeros());
  
  for (int k = 0; k < m1.outerSize(); ++k) {
    for (SpMat::InnerIterator it(m1,k); it; ++it) {
      m.insert(it.row(), it.col()) = it.value();
    }
  }
  
  int n_rows = (int)m1.rows(), n_cols = (int)m1.cols();
  
  for (int k = 0; k < m2.outerSize(); ++k) {
    for (SpMat::InnerIterator it(m2,k); it; ++it) {
      m.insert(it.row() + n_rows, it.col() + n_cols) = it.value();
    }
  }
  return m;
}

SparseMatrixd copy_sparse_matrix(SpMat m) {
  SparseMatrixd m1((int)m.rows());
  m1.zero();
  for (int k = 0; k < m.outerSize(); ++k) {
    for (SpMat::InnerIterator it(m,k); it; ++it) {
      m1.set_element((int)it.row(), (int)it.col(), it.value());
    }
  }
  return m1;
}


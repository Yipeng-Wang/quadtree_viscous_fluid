#ifndef FLUIDSIM_H
#define FLUIDSIM_H

#include "array2.h"
#include "vec.h"
#include "pcgsolver/sparse_matrix.h"
#include "pcgsolver/pcg_solver.h"
#include "viscositysolver.h"

#include <vector>

static Vec2i cell_offset[] = { Vec2i(-1,0), Vec2i(1,0), Vec2i(0,-1), Vec2i(0,1) }; // cell to cell offset

enum marked { UNVISITED, VISITED, FINISHED }; // BFS marker type.

typedef Array2<marked, Array1<marked>> Array2m; // define the type for the redistancing marker.

class FluidSim {
  
public:
  void initialize(float width, int ni_, int nj_);
  void set_boundary(float(*phi)(const Vec2f&));
  void advance(float dt);
  
  // Grid dimensions
  int ni, nj;
  float width;
  float dx;
  
  // Fluid velocity
  Array2f u, v;
  Array2f temp_u, temp_v;
  
  // Static geometry representation
  Array2f nodal_solid_phi;
  
  // Data for pressure solve and extrapolation
  Array2c u_valid, v_valid;
  Array2f liquid_phi; //extracted from particles
  Array2f u_weights, v_weights;
  
  // Data for viscosity solve, c_vol is cell volumes and n_vol is node volumes.
  Array2f u_vol, v_vol, c_vol, n_vol;
  Array2f viscosity;
  
  std::vector<Vec2f> particles; //For marker particle simulation, it stores the coordinates of particles.
  float particle_radius;
  
  // Data arrays for extrapolation
  Array2c valid, old_valid;
  
  // Data for redistancing.
  float redis_depth;
  
  //Solver data
  PCGSolver<double> solver;
  SparseMatrixd matrix;
  std::vector<double> rhs;
  std::vector<double> pressure;
  
  SparseMatrixd vmatrix;
  std::vector<double> vrhs;
  std::vector<double> velocities;
  
  Vec2f get_velocity(const Vec2f& position);
  void add_particle(const Vec2f& position);
  
  static bool solve_on_quadtree;
  
private:
  
  Vec2f trace_rk2(const Vec2f& position, float dt);
  
  void advect_particles(float dt);
  
  void compute_phi();
  // Redistance to meaningful signed distance.
  void redistance(size_t);
  Vec2f interface_search (Vec2f grid_pos, size_t iter_limit);
  void fast_marching(Array2m marked_cells);
  
  float cfl();
  
  //fluid velocity operations
  void advect(float dt);
  void add_force(float dt);
  
  void apply_projection(float dt);
  void compute_pressure_weights();
  void solve_pressure(float dt);
  
  int u_ind(int i, int j);
  int v_ind(int i, int j);
  void apply_viscosity(float dt);
  
  void apply_viscosity_quadtree(float dt);
  
  void compute_viscosity_weights();
  void solve_viscosity(float dt);
  
  void constrain_velocity();
};

#endif

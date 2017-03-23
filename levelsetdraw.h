//
//  levelsetdraw.h
//  SurfaceFluid
//  2d level set plot, create surface edges and vertices in two vectors.
//  use draw_segmentset2d in openglutils.h to plot the surface directly.
//
//  Created by Yipeng Wang on 2017-03-21.
//

#ifndef levelsetdraw_h
#define levelsetdraw_h

#include "array2.h"
#include "vec.h"

static const int mc_template[16][4] =
{ { -1, -1, -1, -1 },
    { 3, 0, -1, -1 },
    { 0, 1, -1, -1 },
    { 3, 1, -1, -1 },
    
    { 1, 2, -1, -1 },
    { 3, 0, 1, 2 },
    { 0, 2, -1, -1 },
    { 3, 2, -1, -1 },
    
    { 2, 3, -1, -1 },
    { 2, 0, -1, -1 },
    { 0, 1, 2, 3 },
    { 2, 1, -1, -1 },
    
    { 1, 3, -1, -1 },
    { 1, 0, -1, -1 },
    { 0, 3, -1, -1 },
    { -1, -1, -1, -1 } };

static const Vec2i node_in_cell[4] =
{ Vec2i(0,0), Vec2i(1,0), Vec2i(1,1), Vec2i(0,1) };

static const Vec2i edge_node_map[4] =
{ Vec2i(0,1), Vec2i(1,2), Vec2i(2,3), Vec2i(3,0) };

class LevelSetDraw {
    
public:
    mutable std::vector<Vec2f> verts;
    mutable std::vector<Vec2i> edges;
    
    // FluidQuadTree constructor.
    LevelSetDraw(Array2f liquid_phi_, float dx_) : liquid_phi(liquid_phi_), dx(dx_){};
    
    void extract_mesh() const {
        int ni = liquid_phi.ni, nj = liquid_phi.nj;
        
        // Run marching squares loop
        for (int x = 0; x < ni - 1; ++x)
            for (int y = 0; y < nj - 1; ++y) {
                Vec2i idx(x, y);
                int mc_idx = 0;
                
                for (int i = 0; i < 4; ++i) {
                    Vec2i coord = Vec2i(x, y) + node_in_cell[i];
                    if (liquid_phi(coord[0], coord[1]) <= 0.) {
                        mc_idx += (1 << i);
                    }
                }
                
                // Connect edges using the marching squares template
                for (int e = 0; e < 4 && mc_template[mc_idx][e] >= 0; e += 2) {
                    // Find first vertex
                    int edge = mc_template[mc_idx][e];
                    Vec2i n0_idx = idx + node_in_cell[edge_node_map[edge][0]];
                    Vec2i n1_idx = idx + node_in_cell[edge_node_map[edge][1]];
                    
                    Vec2f v0 = interp_interface(n0_idx, n1_idx);
                    
                    // Find second vertex
                    edge = mc_template[mc_idx][e + 1];
                    n0_idx = idx + node_in_cell[edge_node_map[edge][0]];
                    n1_idx = idx + node_in_cell[edge_node_map[edge][1]];
                    
                    Vec2f v1 = interp_interface(n0_idx, n1_idx);
                    
                    // Store vertices
                    Vec2f v0_ws((v0[0] + 0.5) * dx, (v0[1] + 0.5) * dx);
                    Vec2f v1_ws((v1[0] + 0.5) * dx, (v1[1] + 0.5) * dx);
                    
                    verts.push_back(v0_ws);
                    verts.push_back(v1_ws);
                    
                    edges.push_back(Vec2i((int)verts.size() - 2, (int)verts.size() - 1));
                }
            }
    };
    
    Vec2f interp_interface(const Vec2i& i0, const Vec2i& i1) const {
        assert(liquid_phi(i0[0], i0[1]) * liquid_phi(i1[0], i1[1]) <= 0.0);
        //Find weight to zero isosurface
        double s = liquid_phi(i0[0], i0[1]) / (liquid_phi(i0[0], i0[1]) - liquid_phi(i1[0], i1[1]));
        
        if (s < 0.0) s = 0.0;
        else if (s > 1.0) s = 1.0;
        
        Vec2f dis = Vec2f(i1) - Vec2f(i0);
        return Vec2f(i0[0] + dis[0]*s, i0[1] + dis[1]*s);
    };
    
private:
    Array2f liquid_phi;
    float dx;
};

#endif

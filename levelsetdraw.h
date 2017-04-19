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
    
    void extract_mesh() const;
    
private:
    Array2f liquid_phi;
    float dx;
    
    Vec2f interp_interface(const Vec2i& i0, const Vec2i& i1) const;
};

#endif

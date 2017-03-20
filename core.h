//
//  core.h
//  SurfaceFluid
//
//  Created by Yipeng Wang on 2017-03-18.
//  Copyright © 2017 Yipeng Wang. All rights reserved.
//

#ifndef core_h
#define core_h

#include "array2.h"
#include "vec.h"

static Vec2i cell_offset[] = { Vec2i(-1,0), Vec2i(1,0), Vec2i(0,-1), Vec2i(0,1) }; // cell to cell offset

// BFS markers
enum marked { UNVISITED, VISITED, FINISHED };

// define the type for the marker. 
typedef Array2<marked, Array1<marked>>          Array2m;

#endif /* core_h */

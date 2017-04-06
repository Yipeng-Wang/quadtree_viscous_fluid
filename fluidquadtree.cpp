#include "fluidquadtree.h"

#include <vector>
#include <stack>
#include <set>

FluidQuadTree::FluidQuadTree(float width, Array2f liquid_phi_) : domain_width(width), liquid_phi(liquid_phi_) {
    
    ni = liquid_phi.ni;
    nj = liquid_phi.nj;
    dx = width / (float) ni;
    
    // Test on regular grid.
    //max_depth = 1;
    
    // Test 2 level tree.
    max_depth = 2;
    
    cell_markers.resize(max_depth);
    cellInds.resize(max_depth);
    
    int ni_tmp = ni, nj_tmp = nj;
    for (int depth = max_depth-1; depth >= 0; --depth) {
        cell_markers[depth].resize(ni_tmp, nj_tmp);
        cell_markers[depth].assign(1);
        
        cellInds[depth].resize(ni_tmp, nj_tmp);
        cellInds[depth].assign(-1);
                
        ni_tmp /= 2;
    }
    
    // Generate quad-tree, this only works for 2 level trees.
    // TODO: generalize into multi-level tree, a graded method is also required in this case.
    
    int depth = max_depth - 2;
    
    float threshold = -2.5 * dx;
    
    ni_tmp = get_level_dims(depth);
    for (int i = 0; i < ni_tmp; ++i) for (int j = 0; j < ni_tmp; ++j) {
        
        float phi = 0;
        int i_ind = 2 * i;
        int j_ind = 2 * j;
        for (int offi = 0; offi < 2; ++offi) {
            for (int offj = 0; offj < 2; ++offj) {
                phi += liquid_phi(i_ind + offi, j_ind + offj);
            }
        }
        phi *= 0.25f;
        
        if (phi > threshold) {
            continue;
        }
        
        for (int offi = 0; offi < 2; ++offi) {
            for (int offj = 0; offj < 2; ++offj) {
                cell_markers[depth + 1](i_ind + offi, j_ind + offj) = 0;
            }
        }
    }
    reindex();
}

void FluidQuadTree::reindex() {
    //count leaf nodes, and set up indices
    leaf_cell_count = 0;
    for (int depth = 0; depth < max_depth; ++depth) {
        int ni = get_level_dims(depth);
        for (int i = 0; i < ni; ++i) for (int j = 0; j < ni; ++j) {
            Cell c(depth, i, j);
            
            if (is_leaf_cell(c) && is_cell_in_bounds(c)) { //don't index exterior cells
                leaf_cells.push_back(c);
                cellInds[depth](i, j) = leaf_cell_count;
                ++leaf_cell_count;
            }
        }
    }
    
    active_cell_count = leaf_cell_count;
    // Add extra indices for non-leaf nodes (for alternate axis-aligned solution - *not used for core method*)
    for (int depth = 0; depth < max_depth; ++depth) {
        int ni = get_level_dims(depth);
        for (int i = 0; i < ni; ++i) for (int j = 0; j < ni; ++j) {
            Cell c(depth, i, j);
            if (is_cell_active(c) && !is_leaf_cell(c)) {
                cellInds[depth](i, j) = active_cell_count;
                ++active_cell_count;
            }
        }
    }
    
    // Set up velocity faces in the right places.
    for (int depth = 0; depth < max_depth; ++depth) {
        int ni = get_level_dims(depth);
        for (int i = 0; i < ni; ++i) for (int j = 0; j < ni; ++j) {
            Cell c(depth, i, j);
            if (!is_leaf_cell(c)) {
                continue;
            }
            
            //right
            if (i == ni - 1 || (i < ni - 1 && is_cell_active(right(c)))) {
                Face new_face(depth, i, j, RIGHT);
                qt_faces.push_back(new_face);
            }
            
            //top
            if (j == ni - 1 || (j < ni - 1 && is_cell_active(above(c)))) {
                Face new_face(depth, i, j, TOP);
                qt_faces.push_back(new_face);
            }
            
            //left
            if (i == 0 || (i > 0 && is_cell_active(left(c)) && !is_leaf_cell(left(c)))) {
                //only keep if left neighbour is finer
                Face new_face(depth, i, j, LEFT);
                qt_faces.push_back(new_face);
            }
            
            //bottom
            if (j == 0 || (j > 0 && is_cell_active(below(c)) && !is_leaf_cell(below(c)))) {
                //only keep if bottom neighbour is finer
                Face new_face(depth, i, j, BOTTOM);
                qt_faces.push_back(new_face);
            }
        }
    }
    
    face_count = (int)qt_faces.size();
    
    // Create the cell to face map.
    cell_to_face_map.resize(leaf_cell_count, 4);
    
    for (int f = 0; f < face_count; ++f) {
        Face& face = qt_faces[f];
        
        int forward_i, forward_j;
        int backward_i, backward_j;
        
        FacePosition for_pos;
        FacePosition back_pos;
        
        switch (face.position) {
            case RIGHT:
                forward_i = face.i + 1;
                backward_i = face.i;
                forward_j = face.j;
                backward_j = face.j;
                for_pos = LEFT;
                back_pos = RIGHT;
                break;
            case LEFT:
                forward_i = face.i;
                backward_i = face.i - 1;
                forward_j = face.j;
                backward_j = face.j;
                for_pos = LEFT;
                back_pos = RIGHT;
                break;
            case TOP:
                forward_i = face.i;
                backward_i = face.i;
                forward_j = face.j + 1;
                backward_j = face.j;
                for_pos = BOTTOM;
                back_pos = TOP;
                break;
            case BOTTOM:
                forward_i = face.i;
                backward_i = face.i;
                forward_j = face.j;
                backward_j = face.j - 1;
                for_pos = BOTTOM;
                back_pos = TOP;
                break;
        }
        
        //find all the cells touching this face
        std::stack<std::tuple<Cell, FacePosition>> cells_to_process;
        Cell first_cell(face.depth, forward_i, forward_j);
        Cell second_cell(face.depth, backward_i, backward_j);
        
        if (is_cell_in_bounds(first_cell)) {
            cells_to_process.push(std::make_tuple(first_cell, for_pos));
        }
        
        if (is_cell_in_bounds(second_cell)) {
            cells_to_process.push(std::make_tuple(second_cell, back_pos));
        }
        
        while (cells_to_process.size() > 0) {
            auto& cur_entry = cells_to_process.top();
            cells_to_process.pop();
            Cell& cur_cell = std::get<0>(cur_entry);
            FacePosition& cur_position = std::get<1>(cur_entry);
            
            if (is_leaf_cell(cur_cell)) {
                cell_to_face_map(cellInds[cur_cell.depth](cur_cell.i, cur_cell.j), cur_position) = f;
            } else {
                //push the appropriate 2 children onto the stack
                Cell child0, child1;
                int left_i = 2 * cur_cell.i,
                right_i = 2 * cur_cell.i + 1,
                bottom_j = 2 * cur_cell.j,
                top_j = 2 * cur_cell.j + 1;
                
                switch (cur_position) {
                    case RIGHT:
                        child0 = Cell(cur_cell.depth + 1, right_i, bottom_j);
                        child1 = Cell(cur_cell.depth + 1, right_i, top_j);
                        break;
                    case LEFT:
                        child0 = Cell(cur_cell.depth + 1, left_i, bottom_j);
                        child1 = Cell(cur_cell.depth + 1, left_i, top_j);
                        break;
                    case TOP:
                        child0 = Cell(cur_cell.depth + 1, left_i, top_j);
                        child1 = Cell(cur_cell.depth + 1, right_i, top_j);
                        break;
                    case BOTTOM:
                        child0 = Cell(cur_cell.depth + 1, left_i, bottom_j);
                        child1 = Cell(cur_cell.depth + 1, right_i, bottom_j);
                        break;
                }
                cells_to_process.push(std::make_tuple(child0, cur_position));
                cells_to_process.push(std::make_tuple(child1, cur_position));
            }
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////
    // Set up nodes in the right places
    std::vector<Node> inner_nodes;
    int node_idx = 0;
    node_to_face_map.resize(leaf_cell_count * 4);
    
    cell_to_node_map.resize(leaf_cell_count, 4);
    cell_to_node_map.assign(-1);
    
    for (int depth = 0; depth < max_depth; ++depth) {
        int ni = get_level_dims(depth);
        for (int i = 0; i < ni; ++i) for (int j = 0; j < ni; ++j) {
            Cell c(depth, i, j);
            if (!is_leaf_cell(c)) {
                continue;
            }
            
            // First put in boundary nodes.
            // Nodes on the left boundary.
            if (i == 0 && j != 0) {
                Node left_bound_node(depth, i, j, SW);
                nodes.push_back(left_bound_node);
                int cell_idx = cellInds[left_bound_node.depth](left_bound_node.i, left_bound_node.j);
                node_to_face_map[node_idx].push_back(cell_to_face_map(cell_idx,3));
                // Update cell_to_node_map.
                cell_to_node_map(cell_idx, SW) = node_idx;
                Cell& c = leaf_cells[cell_idx];
                Cell below_c;
                if (is_cell_active(below(c))) {
                    below_c = get_leaf_cell(below(c), NW);
                } else {
                    below_c = get_active_parent(below(c));
                }
                int below_c_idx = cellInds[below_c.depth](below_c.i, below_c.j);
                cell_to_node_map(below_c_idx, NW) = node_idx;
                ++node_idx;
            }
            
            // Nodes on the bottom boundary.
            if (i != 0 && j == 0) {
                Node bottom_bound(depth, i, j, SW);
                nodes.push_back(bottom_bound);
                int cell_idx = cellInds[bottom_bound.depth](bottom_bound.i, bottom_bound.j);
                node_to_face_map[node_idx].push_back(cell_to_face_map(cell_idx, 1));
                // Update cell_to_node_map.
                cell_to_node_map(cell_idx, SW) = node_idx;
                Cell c = leaf_cells[cell_idx];
                Cell left_c;
                if (is_cell_active(left(c))) {
                    left_c = get_leaf_cell(left(c), SE);
                } else {
                    left_c = get_active_parent(left(c));
                }
                int left_c_idx = cellInds[left_c.depth](left_c.i, left_c.j);
                cell_to_node_map(left_c_idx, SE) = node_idx;
                node_idx++;
            }
            
            // Nodes on the right boundary.
            if (i == ni - 1 && j != ni - 1) {
                Node right_bound(depth, i, j, NE);
                nodes.push_back(right_bound);
                int cell_idx = cellInds[right_bound.depth](right_bound.i, right_bound.j);
                node_to_face_map[node_idx].push_back(cell_to_face_map(cell_idx, 2));
                // Update cell_to_node_map.
                cell_to_node_map(cell_idx, NE) = node_idx;
                Cell c = leaf_cells[cell_idx];
                Cell above_c;
                if (is_cell_active(above(c))) {
                    above_c = get_leaf_cell(above(c), SE);
                } else {
                    above_c = get_active_parent(above(c));
                }
                int above_c_idx = cellInds[above_c.depth](above_c.i, above_c.j);
                cell_to_node_map(above_c_idx, SE) = node_idx;
                node_idx++;
            }
            
            // Nodes on the top boundary.
            if (i != ni - 1 && j == ni - 1) {
                Node top_bound(depth, i, j, NE);
                nodes.push_back(top_bound);
                int cell_idx = cellInds[top_bound.depth](top_bound.i, top_bound.j);
                node_to_face_map[node_idx].push_back(cell_to_face_map(cell_idx, 0));
                // Update cell_to_node_map.
                cell_to_node_map(cell_idx, NE) = node_idx;
                Cell c = leaf_cells[cell_idx];
                Cell right_c;
                if (is_cell_active(right(c))) {
                    right_c = get_leaf_cell(right(c), NW);
                } else {
                    right_c = get_active_parent(right(c));
                }
                int right_c_idx = cellInds[right_c.depth](right_c.i, right_c.j);
                cell_to_node_map(right_c_idx, NW) = node_idx;
                ++node_idx;
            }
            
            //collect all inner nodes.
            if (i != ni - 1 && j != ni - 1) {
                // First, add the node at the top right corner of each cell.
                Node top_right_node(depth, i, j, NE);
                inner_nodes.push_back(top_right_node);
                
                // Add nodes at lower right corners when the lower cell is larger.
                if (j != 0 && !is_cell_active(below(c))) {
                    // Suppose the cell c is at the top left corner of the T
                    // junction node.
                    Cell lower_left = get_active_parent(below(c));
                    Cell lower_right = get_active_parent(below(right(c)));
                    if (lower_left == lower_right) {
                        Node t_junction_node(depth, i, j, SE);
                        inner_nodes.push_back(t_junction_node);
                    }
                }
                
                // Add nodes at upper left corners when the left cell is larger.
                if (i != 0 && !is_cell_active(left(c))) {
                    // Suppose the cell c is at the lower right corner of the T
                    // junction node.
                    Cell lower_left = get_active_parent(left(c));
                    Cell upper_left = get_active_parent(left(above(c)));
                    if (lower_left == upper_left) {
                        Node t_junction_node(depth, i, j, NW);
                        inner_nodes.push_back(t_junction_node);
                    }
                }
            }
        }
    }
    
    // Address inner nodes.
    for (Node& node : inner_nodes) {
        // The default position of a node respect to a cell is NE.
        int	LL_i, LL_j, LR_i, LR_j;
        int UL_i, UL_j, UR_i, UR_j;
        
        switch (node.position) {
            case NE:
                LL_i = node.i, LL_j = node.j;
                LR_i = node.i + 1, LR_j = node.j;
                UL_i = node.i, UL_j = node.j + 1;
                UR_i = node.i + 1, UR_j = node.j + 1;
                break;
            case SE:
                LL_i = node.i, LL_j = node.j - 1;
                LR_i = node.i + 1, LR_j = node.j - 1;
                UL_i = node.i, UL_j = node.j;
                UR_i = node.i + 1, UR_j = node.j;
                break;
            case SW:
                LL_i = node.i - 1, LL_j = node.j - 1;
                LR_i = node.i, LR_j = node.j - 1;
                UL_i = node.i - 1, UL_j = node.j;
                UR_i = node.i, UR_j = node.j;
                break;
            case NW:
                LL_i = node.i - 1, LL_j = node.j;
                LR_i = node.i, LR_j = node.j;
                UL_i = node.i - 1, UL_j = node.j + 1;
                UR_i = node.i, UR_j = node.j + 1;
                break;
        }
        
        // Find all the cells touching this node
        Cell LL_cell(node.depth, LL_i, LL_j);
        Cell LR_cell(node.depth, LR_i, LR_j);
        Cell UL_cell(node.depth, UL_i, UL_j);
        Cell UR_cell(node.depth, UR_i, UR_j);
        
        Cell LL_leaf = get_leaf_cell(LL_cell, NE);
        Cell LR_leaf = get_leaf_cell(LR_cell, NW);
        Cell UL_leaf = get_leaf_cell(UL_cell, SE);
        Cell UR_leaf = get_leaf_cell(UR_cell, SW);
        
        nodes.push_back(node);
        
        int LL_leaf_idx = cellInds[LL_leaf.depth](LL_leaf.i, LL_leaf.j);
        int LR_leaf_idx = cellInds[LR_leaf.depth](LR_leaf.i, LR_leaf.j);
        int UL_leaf_idx = cellInds[UL_leaf.depth](UL_leaf.i, UL_leaf.j);
        int UR_leaf_idx = cellInds[UR_leaf.depth](UR_leaf.i, UR_leaf.j);
        // The order the faces is also right, left, top, bottom respect to the node.
        // Update node_to_face_map.
        // Right face.
        if (UR_leaf_idx == LR_leaf_idx) {
            node_to_face_map[node_idx].push_back(-1);
        } else {
            node_to_face_map[node_idx].push_back(cell_to_face_map(UR_leaf_idx, 3));
            if (UL_leaf_idx != UR_leaf_idx) {
                cell_to_node_map(UR_leaf_idx, SW) = node_idx;
            }
            if (LL_leaf_idx != LR_leaf_idx) {
                cell_to_node_map(LR_leaf_idx, NW) = node_idx;
            }
        }
        // Left face.
        if (UL_leaf_idx == LL_leaf_idx) {
            node_to_face_map[node_idx].push_back(-1);
        } else {
            node_to_face_map[node_idx].push_back(cell_to_face_map(LL_leaf_idx, 2));
            if (UL_leaf_idx != UR_leaf_idx) {
                cell_to_node_map(UL_leaf_idx, SE) = node_idx;
            }
            if (LL_leaf_idx != LR_leaf_idx) {
                cell_to_node_map(LL_leaf_idx, NE) = node_idx;
            }
        }
        // Top face.
        if (UL_leaf_idx == UR_leaf_idx) {
            node_to_face_map[node_idx].push_back(-1);
        } else {
            node_to_face_map[node_idx].push_back(cell_to_face_map(UR_leaf_idx, 1));
        }
        // Bottom face.
        if (LL_leaf_idx == LR_leaf_idx) {
            node_to_face_map[node_idx].push_back(-1);
        } else {
            node_to_face_map[node_idx].push_back(cell_to_face_map(LL_leaf_idx, 0));
        }
        ++node_idx;
    }
    node_to_face_map.resize(node_idx);
    node_count = node_idx;
    
    ///////////////////////////////////////////////////////////////////////////
    // Set face_to_node_map and face_to_cell_map.
    for (int i = 0; i < qt_faces.size(); ++i) {
        Face& f = qt_faces[i];
        int n = get_level_dims(f.depth);
        // Skip faces on boundaries.
        if ((f.i == 0 && f.position == LEFT) || (f.i == n - 1 && f.position == RIGHT) || (f.j == 0 && f.position == BOTTOM) || (f.j == n - 1 && f.position == TOP)) {
            continue;
        }
        
        int owner_cell_idx = cellInds[f.depth](f.i, f.j);
        Cell& owner_cell = leaf_cells[owner_cell_idx];
        
        // Set face_to_node_map. Every face has two nodes.
        if (f.position == LEFT || f.position == RIGHT) {
            DiagonalDirection top_position = NW;
            DiagonalDirection bottom_position = SW;
            
            if (f.position == RIGHT) {
                top_position = NE;
                bottom_position = SE;
            }
            int top_node_idx = cell_to_node_map(owner_cell_idx, top_position);
            int bottom_node_idx = cell_to_node_map(owner_cell_idx, bottom_position);
            Vec2i node_inds(top_node_idx, bottom_node_idx);
            face_to_node_map[i] = node_inds;
        }
        
        if (f.position == TOP || f.position == BOTTOM) {
            DiagonalDirection right_position = NE;
            DiagonalDirection left_position = NW;
            
            if (f.position == BOTTOM) {
                right_position = SE;
                left_position = SW;
            }
            int right_node_idx = cell_to_node_map(owner_cell_idx, right_position);
            int left_node_idx = cell_to_node_map(owner_cell_idx, left_position);
            Vec2i node_inds(right_node_idx, left_node_idx);
            face_to_node_map[i] = node_inds;
        }
        
        // Set face_to_cell_map.
        std::stack<Cell> stack;
        std::vector<int> fir_cell_ids;
        std::vector<int> sec_cell_ids;
        if (f.position == LEFT) {
            stack.push(left(owner_cell));
            while (!stack.empty()) {
                Cell c = stack.top();
                stack.pop();
                if (is_leaf_cell(c)) {
                    sec_cell_ids.push_back(cellInds[c.depth](c.i, c.j));
                } else {
                    stack.push(get_child(c, SE));
                    stack.push(get_child(c, NE));
                }
            }
            fir_cell_ids.push_back(owner_cell_idx);
        } else if (f.position == RIGHT) {
            stack.push(right(owner_cell));
            while (!stack.empty()) {
                Cell c = stack.top();
                stack.pop();
                if (is_leaf_cell(c)) {
                    fir_cell_ids.push_back(cellInds[c.depth](c.i, c.j));
                } else {
                    stack.push(get_child(c, SW));
                    stack.push(get_child(c, NW));
                    
                }
            }
            sec_cell_ids.push_back(owner_cell_idx);
        } else if (f.position == TOP) {
            stack.push(above(owner_cell));
            while (!stack.empty()) {
                Cell c = stack.top();
                stack.pop();
                if (is_leaf_cell(c)) {
                    fir_cell_ids.push_back(cellInds[c.depth](c.i, c.j));
                } else {
                    stack.push(get_child(c, SW));
                    stack.push(get_child(c, SE));
                }
            }
            sec_cell_ids.push_back(owner_cell_idx);
        } else {
            stack.push(below(owner_cell));
            while (!stack.empty()) {
                Cell c = stack.top();
                stack.pop();
                if (is_leaf_cell(c)) {
                    sec_cell_ids.push_back(cellInds[c.depth](c.i, c.j));
                } else {
                    stack.push(get_child(c, NW));
                    stack.push(get_child(c, NE));
                }
            }
            fir_cell_ids.push_back(owner_cell_idx);
        }
        face_to_cell_map[i] = std::make_pair(fir_cell_ids, sec_cell_ids);
    }
}

Cell FluidQuadTree::above(const Cell& c) {
    return Cell(c.depth, c.i, c.j + 1);
}

Cell FluidQuadTree::below(const Cell& c) {
    return Cell(c.depth, c.i, c.j - 1);
}

Cell FluidQuadTree::right(const Cell& c) {
    return Cell(c.depth, c.i + 1, c.j);
}

Cell FluidQuadTree::left(const Cell& c) {
    return Cell(c.depth, c.i - 1, c.j);
}

Cell FluidQuadTree::child(const Cell& c) {
    return Cell(c.depth + 1, 2 * c.i, 2 * c.j);
}

Cell FluidQuadTree::get_active_parent(Cell c) {
    assert(is_cell_in_bounds(c));
    //Note: This may return the cell itself, and that is indeed desired behavior in some places.
    while (c.depth > 0 && !is_cell_active(c)) {
        c.depth -= 1;
        c.i /= 2;
        c.j /= 2;
    }
    return c;
}

bool FluidQuadTree::is_leaf_cell(const Cell& c) {
    assert(c.depth >= 0 && c.depth < max_depth);
    
    bool children_dont_exist = c.depth == max_depth-1? true : cell_markers[c.depth + 1](2 * c.i, 2 * c.j) == 0;
    return is_cell_active(c) && children_dont_exist;
}

bool FluidQuadTree::is_cell_in_bounds(const Cell& c) {
    int dim_size = get_level_dims(c.depth);
    return c.i >= 0 && c.j >= 0 && c.i < dim_size && c.j < dim_size;
}

bool FluidQuadTree::is_cell_active(const Cell& c) {
    assert(c.depth >= 0 && c.depth < max_depth);
    return cell_markers[c.depth](c.i, c.j) == 1;
}

Vec2f FluidQuadTree::get_cell_centre(const Cell& c) {
    float h = get_cell_width(c.depth);
    float x = (c.i + 0.5) * h;
    float y = (c.j + 0.5) * h;
    return Vec2f(x, y);
}

Vec2f FluidQuadTree::get_face_centre(const Face& f) {
    
    Cell c(f.depth, f.i, f.j);
    float h = get_cell_width(f.depth);
    Vec2f cell_centre = get_cell_centre(c);
    if (f.position == RIGHT)
        return cell_centre + Vec2f(0.5*h, 0);
    else if (f.position == LEFT)
        return cell_centre + Vec2f(-0.5*h, 0);
    else if (f.position == TOP)
        return cell_centre + Vec2f(0, 0.5*h);
    else if (f.position == BOTTOM)
        return cell_centre + Vec2f(0, -0.5*h);
    else {
        assert(false);
        return Vec2f(0, 0);
    }
    
}

Vec2f FluidQuadTree::get_node_position(const Node& n) {
    Cell c(n.depth, n.i, n.j);
    float h = get_cell_width(c.depth);
    Vec2f cell_centre = get_cell_centre(c);
    if (n.position == SW)
        return cell_centre + Vec2f(-0.5*h, -0.5*h);
    else if (n.position == SE)
        return cell_centre + Vec2f(0.5*h, -0.5*h);
    else if (n.position == NW)
        return cell_centre + Vec2f(-0.5*h, 0.5*h);
    else if (n.position == NE)
        return cell_centre + Vec2f(0.5*h, 0.5*h);
    else {
        assert(false);
        return Vec2f(0, 0);
    }
}

int FluidQuadTree::get_level_dims(int level) {
    return cell_markers[level].ni;
}

float FluidQuadTree::get_cell_width(int level) {
    return domain_width / (float)get_level_dims(level);
}

int FluidQuadTree::get_cell_index(const Cell& c) {
    return cellInds[c.depth](c.i, c.j);
}

Cell FluidQuadTree::get_child(const Cell& c, const DiagonalDirection& dir) {
    Cell result(c.depth + 1, c.i * 2, c.j * 2);
    switch (dir) {
        case SW:
            break;
        case SE:
            result.i++;
            break;
        case NW:
            result.j++;
            break;
        case NE:
            result.i++;
            result.j++;
            break;
    }
    return result;
}

Cell FluidQuadTree::get_leaf_cell(const Cell &c, const DiagonalDirection &dir) {
    Cell result = c;
    if (is_cell_active(result)) {
        while (!is_leaf_cell(result)) {
            result = get_child(result, dir);
        }
    } else {
        result = get_active_parent(result);
    }
    return result;
}

Cell FluidQuadTree::get_parent(const Cell& c) {
    Cell result(c.depth - 1, c.i / 2, c.j / 2);
    return result;
}







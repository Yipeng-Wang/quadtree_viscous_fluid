#ifndef POISSON_DOMAIN_H
#define POISSON_DOMAIN_H

#include "vec.h"
#include <limits>

class Domain {
public:
    
    double width;
    Vec2d origin;
    
    Domain(double w, Vec2d o) : origin(o), width(w) {}
    
    virtual double get_width() { return width; }
    virtual Vec2d get_origin() { return origin; }
    
    virtual bool is_in_bounds(const Vec2d& pos) = 0;
    virtual bool find_exit_point(const Vec2d& start, const Vec2d& end, Vec2d& exit_point) = 0;
    
};

class SquareDomain: public Domain {
    
public:
    
    SquareDomain(double width_, Vec2d origin_) : Domain(width_, origin_) {}
    
    bool is_in_bounds(const Vec2d& point) {
        return point[0] >= origin[0] && point[1] >= origin[1] && point[0] <= origin[0] + width && point[1] <= origin[1] + width;
    }
    
    bool find_exit_point(const Vec2d& start, const Vec2d& end, Vec2d& exit_point) {
        
        //for a square we can just check if the endpoint is outside
        //and do the usual test to find the crossing point
        
        bool exits = !is_in_bounds(end);
        if (exits) {
            
            //compute the exit point
            if (end[0] < origin[0]) {
                double t = (origin[0] - end[0]) / (start[0] - end[0]);
                assert(0 < t && t < 1);
                exit_point = lerp(end, start, t);
            }
            else if (end[1] < origin[1]) {
                double t = (origin[1] - end[1]) / (start[1] - end[1]);
                assert(0 < t && t < 1);
                exit_point = lerp(end, start, t);
            }
            else if (end[0] > origin[0] + width) {
                double t = (origin[0] + width - start[0]) / (end[0] - start[0]);
                assert(0 < t && t < 1);
                exit_point = lerp(start, end, t);
            }
            else if (end[1] > origin[1] + width) {
                double t = (origin[1] + width - start[1]) / (end[1] - start[1]);
                assert(0 < t && t < 1);
                exit_point = lerp(start, end, t);
            }
            
        }
        
        return exits;
    }
};

#endif

#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cfloat>

#include "gluvi.h"
#include "fluidsim.h"
#include "openglutils.h"
#include "array2_utils.h"

using namespace std;

//Try changing the grid resolution
int grid_resolution = 100;
float timestep = 0.002;

//Display properties
bool draw_grid = false;
bool draw_particles = true;
bool draw_velocities = false;
bool draw_boundaries = true;

float grid_width = 1;

FluidSim sim;

//Gluvi stuff
//-------------
Gluvi::PanZoom2D cam(-0.1, -0.35, 1.2);
double oldmousetime;
Vec2f oldmouse;
void display();
void mouse(int button, int state, int x, int y);
void drag(int x, int y);
void timer(int junk);

//Boundary definition - several circles in a circular domain.

Vec2f c0(0.5,0.5), c1(0.7,0.5), c2(0.3,0.35), c3(0.5,0.7);
float rad0 = 0.4,  rad1 = 0.1,  rad2 = 0.1,   rad3 = 0.1;
float side = 0.9;

float circle_phi(const Vec2f& position, const Vec2f& centre, float radius) {
    return radius - dist(position, centre);
}

float square_phi(const Vec2f& position, const Vec2f& centre, float side) {
    float x_dis = fabs(position[0]-centre[0]);
    float y_dis = fabs(position[1]-centre[1]);
    float half_side = side * 0.5;
    if (x_dis <= half_side && y_dis <= half_side) {
        return x_dis > y_dis? half_side - x_dis : half_side - y_dis;
    } else if (x_dis >= half_side && y_dis < half_side) {
        return half_side - x_dis;
    } else if (y_dis >= half_side && x_dis < half_side) {
        return half_side - y_dis;
    } else {
        return -sqrt(sqr(x_dis-half_side) + sqr(y_dis-half_side));
    }
}

float boundary_phi(const Vec2f& position) {
    float phi0 = square_phi(position, c0, side);
    
    //   float phi0 = circle_phi(position, c0, rad0);
    //   float phi1 = circle_phi(position, c1, rad1);
    //   float phi2 = circle_phi(position, c2, rad2);
    //   float phi3 = circle_phi(position, c3, rad3);
    
    return phi0;//min(min(phi0,phi1),min(phi2,phi3));
}

//Main testing code
//-------------
int main(int argc, char **argv)
{
    
    //Setup viewer stuff
    Gluvi::init("GFM Free Surface Liquid Solver with Static Variational Boundaries", &argc, argv);
    Gluvi::camera=&cam;
    Gluvi::userDisplayFunc=display;
    Gluvi::userMouseFunc=mouse;
    Gluvi::userDragFunc=drag;
    glClearColor(1,1,1,1);
    
    glutTimerFunc(1000, timer, 0);
    
    //Set up the simulation
    sim.initialize(grid_width, grid_resolution, grid_resolution);
    
    //set up a circle boundary
    sim.set_boundary(boundary_phi);
    
    //Stick some liquid particles in the domain
    int offset = 0;
    for(int i = 0; i < sqr(grid_resolution); ++i) {
        for(int parts = 0; parts < 3; ++parts) {
            float x = randhashf(++offset, 0,1);
            float y = randhashf(++offset, 0,1);
            Vec2f pt(x,y);
            
            //add a column (for buckling) and a beam (for bending) and a disk (for rolling and flowing)
            if (boundary_phi(pt) > 0 && ((pt[0] > 0.42 && pt[0] < 0.46)
                || (pt[0] < 0.36 && pt[1] > 0.45 && pt[1] < 0.5)
                || circle_phi(pt, Vec2f(0.7, 0.65), 0.15) > 0)) {
                sim.add_particle(pt);
            }
            
            
            // filled half of the tank with fluid. Check if the signed distance is correct.
//            if (boundary_phi(pt) > 0 && pt[1] <= 0.5) {
//                sim.add_particle(pt);
//            }
        }
    }
    
    //sim.advance(timestep);
    
    Gluvi::run();
    return 0;
}


void display(void)
{
    
    if(draw_grid) {
        glColor3f(0,0,0);
        glLineWidth(0.5);
        draw_grid2d(Vec2f(0,0), sim.dx, sim.ni, sim.nj);
    }
    
    if(draw_boundaries) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        // draw a circular boundary.
        //draw_circle2d(c0, rad0, 50);
        
        // draw square boundary.
        draw_box2d(c0 - Vec2f(0.5*side, 0.5*side), side, side);
    }
    
    if(draw_particles) {
        glColor3f(0,0,1);
        glPointSize(3);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        draw_points2d(sim.particles);
        for(unsigned int p = 0; p < sim.particles.size(); ++p) {
            draw_circle2d(sim.particles[p], sim.particle_radius, 20);
        }
    }
    
    if(draw_velocities) {
        glColor3f(1,0,0);
        for(int j = 0;j < sim.nj; ++j) for(int i = 0; i < sim.ni; ++i) {
            Vec2f pos((i+0.5)*sim.dx,(j+0.5)*sim.dx);
            draw_arrow2d(pos, pos + 0.01f*sim.get_velocity(pos), 0.01*sim.dx);
        }
    }
}

void mouse(int button, int state, int x, int y)
{
    Vec2f newmouse;
    cam.transform_mouse(x, y, newmouse.v);
    //double newmousetime=get_time_in_seconds();
    
    oldmouse=newmouse;
    //oldmousetime=newmousetime;
}

void drag(int x, int y)
{
    Vec2f newmouse;
    cam.transform_mouse(x, y, newmouse.v);
    //double newmousetime=get_time_in_seconds();
    
    oldmouse=newmouse;
    //oldmousetime=newmousetime;
}


void timer(int junk)
{
    
    sim.advance(timestep);
    
    glutPostRedisplay();
    glutTimerFunc(1, timer, 0);
    
}

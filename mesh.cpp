#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "mesh.h"

int count_points1;
int count_internal_points1;
int count_internal_quadrangles1;
int count_internal_surfaces1;
int count_boundary_quadrangles1;
int count_boundary_surfaces1;
int count_elements1;

mesh::point1* points1;
mesh::internal_surface1* isurfaces1;
mesh::boundary_surface1* bsurfaces1;
mesh::hexahedron1* hexahedrons1;

int count_points2;
int count_internal_points2;
int count_internal_quadrangles2;
int count_internal_surfaces2;
int count_boundary_quadrangles2;
int count_boundary_surfaces2;
int count_elements2;

mesh::point2* points2;
mesh::internal_surface2* isurfaces2;
mesh::boundary_surface2* bsurfaces2;
mesh::hexahedron2* hexahedrons2;

bool mesh::LoadMesh(std::string p_mesh_file)
{
    std::ifstream ifile(p_mesh_file.c_str());
    if (!ifile)
        return false;
    std::stringstream ar;
    std::string str;

    std::getline(ifile, str);
    ar.str(str);
    ar >> count_points1 >> count_internal_points1;
    ar.clear();
    points1 = new point1[count_points1];
    for (int i = 0, n = 0; i < count_points1; ++i)
    {
        std::getline(ifile, str);
        ar.str(str);
        ar >> points1[i].front;
        ar >> points1[i].right;
        ar >> points1[i].up;
        ar >> points1[i].back;
        ar >> points1[i].left;
        ar >> points1[i].down;
        ar >> points1[i].x;
        ar >> points1[i].y;
        ar >> points1[i].z;
        ar >> points1[i].distance;
        ar >> points1[i].entry_or_boundary_condition;
        ar >> points1[i].surface;
        ar >> points1[i].is_node_on_real_border;
        ar >> points1[i].is_node_on_hard_border;
        if (points1[i].is_node_on_real_border)
            points1[i].output_number = ++n;
        switch (points1[i].entry_or_boundary_condition)
        {
        case 10:
            points1[i].type = FreedomExit;
            break;
        case 20:
            points1[i].type = HardBorder;
            break;
        case 30:
            points1[i].type = Symmetry;
            break;
        case 40:
            points1[i].type = Custom;
            break;
        default:
            break;
        }
        ar.clear();
    }

    std::getline(ifile, str);
    ar.str(str);
    ar >> count_internal_quadrangles1 >> count_internal_surfaces1;
    ar.clear();
    isurfaces1 = new internal_surface1[count_internal_surfaces1];
    int temp = 0;
    for (int i = 0; i < count_internal_surfaces1; ++i)
    {
        std::getline(ifile, str);
        ar.str(str);
        ar >> isurfaces1[i].A;
        ar >> isurfaces1[i].B;
        ar >> isurfaces1[i].C;
        ar >> isurfaces1[i].D;
        ar >> isurfaces1[i].area;
        ar >> isurfaces1[i].x;
        ar >> isurfaces1[i].y;
        ar >> isurfaces1[i].z;
        ar >> isurfaces1[i].distance;
        ar >> isurfaces1[i].nx;
        ar >> isurfaces1[i].ny;
        ar >> isurfaces1[i].nz;
        ar >> isurfaces1[i].element1;
        ar >> temp;
        ar >> isurfaces1[i].element2;
        ar >> temp;
        ar >> isurfaces1[i].is_two_elements;
        ar.clear();
    }

    std::getline(ifile, str);
    ar.str(str);
    ar >> count_boundary_quadrangles1 >> count_boundary_surfaces1;
    ar.clear();
    bsurfaces1 = new boundary_surface1[count_boundary_surfaces1];
    for (int i = 0; i < count_boundary_surfaces1; ++i)
    {
        std::getline(ifile, str);
        ar.str(str);
        ar >> bsurfaces1[i].A;
        ar >> bsurfaces1[i].B;
        ar >> bsurfaces1[i].C;
        ar >> bsurfaces1[i].D;
        ar >> bsurfaces1[i].area;
        ar >> bsurfaces1[i].x;
        ar >> bsurfaces1[i].y;
        ar >> bsurfaces1[i].z;
        ar >> bsurfaces1[i].distance;
        ar >> bsurfaces1[i].nx;
        ar >> bsurfaces1[i].ny;
        ar >> bsurfaces1[i].nz;
        ar >> bsurfaces1[i].element1;
        ar >> temp;
        ar >> temp;
        ar >> bsurfaces1[i].boundary_condition;
        switch (bsurfaces1[i].boundary_condition)
        {
        case 10:
            bsurfaces1[i].type = FreedomExit;
            break;
        case 20:
            bsurfaces1[i].type = HardBorder;
            break;
        case 30:
            bsurfaces1[i].type = Symmetry;
            break;
        case 40:
            bsurfaces1[i].type = Custom;
            break;
        default:
            break;
        }
        ar.clear();
    }

    std::getline(ifile, str);
    ar.str(str);
    ar >> count_elements1;
    ar.clear();
    hexahedrons1 = new hexahedron1[count_elements1];
    for (int i = 0; i < count_elements1; ++i)
    {
        std::getline(ifile, str);
        ar.str(str);
        ar >> hexahedrons1[i].A;
        ar >> hexahedrons1[i].B;
        ar >> hexahedrons1[i].C;
        ar >> hexahedrons1[i].D;
        ar >> hexahedrons1[i].E;
        ar >> hexahedrons1[i].F;
        ar >> hexahedrons1[i].G;
        ar >> hexahedrons1[i].H;
        ar >> hexahedrons1[i].volume;
        ar >> temp;
        ar.clear();
    }

    std::getline(ifile, str);
    ar.str(str);
    ar >> count_points2 >> count_internal_points2;
    ar.clear();
    points2 = new point2[count_points2];
    for (int i = 0; i < count_points2; ++i)
    {
        std::getline(ifile, str);
        ar.str(str);
        ar >> points2[i].neighbour;
        ar >> points2[i].x;
        ar >> points2[i].y;
        ar >> points2[i].z;
        ar.clear();
    }

    std::getline(ifile, str);
    ar.str(str);
    ar >> count_internal_quadrangles2 >> count_internal_surfaces2;
    ar.clear();
    isurfaces2 = new internal_surface2[count_internal_surfaces2];
    for (int i = 0; i < count_internal_surfaces2; ++i)
    {
        std::getline(ifile, str);
        ar.str(str);
        ar >> isurfaces2[i].A;
        ar >> isurfaces2[i].B;
        ar >> isurfaces2[i].C;
        ar >> isurfaces2[i].D;
        ar >> isurfaces2[i].area;
        ar >> isurfaces2[i].nx;
        ar >> isurfaces2[i].ny;
        ar >> isurfaces2[i].nz;
        ar >> isurfaces2[i].element1;
        ar >> isurfaces2[i].element2;
        ar.clear();
    }

    std::getline(ifile, str);
    ar.str(str);
    ar >> count_boundary_quadrangles2 >> count_boundary_surfaces2;
    ar.clear();
    bsurfaces2 = new boundary_surface2[count_boundary_surfaces2];
    for (int i = 0; i < count_boundary_surfaces2; ++i)
    {
        std::getline(ifile, str);
        ar.str(str);
        ar >> bsurfaces2[i].A;
        ar >> bsurfaces2[i].B;
        ar >> bsurfaces2[i].C;
        ar >> bsurfaces2[i].D;
        ar >> bsurfaces2[i].area;
        ar >> bsurfaces2[i].nx;
        ar >> bsurfaces2[i].ny;
        ar >> bsurfaces2[i].nz;
        ar >> bsurfaces2[i].element1;
        ar >> bsurfaces2[i].boundary_condition;
        switch (bsurfaces2[i].boundary_condition)
        {
        case 10:
            bsurfaces2[i].type = FreedomExit;
            break;
        case 20:
            bsurfaces2[i].type = HardBorder;
            break;
        case 30:
            bsurfaces2[i].type = Symmetry;
            break;
        case 40:
            bsurfaces2[i].type = Custom;
            break;
        default:
            break;
        }
        ar.clear();
    }

    std::getline(ifile, str);
    ar.str(str);
    ar >> count_elements2;
    ar.clear();
    hexahedrons2 = new hexahedron2[count_elements2];
    for (int i = 0; i < count_elements2; ++i)
    {
        std::getline(ifile, str);
        ar.str(str);
        ar >> hexahedrons2[i].volume;
        ar.clear();
    }
    ifile.close();
    return true;
}

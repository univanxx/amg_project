#ifndef EULER_H
#define EULER_H

#include <map>
#include <string>
#include "main.h"
#include "mesh.h"
#include "vector"

extern int count_points1;
extern int count_internal_points1;
extern int count_internal_quadrangles1;
extern int count_internal_surfaces1;
extern int count_boundary_quadrangles1;
extern int count_boundary_surfaces1;
extern int count_elements1;

extern mesh::point1* points1;
extern mesh::internal_surface1* isurfaces1;
extern mesh::boundary_surface1* bsurfaces1;
extern mesh::hexahedron1* hexahedrons1;

extern int count_points2;
extern int count_internal_points2;
extern int count_internal_quadrangles2;
extern int count_internal_surfaces2;
extern int count_boundary_quadrangles2;
extern int count_boundary_surfaces2;
extern int count_elements2;

extern mesh::point2* points2;
extern mesh::internal_surface2* isurfaces2;
extern mesh::boundary_surface2* bsurfaces2;
extern mesh::hexahedron2* hexahedrons2;


class Euler
{
public:

    void SaveMeshInGMSHFile();
    void SaveSolutionInGMSHFile();

    // консервативные переменные
    double* m_u0, * m_u1, * m_u2, * m_u3, * m_u4;
    // примитивные переменные
    double* m_rho, * m_w1, * m_w2, * m_w3, * pressure;
    // универсальная газовая постоянная [Дж / (кг * К)]
    double m_R = 286.7;

protected:

    // пересчитывает консервативные переменные в примитивные
    void Decoding();

    // gmsh-файл в который сохраняется решение
    std::string m_gmsh_file = "C:/Users/wchhi/source/repos/EulerProj/EulerProj/results/res.msh";
    // ДЛЯ ЗАДАЧИ: время = 1, шаг = 0.1
    double m_time_moment = 1;
    uint32_t m_step = 0.1;
};

#endif // EULER_H



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
    // консервативные переменные
    double* m_u0, * m_u1, * m_u2, * m_u3, * m_u4;
    // примитивные переменные
    double* m_rho, * m_w1, * m_w2, * m_w3, * pressure;
    // универсальная газовая постоянная [Дж / (кг * К)]
    double m_R = 286.7;
    // ГУ
    double* p_lim, * ro_lim, * w1_lim, * w2_lim, * w3_lim;
	void SaveMeshInGMSHFile();
	void SaveSolutionInGMSHFile();
	// Конструктор
	Euler()
	{
		m_u0 = new double[count_elements1];
		m_u1 = new double[count_elements1];
		m_u2 = new double[count_elements1];
		m_u3 = new double[count_elements1];
		m_u4 = new double[count_elements1];

        p_lim = new double[count_elements1];
        ro_lim = new double[count_elements1];
        w1_lim = new double[count_elements1];
        w2_lim = new double[count_elements1];
        w3_lim = new double[count_elements1];

        int count_real_border_points = 0; // количество точек на реальной поверхности
        for (int i = 0; i < count_elements1; ++i)
            if (points1[i].is_node_on_real_border)
                ++count_real_border_points;

        m_rho = new double[count_real_border_points];
        m_w1 = new double[count_real_border_points];
        m_w2 = new double[count_real_border_points];
        m_w3 = new double[count_real_border_points];
        pressure = new double[count_real_border_points];

		for (int j = 0; j < count_elements1; ++j)
		{
			m_u0[j] = 0.0;
			m_u1[j] = 0.0;
			m_u2[j] = 0.0;
			m_u3[j] = 0.0;
			m_u4[j] = 0.0;

            ro_lim[j] = 0.0;
            p_lim[j] = 0.0;
            w1_lim[j] = 0.0;
            w2_lim[j] = 0.0;
            w3_lim[j] = 0.0;
		}

        for (int j = 0; j < count_real_border_points; ++j)
        {
            m_rho[j] = 0.0;
            m_w1[j] = 0.0;
            m_w2[j] = 0.0;
            m_w3[j] = 0.0;
            pressure[j] = 0.0;
        }
	}
    ~Euler() 
    {
        delete[] m_rho;
        delete[] m_w1;
        delete[] m_w2;
        delete[] m_w3;
        delete[] pressure;

        delete[] m_u0;
        delete[] m_u1;
        delete[] m_u2;
        delete[] m_u3;
        delete[] m_u4;

        delete[] ro_lim;
        delete[] w1_lim;
        delete[] w2_lim;
        delete[] w3_lim;
        delete[] p_lim;
    }
protected:

    // пересчитывает консервативные переменные в примитивные
    void Decoding();

    // gmsh-файл в который сохраняется решение
    std::string m_gmsh_file = "C:/Users/wchhi/source/repos/EulerProj/EulerProj/results/res.msh";
	//std::string m_gmsh_file = "C:/Users/Asus/Documents/Visual Studio 2013/Projects/EulerProject/results/res.msh";
    // ДЛЯ ЗАДАЧИ: время = 1, шаг = ?
    double m_time_moment = 1;
    int m_step = -1;
};

#endif // EULER_H



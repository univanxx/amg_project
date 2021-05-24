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
    // êîíñåðâàòèâíûå ïåðåìåííûå
    double* m_u0, * m_u1, * m_u2, * m_u3, * m_u4;
    // ïðèìèòèâíûå ïåðåìåííûå
    double* m_rho, * m_w1, * m_w2, * m_w3, * pressure;
    // óíèâåðñàëüíàÿ ãàçîâàÿ ïîñòîÿííàÿ [Äæ / (êã * Ê)]
    double m_R = 286.7;

	void SaveMeshInGMSHFile();
	void SaveSolutionInGMSHFile();
	// Êîíñòðóêòîð
	Euler()
	{
		m_u0 = new double[count_elements1];
		m_u1 = new double[count_elements1];
		m_u2 = new double[count_elements1];
		m_u3 = new double[count_elements1];
		m_u4 = new double[count_elements1];

		m_rho = new double[count_elements1];
		m_w1 = new double[count_elements1];
		m_w2 = new double[count_elements1];
		m_w3 = new double[count_elements1];
		pressure = new double[count_elements1];

		for (int j = 0; j < count_elements1; ++j)
		{
			m_u0[j] = 0.0;
			m_u1[j] = 0.0;
			m_u2[j] = 0.0;
			m_u3[j] = 0.0;
			m_u4[j] = 0.0;

			m_rho[j] = 0.0;
			m_w1[j] = 0.0;
			m_w2[j] = 0.0;
			m_w3[j] = 0.0;
			pressure[j] = 0.0;
		}
	}
  
    // êîíñåðâàòèâíûå ïåðåìåííûå
    double* m_u0, * m_u1, * m_u2, * m_u3, * m_u4;
    // ïðèìèòèâíûå ïåðåìåííûå
    double* m_rho, * m_w1, * m_w2, * m_w3, * pressure;
    // óíèâåðñàëüíàÿ ãàçîâàÿ ïîñòîÿííàÿ [Äæ / (êã * Ê)]
    double m_R = 286.7;

protected:

    // ïåðåñ÷èòûâàåò êîíñåðâàòèâíûå ïåðåìåííûå â ïðèìèòèâíûå
    void Decoding();

    // gmsh-ôàéë â êîòîðûé ñîõðàíÿåòñÿ ðåøåíèå
    std::string m_gmsh_file = "C:/Users/wchhi/source/repos/EulerProj/EulerProj/results/res.msh";
	  //std::string m_gmsh_file = "C:/Users/Asus/Documents/Visual Studio 2013/Projects/EulerProject/results/res.msh";
    // ÄËß ÇÀÄÀ×È: âðåìÿ = 1, øàã = 0.05
    double m_time_moment = 1;
    double m_step = 0.05;
    // ÄËß ÇÀÄÀ×È: âðåìÿ = 1, øàã = 0.1
    double m_time_moment = 1;
    uint32_t m_step = 0.1;

};

#endif // EULER_H



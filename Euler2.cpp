//#include "stdafx.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include "Euler2.h"
// #include "vector"

void Euler::Decoding()
{
    for (int j = 0; j < count_elements1; ++j)
    {
        m_rho[j] = m_u0[j];
        m_w1[j] = m_u1[j] / m_rho[j];
        m_w2[j] = m_u2[j] / m_rho[j];
        m_w3[j] = m_u3[j] / m_rho[j];
        pressure[j] = 0.4 * (m_u4[j] - m_rho[j] * 0.5 * (m_w1[j] * m_w1[j] + m_w2[j] * m_w2[j] + m_w3[j] * m_w3[j]));
    }
}

void Euler::SaveMeshInGMSHFile()
{
    std::ostringstream ar;
    int N = 0; // количество точек на реальной поверхности
    for (int i = 0; i < count_points1; ++i)
        if (points1[i].is_node_on_real_border)
            ++N;
    ar << "$MeshFormat" << std::endl;
    ar << "2.0" << " " << 0 << " " << sizeof(double) << std::endl;
    ar << "$EndMeshFormat" << std::endl;
    ar << "$Nodes" << std::endl;
    ar << N << std::endl;
    for (int j = 0, n = 0; j < count_points1; ++j)
    {
        if (points1[j].is_node_on_real_border)
        {
            ar << ++n << " ";
            ar << points1[j].x << " " << points1[j].y << " " << points1[j].z << '\n';
        }
    }
    ar << "$EndNodes" << std::endl;
    ar << "$Elements" << std::endl;
    ar << count_boundary_surfaces2 << std::endl;
    double xA, yA, zA, xB, yB, zB, xC, yC, zC, x1, y1, z1, x2, y2, z2;
    double nx, ny, nz;
    for (int j = 0; j < count_boundary_surfaces2; ++j)
    {
        xA = points1[bsurfaces2[j].A].x;
        yA = points1[bsurfaces2[j].A].y;
        zA = points1[bsurfaces2[j].A].z;

        xB = points1[bsurfaces2[j].B].x;
        yB = points1[bsurfaces2[j].B].y;
        zB = points1[bsurfaces2[j].B].z;

        xC = points1[bsurfaces2[j].C].x;
        yC = points1[bsurfaces2[j].C].y;
        zC = points1[bsurfaces2[j].C].z;

        x1 = xB - xA;
        y1 = yB - yA;
        z1 = zB - zA;

        x2 = xC - xB;
        y2 = yC - yB;
        z2 = zC - zB;

        nx = bsurfaces2[j].nx;
        ny = bsurfaces2[j].ny;
        nz = bsurfaces2[j].nz;

        // запись четырёхугольников
        if (j < count_boundary_quadrangles2)
        {
            if (nx * (y1 * z2 - y2 * z1) - ny * (x1 * z2 - x2 * z1) + nz * (x1 * y2 - x2 * y1) > 0.0)
                ar << j + 1 << " " << 3 << " " << 0 << " " << points1[bsurfaces2[j].A].output_number << " " << points1[bsurfaces2[j].B].output_number << \
                " " << points1[bsurfaces2[j].C].output_number << " " << points1[bsurfaces2[j].D].output_number << std::endl;
            else
                ar << j + 1 << " " << 3 << " " << 0 << " " << points1[bsurfaces2[j].A].output_number << " " << points1[bsurfaces2[j].D].output_number << \
                " " << points1[bsurfaces2[j].C].output_number << " " << points1[bsurfaces2[j].B].output_number << std::endl;
        }
        // запись треугольников
        else
        {
            if (nx * (y1 * z2 - y2 * z1) - ny * (x1 * z2 - x2 * z1) + nz * (x1 * y2 - x2 * y1) > 0.0)
                ar << j + 1 << " " << 2 << " " << 0 << " " << points1[bsurfaces2[j].A].output_number << " " << points1[bsurfaces2[j].B].output_number << \
                " " << points1[bsurfaces2[j].C].output_number << std::endl;
            else
                ar << j + 1 << " " << 2 << " " << 0 << " " << points1[bsurfaces2[j].A].output_number << " " << points1[bsurfaces2[j].C].output_number << \
                " " << points1[bsurfaces2[j].B].output_number << std::endl;
        }
    }
    ar << "$EndElements" << std::endl;
    std::ofstream ofst(m_gmsh_file.c_str());
    ofst << ar.str();
    ofst.close();
}

void Euler::SaveSolutionInGMSHFile()
{
    Decoding();
    static int N = 0; // количество точек на реальной поверхности
    if (!N)
        for (int i = 0; i < count_points1; ++i)
            if (points1[i].is_node_on_real_border)
                ++N;
    std::ostringstream ar;
    ar << "$NodeData" << std::endl;
    ar << 1 << std::endl;
    ar << "\"" << "density" << "\"" << std::endl;
    ar << 1 << std::endl;
    ar << m_time_moment << std::endl;
    ar << 3 << std::endl;
    ar << m_step << std::endl;
    ar << 1 << std::endl;
    ar << N << std::endl;
    for (int j = 0; j < N; j++)
    {
        ar << j + 1 << " ";
        ar << m_rho[j] << '\n';
    }
    ar << "$EndNodeData" << std::endl;

    ar << "$NodeData" << std::endl;
    ar << 1 << std::endl;
    ar << "\"" << "Vx" << "\"" << std::endl;
    ar << 1 << std::endl;
    ar << m_time_moment << std::endl;
    ar << 3 << std::endl;
    ar << m_step << std::endl;
    ar << 1 << std::endl;
    ar << N << std::endl;
    for (int j = 0; j < N; j++)
    {
        ar << j + 1 << " ";
        ar << m_w1[j] << '\n';
    }
    ar << "$EndNodeData" << std::endl;

    ar << "$NodeData" << std::endl;
    ar << 1 << std::endl;
    ar << "\"" << "Vy" << "\"" << std::endl;
    ar << 1 << std::endl;
    ar << m_time_moment << std::endl;
    ar << 3 << std::endl;
    ar << m_step << std::endl;
    ar << 1 << std::endl;
    ar << N << std::endl;
    for (int j = 0; j < N; j++)
    {
        ar << j + 1 << " ";
        ar << m_w2[j] << '\n';
    }
    ar << "$EndNodeData" << std::endl;

    ar << "$NodeData" << std::endl;
    ar << 1 << std::endl;
    ar << "\"" << "Vz" << "\"" << std::endl;
    ar << 1 << std::endl;
    ar << m_time_moment << std::endl;
    ar << 3 << std::endl;
    ar << m_step << std::endl;
    ar << 1 << std::endl;
    ar << N << std::endl;
    for (int j = 0; j < N; j++)
    {
        ar << j + 1 << " ";
        ar << m_w3[j] << '\n';
    }
    ar << "$EndNodeData" << std::endl;

    ar << "$NodeData" << std::endl;
    ar << 1 << std::endl;
    ar << "\"" << "pressure" << "\"" << std::endl;
    ar << 1 << std::endl;
    ar << m_time_moment << std::endl;
    ar << 3 << std::endl;
    ar << m_step << std::endl;
    ar << 1 << std::endl;
    ar << N << std::endl;
    for (int j = 0; j < N; j++)
    {
        ar << j + 1 << " ";
        ar << pressure[j] << '\n';
    }
    ar << "$EndNodeData" << std::endl;

    ar << "$NodeData" << std::endl;
    ar << 1 << std::endl;
    ar << "\"" << "temperature" << "\"" << std::endl;
    ar << 1 << std::endl;
    ar << m_time_moment << std::endl;
    ar << 3 << std::endl;
    ar << m_step << std::endl;
    ar << 1 << std::endl;
    ar << N << std::endl;
    for (int j = 0; j < N; j++)
    {
        ar << j + 1 << " ";
        ar << pressure[j] / (m_rho[j] * m_R) << '\n';
    }
    ar << "$EndNodeData" << std::endl;

    std::ofstream ofst(m_gmsh_file.c_str(), std::ios_base::app);
    ofst << ar.str();
    ofst.close();

    delete[] m_rho;
    delete[] m_w1;
    delete[] m_w2;
    delete[] m_w3;
    delete[] pressure;
}

#include <iostream>
using namespace std;
#include "main.h"
#include "mesh.h"
#include "Euler2.h"
#include "vector"

extern mesh::point1* points1;
extern mesh::internal_surface1* isurfaces1;
extern mesh::boundary_surface1* bsurfaces1;
extern mesh::hexahedron1* hexahedrons1;

extern mesh::point2* points2;
extern mesh::internal_surface2* isurfaces2;
extern mesh::boundary_surface2* bsurfaces2;
extern mesh::hexahedron2* hexahedrons2;

int main()
{
    setlocale(LC_ALL, "Russian");
    mesh::LoadMesh("C:/Users/wchhi/Downloads/mymesh.txt");
    cout << "Информация о внутренних рёбрах:" << endl;
    for (int i = 0; i < count_internal_surfaces1; ++i)
    {
        cout << isurfaces1[i].A << " ";
        cout << isurfaces1[i].B << " ";
        cout << isurfaces1[i].C << " ";
        cout << isurfaces1[i].D << " ";
        cout << "Площадь соотв. элемента: " << isurfaces1[i].area << " ";
        cout << "Координаты: " << endl;
        cout << isurfaces1[i].x << " ";
        cout << isurfaces1[i].y << " ";
        cout << isurfaces1[i].z << " ";
        //cout << isurfaces1[i].distance << " ";
        cout << isurfaces1[i].nx << " ";
        cout << isurfaces1[i].ny << " ";
        cout << isurfaces1[i].nz << " ";
        cout << "Элемент, которому принадлежит эта граница:" << isurfaces1[i].element1 << endl;
        cout << "Элемент, которому принадлежит эта граница (соседний):" << isurfaces1[i].element2 << endl;
        // Задача распада-разрыва - 1 и 0 отвечают за области
        cout << "Номер ГУ или НУ элемента, к которому принадлежит эта грань:" << points1[isurfaces1[i].element1].entry_or_boundary_condition << endl;
    }
    cout << "Информация о внешних рёбрах:" << endl;
    for (int i = 0; i < 10; ++i)
    {
        cout << bsurfaces1[i].A << " ";
        cout << bsurfaces1[i].B << " ";
        cout << bsurfaces1[i].C << " ";
        cout << bsurfaces1[i].D << " ";
        cout << "Площадь соотв. элемента: " << bsurfaces1[i].area << " ";
        cout << "Координаты: " << endl;
        cout << bsurfaces1[i].x << " ";
        cout << bsurfaces1[i].y << " ";
        cout << bsurfaces1[i].z << " ";
        //cout << bsurfaces1[i].distance << " ";
        cout << bsurfaces1[i].nx << " ";
        cout << bsurfaces1[i].ny << " ";
        cout << bsurfaces1[i].nz << " ";
        cout << "Элемент, которому принадлежит эта граница:" << bsurfaces1[i].element1 << endl;
        cout << "Номер ГУ или НУ элемента, к которому принадлежит эта грань:" << points1[bsurfaces1[i].element1].entry_or_boundary_condition << endl;
        cout << bsurfaces1[i].boundary_condition << " ";
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
        cout << bsurfaces1[i].type << endl;
    }


    cout << "Начальные условия:" << endl;
    vector<vector<double> > primitive(count_points1, vector<double>(5));
    //vector<vector<vector<double>> > zero_cond(count_points1, vector<vector<double>>(3, vector<double>(3)));
    for (int i = 0; i < count_points1; ++i)
    {
        if (points1[i].entry_or_boundary_condition == 0)
        {
            primitive[i][0] = 1;  // ro
            primitive[i][1] = 2;  // u
            primitive[i][2] = 2;  // v 
            primitive[i][3] = 2;  // w
            primitive[i][4] = 0.5;  // p
        }
        else if (points1[i].entry_or_boundary_condition == 1)
        {
            primitive[i][0] = 0;
            primitive[i][1] = 1;
            primitive[i][2] = 1;
            primitive[i][3] = 1;
            primitive[i][4] = 2;
        }
    }

    // Решение
    vector<vector<vector<double>> > U(10, vector<vector<double>>(count_points1, vector<double>(5)));
    for (int i = 0; i < 10; ++i)
    {
        for (int j = 0; j < count_points1; ++j)
        {
            U[i][j][0] = primitive[j][0];
            U[i][j][1] = primitive[j][1] / U[i][j][0];
            U[i][j][2] = primitive[j][2] / U[i][j][0];
            U[i][j][3] = primitive[j][3] / U[i][j][0];
            U[i][j][4] = 4*(primitive[j][4] - 0.5*(primitive[j][1]*primitive[j][1] + primitive[j][2]*primitive[j][2] + primitive[j][3]*primitive[j][3]))/10;
        }
    }

    Euler efvm;
    efvm.SaveMeshInGMSHFile();
    delete[] points1;
    delete[] isurfaces1;
    delete[] bsurfaces1;
    delete[] hexahedrons1;
    delete[] points2;
    delete[] isurfaces2;
    delete[] bsurfaces2;
    delete[] hexahedrons2;
    system("pause");
    return 0;
}


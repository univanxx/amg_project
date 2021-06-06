// #include "stdafx.h"
#include <iostream>
using namespace std;
#include "main.h"
#include "mesh.h"
#include "Euler2.h"
#include "Vectors.h"
#include "Matrices.h"
#include <math.h>
#include <algorithm>

extern mesh::point1* points1;
extern mesh::internal_surface1* isurfaces1;
extern mesh::boundary_surface1* bsurfaces1;
extern mesh::hexahedron1* hexahedrons1;

extern mesh::point2* points2;
extern mesh::internal_surface2* isurfaces2;
extern mesh::boundary_surface2* bsurfaces2;
extern mesh::hexahedron2* hexahedrons2;

// Получение вектора примитивных переменных
Vector get_primitives(Vector& conservatives)
{
	Vector primitive(5);
	primitive[0] = conservatives[0];
	primitive[1] = conservatives[1] / primitive[0];
	primitive[2] = conservatives[2] / primitive[0];
	primitive[3] = conservatives[3] / primitive[0];
	primitive[4] = 0.4 * (conservatives[4] - primitive[0] * 0.5 * (primitive[1] * primitive[1] + primitive[2] * primitive[2] + primitive[3] * primitive[3]));
	return primitive;
}

// Получение вектора потоков
Vector get_fluxes(Vector& primitives)
{
	Vector F(5);
	F[0] = primitives[0] * primitives[1];
	F[1] = primitives[0] * primitives[1] * primitives[1] + primitives[4];
	F[2] = primitives[0] * primitives[1] * primitives[2];
	F[3] = primitives[0] * primitives[1] * primitives[3];
	F[4] = primitives[1] * (primitives[4] + 2.5 * primitives[4] + primitives[0] * 0.5 * (primitives[1] * primitives[1] + primitives[2] * primitives[2] + primitives[3] * primitives[3]));
	return F;
}


int main()
{
    mesh::LoadMesh("C:/Users/wchhi/source/repos/EulerProj/EulerProj/meshes/mymesh.txt");
	//mesh::LoadMesh("C:/Users/Asus/Documents/Visual Studio 2013/Projects/EulerProject/meshes/mymesh.txt");
	// Шаг по времени
    double tau = 0.001;
	const double epsilon = 1e-10;
	// Начальные условия
    Matrix primitive(count_elements1, 5);
    for (int i = 0; i < count_elements1; ++i)
    {
        if (points1[i].entry_or_boundary_condition == 0)
        {
            primitive.setElement(i, 0, 0.125);  // ro
            primitive.setElement(i, 1, 0.0);  // u
            primitive.setElement(i, 2, 0.0);  // v
            primitive.setElement(i, 3, 0.0);  // w
            primitive.setElement(i, 4, 0.1);  // p
        }
        else if (points1[i].entry_or_boundary_condition == 1)
        {
            primitive.setElement(i, 0, 1.0);
            primitive.setElement(i, 1, 0.0);
            primitive.setElement(i, 2, 0.0);
            primitive.setElement(i, 3, 0.0);
            primitive.setElement(i, 4, 1.0);
        }
    }

    // Граничные условия для сверхзвуковой границы входа
	Vector v_lim(3);
	v_lim[0] = 0;
	v_lim[1] = 0;
	v_lim[2] = 0;
	double ro_lim = 0.125;
	double p_lim = 0.1;

	// Решение
	Euler efvm;

	// Консервативные переменные
    Matrix U(count_elements1, 5);
	// Время - 0.25 секунды, шаг - 0.001 -> всего 250 итераций
    for (int i = 0; i < 250; ++i)
	{
		cout << endl << i + 1 << "TIME STEP" << endl;
		// Переход к консервативным переменным
		for (int j = 0; j < count_elements1; ++j)
		{
			U.setElement(j, 0, primitive.getElements()[j][0]);
			U.setElement(j, 1, primitive.getElements()[j][1] * primitive.getElements()[j][0]);
			U.setElement(j, 2, primitive.getElements()[j][2] * primitive.getElements()[j][0]);
			U.setElement(j, 3, primitive.getElements()[j][3] * primitive.getElements()[j][0]);
			U.setElement(j, 4, 2.5 * primitive.getElements()[j][4] + primitive.getElements()[j][0] * 0.5 * (primitive.getElements()[j][1] * primitive.getElements()[j][1] + primitive.getElements()[j][2] * primitive.getElements()[j][2] + primitive.getElements()[j][3] * primitive.getElements()[j][3]));
		}
		// Вектор для обновлённых значений U
        Matrix U_new = U;
		// Проходимся по внутренним элементам
		for (int j = 0; j < count_internal_surfaces1; ++j)
		{
            Matrix T(5, 5);
            if (abs(isurfaces1[j].nz - 1.0) < epsilon)
			{
				T.setElement(0, 0, 1.0);
				T.setElement(1, 3, 1.0);
				T.setElement(2, 2, 1.0);
				T.setElement(4, 4, 1.0);
				T.setElement(3, 1, -1.0);
			}
			else if (abs(isurfaces1[j].nz + 1.0) < epsilon)
			{
				T.setElement(0, 0, 1.0);
				T.setElement(1, 3, -1.0);
				T.setElement(2, 2, -1.0);
				T.setElement(4, 4, 1.0);
				T.setElement(3, 1, 1.0);
			}
            else //if (abs(isurfaces1[j].nz*isurfaces1[j].nz - 1.0) > epsilon)
			{
				T.setElement(0, 0, 1.0);
				T.setElement(1, 1, isurfaces1[j].nx);
				T.setElement(1, 2, isurfaces1[j].ny);
				T.setElement(1, 3, isurfaces1[j].nz);
				T.setElement(2, 1, -isurfaces1[j].ny / (sqrt(1.0 - isurfaces1[j].nz * isurfaces1[j].nz)));
				T.setElement(2, 2, isurfaces1[j].nx / (sqrt(1.0 - isurfaces1[j].nz * isurfaces1[j].nz)));
				T.setElement(3, 1, -isurfaces1[j].nx * isurfaces1[j].nz / (sqrt(1.0 - isurfaces1[j].nz * isurfaces1[j].nz)));
				T.setElement(3, 2, -isurfaces1[j].ny * isurfaces1[j].nz / (sqrt(1.0 - isurfaces1[j].nz * isurfaces1[j].nz)));
				T.setElement(3, 3, sqrt(1.0 - isurfaces1[j].nz * isurfaces1[j].nz));
				T.setElement(4, 4, 1.0);
			}
			// Нахождение потока по методу HLL. Переход в локальную систему координат
            Vector U_left = T * U[isurfaces1[j].element1];
            Vector U_right = T * U[isurfaces1[j].element2];

            Vector primitive_left = get_primitives(U_left);
            Vector primitive_right = get_primitives(U_right);

			// Скорость звука
			double a_left = sqrt(1.4 * (primitive_left[4] / primitive_left[0]));
			double a_right = sqrt(1.4 *(primitive_right[4] / primitive_right[0]));

			// Нормальная компонента скорости 
			double absVec_L = abs(primitive_left[1]);
			double absVec_R = abs(primitive_right[1]);
			// Скорости волн 
			double S_r = max(absVec_L + a_left, absVec_R + a_right);
			double S_l = -S_r;

            Vector F_left = get_fluxes(primitive_left);
            Vector F_right = get_fluxes(primitive_right);

            Vector F_hll(5);
			if (S_l >= 0.0)
			{
				F_hll = F_left;
			}
			else if (S_l <= 0.0 && S_r >= 0.0)
			{
				F_hll = (F_left * S_r - F_right * S_l + (U_right - U_left) * S_l * S_r) / (S_r - S_l);
			}
			else if (S_r <= 0.0)
			{
				F_hll = F_right;
			}
			// Переход из локальной системы координат в главную
			
			T = T.transpose();
            Vector F_res = T * F_hll;

			// Обноваление значений
			for (int k = 0; k < 5; ++k)
			{
				U_new.substractElement(isurfaces1[j].element1, k, (F_res[k] * (tau * isurfaces1[j].area / hexahedrons1[isurfaces1[j].element1].volume)));
				U_new.addElement(isurfaces1[j].element2, k, (F_res[k] * (tau * isurfaces1[j].area / hexahedrons1[isurfaces1[j].element2].volume)));
			}
		}
		// Проходимся по граничным граням
		for (int j = 0; j < count_boundary_surfaces1; ++j)
		{
            Vector F(5);

            Vector U_b = U[bsurfaces1[j].element1];
            Vector primitive_b = get_primitives(U_b);

			if (bsurfaces1[j].type == FreedomExit)
			{
				//cout << "сверхзвуковая граница выхода" << endl;
				double v_n = primitive_b[1] * bsurfaces1[j].nx + primitive_b[2] * bsurfaces1[j].ny + primitive_b[3] * bsurfaces1[j].nz;
				F[0] = primitive_b[0] * v_n;
				F[1] = primitive_b[0] * primitive_b[1] * v_n + primitive_b[4] * bsurfaces1[j].nx;
				F[2] = primitive_b[0] * primitive_b[2] * v_n + primitive_b[4] * bsurfaces1[j].ny;
				F[3] = primitive_b[0] * primitive_b[3] * v_n + primitive_b[4] * bsurfaces1[j].nz;
                F[4] = primitive_b[0] * v_n * (primitive_b[4] + 2.5 * primitive_b[4] + 0.5 * primitive_b[0] * (primitive_b[1] * primitive_b[1] + primitive_b[2] * primitive_b[2] + primitive_b[3] * primitive_b[3]));
			}
			else if (bsurfaces1[j].type == HardBorder)
			{
				//cout << "поверхность, представляющая собой твёрдую непроницаемую стенку" << endl;
				F[0] = 0.0;
				F[1] = primitive_b[4] * bsurfaces1[j].nx;
				F[2] = primitive_b[4] * bsurfaces1[j].ny;
				F[3] = primitive_b[4] * bsurfaces1[j].nz;
				F[4] = 0.0;
			}
			else if (bsurfaces1[j].type == Symmetry)
			{
				//cout << "плоскость симметрии" << endl;
				F[0] = 0.0;
				F[1] = primitive_b[4] * bsurfaces1[j].nx;
				F[2] = primitive_b[4] * bsurfaces1[j].ny;
				F[3] = primitive_b[4] * bsurfaces1[j].nz;
				F[4] = 0.0;
			}
			else //(bsurfaces1[j].type == Custom)
			{
				//cout << "Сверхзвуковая граница входа" << endl;
				double v_n_lim = v_lim[0] * bsurfaces1[j].nx + v_lim[1] * bsurfaces1[j].ny + v_lim[2] * bsurfaces1[j].nz;
				F[0] = ro_lim * v_n_lim;
				F[1] = ro_lim * v_lim[0] * v_n_lim + p_lim * bsurfaces1[j].nx;
				F[2] = ro_lim * v_lim[1] * v_n_lim + p_lim * bsurfaces1[j].ny;
				F[3] = ro_lim * v_lim[2] * v_n_lim + p_lim * bsurfaces1[j].nz;
				F[4] = ro_lim * v_n_lim * (p_lim + 2.5 * p_lim + 0.5 * ro_lim * (v_lim[0] * v_lim[0] + v_lim[1] * v_lim[1] + v_lim[2] * v_lim[2]));
			}

			for (int k = 0; k < 5; ++k)
			{
				U_new.substractElement(bsurfaces1[j].element1, k, (F[k] * (tau* bsurfaces1[j].area / hexahedrons1[bsurfaces1[j].element1].volume)));
			}
		}
		// Переход к примитивным элементам
		for (int j = 0; j < count_elements1; ++j)
		{
			primitive.setElement(j, 0, U_new.getElements()[j][0]);
			primitive.setElement(j, 1, U_new.getElements()[j][1] / primitive.getElements()[j][0]);
			primitive.setElement(j, 2, U_new.getElements()[j][2] / primitive.getElements()[j][0]);
			primitive.setElement(j, 3, U_new.getElements()[j][3] / primitive.getElements()[j][0]);
			primitive.setElement(j, 4, 0.4 * (U_new.getElements()[j][4] - primitive.getElements()[j][0] * 0.5 * (primitive.getElements()[j][1] * primitive.getElements()[j][1] + primitive.getElements()[j][2] * primitive.getElements()[j][2] + primitive.getElements()[j][3] * primitive.getElements()[j][3])));
		}
        if (i == 249)
		{
			cout << "SOLUTION FOUND!" << endl;
			for (int j = 0; j < count_elements1; ++j)
			{
				efvm.m_u0[j] = U_new[j][0];
				efvm.m_u1[j] = U_new[j][1];
				efvm.m_u2[j] = U_new[j][2];
				efvm.m_u3[j] = U_new[j][3];
				efvm.m_u4[j] = U_new[j][4];
			}
		}
	}
	efvm.SaveSolutionInGMSHFile();
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


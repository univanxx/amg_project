// #include "stdafx.h"
#include <iostream>
using namespace std;
#include "main.h"
#include "mesh.h"
#include "Euler2.h"
#include "Vectors.h"
#include "Matrices.h"
#include "Tensors.h"
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
Vector get_primitives(Vector conservatives)
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
Vector get_fluxes(Vector primitives)
{
	Vector F(5);
	F[0] = primitives[0] * primitives[1];
	F[1] = primitives[0] * primitives[1] * primitives[1] + primitives[4];
	F[2] = primitives[0] * primitives[1] * primitives[2];
	F[3] = primitives[0] * primitives[1] * primitives[3];
	F[4] = primitives[1] * (primitives[4] + 10.0 / 4 * primitives[4] + primitives[0] * 0.5 * (primitives[1] * primitives[1] + primitives[2] * primitives[2] + primitives[3] * primitives[3]));
	return F;
}

Vector F_hll(int j, double epsilon, Matrix U)
{
	Matrix T(5, 5);
	if (fabs(isurfaces1[j].nz - 1.0) < epsilon)
	{
		T.setElement(0, 0, 1.0);
		T.setElement(1, 3, 1.0);
		T.setElement(2, 2, 1.0);
		T.setElement(4, 4, 1.0);
		T.setElement(3, 1, -1.0);
	}
	else if (fabs(isurfaces1[j].nz + 1.0) < epsilon)
	{
		T.setElement(0, 0, 1.0);
		T.setElement(1, 3, -1.0);
		T.setElement(2, 2, -1.0);
		T.setElement(4, 4, 1.0);
		T.setElement(3, 1, 1.0);
	}
	else if (fabs(isurfaces1[j].nz*isurfaces1[j].nz - 1.0) > epsilon)
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
	double a_left = sqrt(1.4 * primitive_left[4] / primitive_left[0]);
	double a_right = sqrt(1.4 * primitive_right[4] / primitive_right[0]);
	// Нормальная компонента скорости 
	double absVec_L = fabs(primitive_left[1]);
	double absVec_R = fabs(primitive_right[1]);
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
	return F_res;
}

Matrix Jacobian(int j, Matrix U)
{
	// Консервативыне переменные элемента1, находим по ребру. Нормаль к ребру = нормаль к элементу1.
	Vector conservatives = U[isurfaces1[j].element1];
	Vector primitives = get_primitives(conservatives);
	double v_n = primitives[1] * isurfaces1[j].nx + primitives[2] * isurfaces1[j].ny + primitives[3] * isurfaces1[j].nz;
	double a = 1.4 * primitives[4] / primitives[0];
	double h = a / 0.4;
	double e_k = 0.5 * (primitives[1] * primitives[1] + primitives[2] * primitives[2] + primitives[3] * primitives[3]);
	double h_0 = h + e_k;

	Matrix A(5, 5);
	A.setElement(0, 0, 0.0);
	A.setElement(1, 0, 0.4 * e_k * isurfaces1[j].nx - primitives[1] * v_n);
	A.setElement(2, 0, 0.4 * e_k * isurfaces1[j].ny - primitives[2] * v_n);
	A.setElement(3, 0, 0.4 * e_k * isurfaces1[j].nz - primitives[3] * v_n);
	A.setElement(4, 0, v_n * (0.4 * e_k - h_0));

	A.setElement(0, 1, isurfaces1[j].nx);
	A.setElement(1, 1, v_n + 0.6 * primitives[1] * isurfaces1[j].nx);
	A.setElement(2, 1, primitives[2] * isurfaces1[j].nx - 0.4 * primitives[1] * isurfaces1[j].ny);
	A.setElement(3, 1, primitives[3] * isurfaces1[j].nx - 0.4 * primitives[1] * isurfaces1[j].nz);
	A.setElement(4, 1, h_0 * isurfaces1[j].nx - 0.4 * primitives[1] * v_n);

	A.setElement(0, 2, isurfaces1[j].ny);
	A.setElement(1, 2, primitives[1] * isurfaces1[j].ny - 0.4 * primitives[2] * isurfaces1[j].nx);
	A.setElement(2, 2, v_n + 0.6 * primitives[2] * isurfaces1[j].ny);
	A.setElement(3, 2, primitives[3] * isurfaces1[j].ny - 0.4 * primitives[2] * isurfaces1[j].nz);
	A.setElement(4, 2, h_0 * isurfaces1[j].ny - 0.4 * primitives[2] * v_n);

	A.setElement(0, 3, isurfaces1[j].nz);
	A.setElement(1, 3, primitives[1] * isurfaces1[j].nz - 0.4 * primitives[3] * isurfaces1[j].nx);
	A.setElement(2, 3, primitives[2] * isurfaces1[j].nz - 0.4 * primitives[3] * isurfaces1[j].ny);
	A.setElement(3, 3, v_n + 0.6 * primitives[3] * isurfaces1[j].nz);
	A.setElement(4, 3, h_0 * isurfaces1[j].nz - 0.4 * primitives[3] * v_n);

	A.setElement(0, 4, 0.0);
	A.setElement(1, 4, 0.4 * isurfaces1[j].nx);
	A.setElement(2, 4, 0.4 * isurfaces1[j].ny);
	A.setElement(3, 4, 0.4 * isurfaces1[j].nz);
	A.setElement(4, 4, 1.4 * v_n);

	return A;
}


int main()
{
	setlocale(LC_ALL, "Russian");
	// mesh::LoadMesh("C:/Users/wchhi/source/repos/EulerProj/EulerProj/meshes/mymesh.txt");
	mesh::LoadMesh("C:/Users/Asus/Documents/Visual Studio 2013/Projects/EulerProject/meshes/mymesh.txt");
	// Шаг по времени
	double tau = 0.05;
	const double epsilon = 1e-10;
	cout << "Начальные условия:" << endl;
	Matrix primitive(count_elements1, 5);
	for (int i = 0; i < count_elements1; ++i)
	{
		if (points1[i].entry_or_boundary_condition == 0)
		{
			primitive.setElement(i, 0, 0.1);  // ro
			primitive.setElement(i, 1, 0.5);  // u
			primitive.setElement(i, 2, 0.6);  // v 
			primitive.setElement(i, 3, 0.7);  // w
			primitive.setElement(i, 4, 1.0);  // p
		}
		else if (points1[i].entry_or_boundary_condition == 1)
		{
			primitive.setElement(i, 0, 1.0);
			primitive.setElement(i, 1, 0.4);
			primitive.setElement(i, 2, 0.6);
			primitive.setElement(i, 3, 0.2);
			primitive.setElement(i, 4, 0.125);
		}
	}

	cout << "Граничные условия для сверхзвуковой границы входа:" << endl;
	Vector v_lim(3);
	v_lim[0] = 1.5;
	v_lim[1] = 1.0;
	v_lim[2] = 0.7;
	double ro_lim = 1.5;
	double p_lim = 2.3;

	// Решение
	Euler efvm;
	// Консервативные переменные
	Matrix U(count_elements1, 5);
	// Время - 1 секунда, шаг - 0.05 -> всего 5 итераций
	for (int i = 0; i < 1; ++i)
	{
		cout << endl << i + 1 << "-Й МОМЕНТ ВРЕМЕНИ" << endl;
		cout << "Переход к консервативным переменным" << endl;
		// Переход к консервативным переменным
		for (int j = 0; j < count_elements1; ++j)
		{
			U.setElement(j, 0, primitive.getElements()[j][0]);
			U.setElement(j, 1, primitive.getElements()[j][1] * primitive.getElements()[j][0]);
			U.setElement(j, 2, primitive.getElements()[j][2] * primitive.getElements()[j][0]);
			U.setElement(j, 3, primitive.getElements()[j][3] * primitive.getElements()[j][0]);
			U.setElement(j, 4, 2.5 * primitive.getElements()[j][4] + primitive.getElements()[j][0] * 0.5 * (primitive.getElements()[j][1] * primitive.getElements()[j][1] + primitive.getElements()[j][2] * primitive.getElements()[j][2] + primitive.getElements()[j][3] * primitive.getElements()[j][3]));
		}
		// Матрица системы - для начала 0 и 1 для элементов и их соседей
		Matrix A(count_elements1, count_elements1);
		// Вектор правой части СЛАУ
		Matrix B(count_elements1, 5);
		for (int j = 0; j < count_internal_surfaces1; ++j)
		{
			A.setElement(isurfaces1[j].element1, isurfaces1[j].element2, 0.5);
			B[isurfaces1[j].element1] = F_hll(j, epsilon, U);
				// Надо установить фиктивный элемент для граничных рёбер
		}
		double **res = A.getElements();
		for (int row = 0; row < count_elements1; ++row)
		{
			for (int col = 0; col < count_elements1; ++col)
			{
				cout << res[row][col] << " ";
			}
			cout << endl;
		}
	}
	/*
	Matrix A(3, 3);
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			A.setElement(row, col, 1.5);
		}
		cout << endl;
	}
	A.setElement(1, 2, 3);
	A.setElement(2, 1, 2);
	cout << A << endl;

	Tensor B(2, 3, 3);
	for (int row = 0; row < 3; ++row)
	{
		for (int col = 0; col < 3; ++col)
		{
			B.setElement(0, col, col, 0.5);
		}
		cout << endl;
	}
	B.setElement(0, 1, 2, 1);
	B.setElement(0, 2, 1, -4);
	cout << B[0] << endl;
	Matrix R = A * B[0];
	cout << R << endl;
		/*

		// Вектор для обновлённых значений U
		Matrix U_new = U;

		cout << "Проходимся по внутренним граням -> по сумме рёбер" << endl;
		// Проходимся по внутренним элементам
		for (int j = 0; j < count_internal_surfaces1; ++j)
		{
			Matrix T(5, 5);
			//cout << isurfaces1[j].nx << isurfaces1[j].ny << isurfaces1[j].nz << endl;
			if (fabs(isurfaces1[j].nz - 1.0) < epsilon)
			{
				T.setElement(0, 0, 1.0);
				T.setElement(1, 3, 1.0);
				T.setElement(2, 2, 1.0);
				T.setElement(4, 4, 1.0);
				T.setElement(3, 1, -1.0);
			}
			else if (fabs(isurfaces1[j].nz + 1.0) < epsilon)
			{
				T.setElement(0, 0, 1.0);
				T.setElement(1, 3, -1.0);
				T.setElement(2, 2, -1.0);
				T.setElement(4, 4, 1.0);
				T.setElement(3, 1, 1.0);
			}
			else if (fabs(isurfaces1[j].nz*isurfaces1[j].nz - 1.0) > epsilon)
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
			//cout << T << endl;
			// Нахождение потока по методу HLL. Переход в локальную систему координат
			Vector U_left = T * U[isurfaces1[j].element1];
			Vector U_right = T * U[isurfaces1[j].element2];
			//cout << T << endl;
			Vector primitive_left = get_primitives(U_left);
			Vector primitive_right = get_primitives(U_right);
			//cout << primitive_left << endl;
			//cout << primitive_right << endl;
			//cout << U_left << endl;
			//cout << U_right << endl;
			// Скорость звука
			double a_left = sqrt(1.4 * primitive_left[4] / primitive_left[0]);
			double a_right = sqrt(1.4 * primitive_right[4] / primitive_right[0]);
			//cout << a_left << " " << a_right << endl;
			// Нормальная компонента скорости 
			double absVec_L = fabs(primitive_left[1]);
			double absVec_R = fabs(primitive_right[1]);
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
			//cout << F_res << endl;
			// Обноваление значений
			for (int k = 0; k < 5; ++k)
			{
				U_new.getElements()[isurfaces1[j].element1][k] -= (F_res[k] * (tau * isurfaces1[isurfaces1[j].element1].area / hexahedrons1[isurfaces1[j].element1].volume));
				U_new.getElements()[isurfaces1[j].element2][k] += (F_res[k] * (tau * isurfaces1[isurfaces1[j].element2].area / hexahedrons1[isurfaces1[j].element2].volume));
			}
		}
		// Проходимся по граничным элементам
		cout << "Проходимся по граничным граням" << endl;
		for (int j = 0; j < count_boundary_surfaces1; ++j)
		{
			Vector F_hll(5);
			Matrix T(5, 5);
			if (fabs(bsurfaces1[j].nz - 1.0) < epsilon)
			{
				T.setElement(0, 0, 1.0);
				T.setElement(1, 3, 1.0);
				T.setElement(2, 2, 1.0);
				T.setElement(4, 4, 1.0);
				T.setElement(3, 1, -1.0);
			}
			else if (fabs(bsurfaces1[j].nz + 1.0) < epsilon)
			{
				T.setElement(0, 0, 1.0);
				T.setElement(1, 3, -1.0);
				T.setElement(2, 2, -1.0);
				T.setElement(4, 4, 1.0);
				T.setElement(3, 1, 1.0);
			}
			else if (fabs(bsurfaces1[j].nz * bsurfaces1[j].nz - 1.0) > epsilon)
			{
				T.setElement(0, 0, 1.0);
				T.setElement(1, 1, bsurfaces1[j].nx);
				T.setElement(1, 2, bsurfaces1[j].ny);
				T.setElement(1, 3, bsurfaces1[j].nz);
				T.setElement(2, 1, -bsurfaces1[j].ny / (sqrt(1.0 - bsurfaces1[j].nz * bsurfaces1[j].nz)));
				T.setElement(2, 2, bsurfaces1[j].nx / (sqrt(1.0 - bsurfaces1[j].nz * bsurfaces1[j].nz)));
				T.setElement(3, 1, -bsurfaces1[j].nx * bsurfaces1[j].nz / (sqrt(1.0 - bsurfaces1[j].nz * bsurfaces1[j].nz)));
				T.setElement(3, 2, -bsurfaces1[j].ny * bsurfaces1[j].nz / (sqrt(1.0 - bsurfaces1[j].nz * bsurfaces1[j].nz)));
				T.setElement(3, 3, sqrt(1.0 - bsurfaces1[j].nz * bsurfaces1[j].nz));
				T.setElement(4, 4, 1.0);
			}
			//cout << bsurfaces1[j].nz << endl;
			//cout << T << endl;
			// Переход в локальную систему координат
			Vector U_rotated = T * U[bsurfaces1[j].element1];
			Vector primitive_rotated = get_primitives(U_rotated);
			//cout << U_rotated << endl;
			if (bsurfaces1[j].type == FreedomExit)
			{
				//cout << "сверхзвуковая граница выхода" << endl;
				double v_n = primitive_rotated[1] * bsurfaces1[j].nx + primitive_rotated[2] * bsurfaces1[j].ny + primitive_rotated[3] * bsurfaces1[j].nz;
				F_hll[0] = primitive_rotated[0] * v_n;
				F_hll[1] = primitive_rotated[0] * primitive_rotated[1] * v_n + primitive_rotated[4] * bsurfaces1[j].nx;
				F_hll[2] = primitive_rotated[0] * primitive_rotated[2] * v_n + primitive_rotated[4] * bsurfaces1[j].ny;
				F_hll[3] = primitive_rotated[0] * primitive_rotated[3] * v_n + primitive_rotated[4] * bsurfaces1[j].nz;
				F_hll[4] = primitive_rotated[0] * v_n * (U[bsurfaces1[j].element1][4] + primitive_rotated[4]);
			}
			else if (bsurfaces1[j].type == HardBorder)
			{
				//cout << "поверхность, представляющая собой твёрдую непроницаемую стенку" << endl;
				F_hll[0] = 0.0;
				F_hll[1] = primitive_rotated[4] * bsurfaces1[j].nx;
				F_hll[2] = primitive_rotated[4] * bsurfaces1[j].ny;
				F_hll[3] = primitive_rotated[4] * bsurfaces1[j].nz;
				F_hll[4] = 0.0;
			}
			else if (bsurfaces1[j].type == Symmetry)
			{
				//cout << "плоскость симметрии" << endl;
				F_hll[0] = 0.0;
				F_hll[1] = primitive_rotated[4] * bsurfaces1[j].nx;
				F_hll[2] = primitive_rotated[4] * bsurfaces1[j].ny;
				F_hll[3] = primitive_rotated[4] * bsurfaces1[j].nz;
				F_hll[4] = 0.0;
			}
			else //(bsurfaces1[j].type == Custom)
			{
				//cout << "Сверхзвуковая граница входа" << endl;
				double v_n_lim = v_lim[0] * bsurfaces1[j].nx + v_lim[1] * bsurfaces1[j].ny + v_lim[2] * bsurfaces1[j].nz;
				F_hll[0] = ro_lim * v_n_lim;
				F_hll[1] = ro_lim * v_lim[0] * v_n_lim + p_lim * bsurfaces1[j].nx;
				F_hll[2] = ro_lim * v_lim[1] * v_n_lim + p_lim * bsurfaces1[j].ny;
				F_hll[3] = ro_lim * v_lim[2] * v_n_lim + p_lim * bsurfaces1[j].nz;
				F_hll[4] = ro_lim * v_n_lim * (p_lim + 10.0 / 4 * p_lim + 0.5 * ro_lim * (v_lim[0] * v_lim[0] + v_lim[1] * v_lim[1] + v_lim[2] * v_lim[2]));
			}
			T = T.transpose();
			Vector F_res = T * F_hll;
			//cout << F_res << endl;
			for (int k = 0; k < 5; ++k)
			{
				U_new.getElements()[bsurfaces1[j].element1][k] -= (F_res[k] * (tau * bsurfaces1[bsurfaces1[j].element1].area / hexahedrons1[bsurfaces1[j].element1].volume));
			}
		}

		cout << "Переход к примитивным переменным" << endl;
		// Переход к примитивным элементам
		for (int j = 0; j < count_elements1; ++j)
		{
			primitive.getElements()[j][0] = U_new[j][0];
			primitive.getElements()[j][1] = U_new[j][1] / primitive.getElements()[j][0];
			primitive.getElements()[j][2] = U_new[j][2] / primitive.getElements()[j][0];
			primitive.getElements()[j][3] = U_new[j][3] / primitive.getElements()[j][0];
			primitive.getElements()[j][4] = 0.4 * (U_new[j][4] - primitive.getElements()[j][0] * 0.5 *(primitive.getElements()[j][1] * primitive.getElements()[j][1] + primitive.getElements()[j][2] * primitive.getElements()[j][2] + primitive.getElements()[j][3] * primitive.getElements()[j][3]));
			//cout << U.getElements()[j][0] << U.getElements()[j][1] << U.getElements()[j][2] << U.getElements()[j][3] << U.getElements()[j][4] << endl;
			//cout << U_new.getElements()[j][0] << U_new.getElements()[j][1] << U_new.getElements()[j][2] << U_new.getElements()[j][3] << U_new.getElements()[j][4] << endl;
			//cout << primitive.getElements()[j][0] << primitive.getElements()[j][1] << primitive.getElements()[j][2] << primitive.getElements()[j][3] << primitive.getElements()[j][4] << endl;
			//cout << "--------------------------------------------------------------------------" << endl;
		}
		if (i == 4)
		{
			cout << "Решение найдено!" << endl;
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
	*/
	system("pause");
	return 0;
}
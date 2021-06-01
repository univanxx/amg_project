// #include "stdafx.h"
#include <iostream>
using namespace std;
#include "main.h"
#include "mesh.h"
#include "Euler2.h"
#include "Vectors.h"
#include "Matrices.h"
#include "Tensors.h"
#include "Tensor4.h"
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

Vector F_hll(int j, double epsilon, Matrix &U)
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

Matrix R(int j, const Matrix &U, double epsilon)
{
	// Консервативыне переменные элемента1, находим по ребру. Нормаль к ребру = нормаль к элементу1.
	Vector conservatives(5);
	conservatives = U.getElements()[isurfaces1[j].element1];
	Vector primitives = get_primitives(conservatives);
	double v_n = primitives[1] * isurfaces1[j].nx + primitives[2] * isurfaces1[j].ny + primitives[3] * isurfaces1[j].nz;
	double a = 1.4 * primitives[4] / primitives[0];
	double h = a / 0.4;
	double e_k = 0.5 * (primitives[1] * primitives[1] + primitives[2] * primitives[2] + primitives[3] * primitives[3]);
	double h_0 = h + e_k;
	// Правые собственные вектора Якобиана
	Matrix R(5, 5);
	if (abs(abs(isurfaces1[j].nx) - 1.0) < epsilon)
	{
		R.setElement(0, 0, 1);
		R.setElement(1, 0, primitives[1] - a * isurfaces1[j].nx);
		R.setElement(2, 0, primitives[2] - a * isurfaces1[j].ny);
		R.setElement(3, 0, primitives[3] - a * isurfaces1[j].nz);
		R.setElement(4, 0, h_0 - a * v_n);

		R.setElement(0, 1, 1);
		R.setElement(1, 1, primitives[1]);
		R.setElement(2, 1, primitives[2]);
		R.setElement(3, 1, primitives[3]);
		R.setElement(4, 1, e_k);

		R.setElement(0, 2, 1);
		R.setElement(1, 2, primitives[1] + a * isurfaces1[j].nx);
		R.setElement(2, 2, primitives[2] + a * isurfaces1[j].ny);
		R.setElement(3, 2, primitives[3] + a * isurfaces1[j].nz);
		R.setElement(4, 2, h_0 + a * v_n);

		R.setElement(0, 3, 0);
		R.setElement(1, 3, isurfaces1[j].ny);
		R.setElement(2, 3, -isurfaces1[j].nx);
		R.setElement(3, 3, 0);
		R.setElement(4, 3, primitives[1] * isurfaces1[j].ny - primitives[2] * isurfaces1[j].nx);

		R.setElement(0, 4, 0);
		R.setElement(1, 4, -isurfaces1[j].nz);
		R.setElement(2, 4, 0);
		R.setElement(3, 4, isurfaces1[j].nx);
		R.setElement(4, 4, primitives[3] * isurfaces1[j].nx - primitives[1] * isurfaces1[j].nz);
		return R;
	}
	else if (abs(abs(isurfaces1[j].ny) - 1.0) < epsilon)
	{
		R.setElement(0, 0, 1);
		R.setElement(1, 0, primitives[1] - a * isurfaces1[j].nx);
		R.setElement(2, 0, primitives[2] - a * isurfaces1[j].ny);
		R.setElement(3, 0, primitives[3] - a * isurfaces1[j].nz);
		R.setElement(4, 0, h_0 - a * v_n);

		R.setElement(0, 1, 1);
		R.setElement(1, 1, primitives[1]);
		R.setElement(2, 1, primitives[2]);
		R.setElement(3, 1, primitives[3]);
		R.setElement(4, 1, e_k);

		R.setElement(0, 2, 1);
		R.setElement(1, 2, primitives[1] + a * isurfaces1[j].nx);
		R.setElement(2, 2, primitives[2] + a * isurfaces1[j].ny);
		R.setElement(3, 2, primitives[3] + a * isurfaces1[j].nz);
		R.setElement(4, 2, h_0 + a * v_n);

		R.setElement(0, 3, 0);
		R.setElement(1, 3, isurfaces1[j].ny);
		R.setElement(2, 3, -isurfaces1[j].nx);
		R.setElement(3, 3, 0);
		R.setElement(4, 3, primitives[1] * isurfaces1[j].ny - primitives[2] * isurfaces1[j].nx);

		R.setElement(0, 4, 0);
		R.setElement(1, 4, 0);
		R.setElement(2, 4, isurfaces1[j].nz);
		R.setElement(3, 4, -isurfaces1[j].ny);
		R.setElement(4, 4, primitives[2] * isurfaces1[j].nz - primitives[3] * isurfaces1[j].ny);
		return R;
	}
	else if (abs(abs(isurfaces1[j].nz) - 1.0) < epsilon)
	{
		R.setElement(0, 0, 1);
		R.setElement(1, 0, primitives[1] - a * isurfaces1[j].nx);
		R.setElement(2, 0, primitives[2] - a * isurfaces1[j].ny);
		R.setElement(3, 0, primitives[3] - a * isurfaces1[j].nz);
		R.setElement(4, 0, h_0 - a * v_n);

		R.setElement(0, 1, 1);
		R.setElement(1, 1, primitives[1]);
		R.setElement(2, 1, primitives[2]);
		R.setElement(3, 1, primitives[3]);
		R.setElement(4, 1, e_k);

		R.setElement(0, 2, 1);
		R.setElement(1, 2, primitives[1] + a * isurfaces1[j].nx);
		R.setElement(2, 2, primitives[2] + a * isurfaces1[j].ny);
		R.setElement(3, 2, primitives[3] + a * isurfaces1[j].nz);
		R.setElement(4, 2, h_0 + a * v_n);

		R.setElement(0, 3, 0);
		R.setElement(1, 3, -isurfaces1[j].nz);
		R.setElement(2, 3, 0);
		R.setElement(3, 3, isurfaces1[j].nx);
		R.setElement(4, 3, primitives[3] * isurfaces1[j].nx - primitives[1] * isurfaces1[j].nz);

		R.setElement(0, 4, 0);
		R.setElement(1, 4, 0);
		R.setElement(2, 4, isurfaces1[j].nz);
		R.setElement(3, 4, -isurfaces1[j].ny);
		R.setElement(4, 4, primitives[2] * isurfaces1[j].nz - primitives[3] * isurfaces1[j].ny);
		return R;
	}

}

Matrix L(int j, const Matrix &U, double epsilon)
{
	// Консервативыне переменные элемента1, находим по ребру. Нормаль к ребру = нормаль к элементу1.
	Vector conservatives(5);
	conservatives = U.getElements()[isurfaces1[j].element1];
	Vector primitives = get_primitives(conservatives);
	double v_n = primitives[1] * isurfaces1[j].nx + primitives[2] * isurfaces1[j].ny + primitives[3] * isurfaces1[j].nz;
	double a = 1.4 * primitives[4] / primitives[0];
	double h = a / 0.4;
	double e_k = 0.5 * (primitives[1] * primitives[1] + primitives[2] * primitives[2] + primitives[3] * primitives[3]);
	double h_0 = h + e_k;
	// Левые собственные вектора Якобиана
	Matrix L(5, 5);
	if (abs(abs(isurfaces1[j].nx) - 1.0) < epsilon)
	{
		L.setElement(0, 0, (0.4 * e_k + a * v_n) / (2 * a * a));
		L.setElement(0, 1, (-0.4 * primitives[1] - a * isurfaces1[j].nx) / (2 * a * a));
		L.setElement(0, 2, (-0.4 * primitives[2] - a * isurfaces1[j].ny) / (2 * a * a));
		L.setElement(0, 3, (-0.4 * primitives[3] - a * isurfaces1[j].nz) / (2 * a * a));
		L.setElement(0, 4, 0.4 / (2 * a * a));

		L.setElement(1, 0, (a * a - 0.4 * e_k) / (a * a));
		L.setElement(1, 1, (0.4 * primitives[1]) / (a * a));
		L.setElement(1, 2, (0.4 * primitives[2]) / (a * a));
		L.setElement(1, 3, (0.4 * primitives[3]) / (a * a));
		L.setElement(1, 4, 0.4 / (a * a));

		L.setElement(2, 0, (0.4 * e_k - a * v_n) / (2 * a * a));
		L.setElement(2, 1, (-0.4 * primitives[1] + a * isurfaces1[j].nx) / (2 * a * a));
		L.setElement(2, 2, (-0.4 * primitives[2] + a * isurfaces1[j].ny) / (2 * a * a));
		L.setElement(2, 3, (-0.4 * primitives[3] + a * isurfaces1[j].nz) / (2 * a * a));
		L.setElement(2, 4, 0.4 / (2 * a * a));

		L.setElement(3, 0, (primitives[2] - v_n * isurfaces1[j].ny) / isurfaces1[j].nx);
		L.setElement(3, 1, isurfaces1[j].ny);
		L.setElement(3, 2, (isurfaces1[j].ny * isurfaces1[j].ny - 1) / isurfaces1[j].nx);
		L.setElement(3, 3, (isurfaces1[j].ny * isurfaces1[j].nz) / isurfaces1[j].nx);
		L.setElement(3, 4, 0);

		L.setElement(4, 0, (-primitives[3] + v_n * isurfaces1[j].nz) / isurfaces1[j].nx);
		L.setElement(4, 1, -isurfaces1[j].nz);
		L.setElement(4, 2, (-isurfaces1[j].ny * isurfaces1[j].nz) / isurfaces1[j].nx);
		L.setElement(4, 3, (1 - isurfaces1[j].nz * isurfaces1[j].nz) / isurfaces1[j].nx);
		L.setElement(4, 4, 0);
		return L;
	}
	else if (abs(abs(isurfaces1[j].ny) - 1.0) < epsilon)
	{
		L.setElement(0, 0, (0.4 * e_k + a * v_n) / (2 * a * a));
		L.setElement(0, 1, (-0.4 * primitives[1] - a * isurfaces1[j].nx) / (2 * a * a));
		L.setElement(0, 2, (-0.4 * primitives[2] - a * isurfaces1[j].ny) / (2 * a * a));
		L.setElement(0, 3, (-0.4 * primitives[3] - a * isurfaces1[j].nz) / (2 * a * a));
		L.setElement(0, 4, 0.4 / (2 * a * a));

		L.setElement(1, 0, (a * a - 0.4 * e_k) / (a * a));
		L.setElement(1, 1, (0.4 * primitives[1]) / (a * a));
		L.setElement(1, 2, (0.4 * primitives[2]) / (a * a));
		L.setElement(1, 3, (0.4 * primitives[3]) / (a * a));
		L.setElement(1, 4, 0.4 / (a * a));

		L.setElement(2, 0, (0.4 * e_k - a * v_n) / (2 * a * a));
		L.setElement(2, 1, (-0.4 * primitives[1] + a * isurfaces1[j].nx) / (2 * a * a));
		L.setElement(2, 2, (-0.4 * primitives[2] + a * isurfaces1[j].ny) / (2 * a * a));
		L.setElement(2, 3, (-0.4 * primitives[3] + a * isurfaces1[j].nz) / (2 * a * a));
		L.setElement(2, 4, 0.4 / (2 * a * a));

		L.setElement(3, 0, (-primitives[1] + v_n * isurfaces1[j].nx) / isurfaces1[j].ny);
		L.setElement(3, 1, (-isurfaces1[j].nx * isurfaces1[j].nx + 1) / isurfaces1[j].ny);
		L.setElement(3, 2, -isurfaces1[j].nx);
		L.setElement(3, 3, (-isurfaces1[j].nx * isurfaces1[j].nz) / isurfaces1[j].ny);
		L.setElement(3, 4, 0);

		L.setElement(4, 0, (primitives[3] - v_n * isurfaces1[j].nz) / isurfaces1[j].ny);
		L.setElement(4, 1, (isurfaces1[j].nx * isurfaces1[j].nz) / isurfaces1[j].ny);
		L.setElement(4, 2, isurfaces1[j].nz);
		L.setElement(4, 3, (-1 + isurfaces1[j].nz * isurfaces1[j].nz) / isurfaces1[j].ny);
		L.setElement(4, 4, 0);
		return L;
	}
	else if (abs(abs(isurfaces1[j].nz) - 1.0) < epsilon)
	{
		L.setElement(0, 0, (0.4 * e_k + a * v_n) / (2 * a * a));
		L.setElement(0, 1, (-0.4 * primitives[1] - a * isurfaces1[j].nx) / (2 * a * a));
		L.setElement(0, 2, (-0.4 * primitives[2] - a * isurfaces1[j].ny) / (2 * a * a));
		L.setElement(0, 3, (-0.4 * primitives[3] - a * isurfaces1[j].nz) / (2 * a * a));
		L.setElement(0, 4, 0.4 / (2 * a * a));

		L.setElement(1, 0, (a * a - 0.4 * e_k) / (a * a));
		L.setElement(1, 1, (0.4 * primitives[1]) / (a * a));
		L.setElement(1, 2, (0.4 * primitives[2]) / (a * a));
		L.setElement(1, 3, (0.4 * primitives[3]) / (a * a));
		L.setElement(1, 4, 0.4 / (a * a));

		L.setElement(2, 0, (0.4 * e_k - a * v_n) / (2 * a * a));
		L.setElement(2, 1, (-0.4 * primitives[1] + a * isurfaces1[j].nx) / (2 * a * a));
		L.setElement(2, 2, (-0.4 * primitives[2] + a * isurfaces1[j].ny) / (2 * a * a));
		L.setElement(2, 3, (-0.4 * primitives[3] + a * isurfaces1[j].nz) / (2 * a * a));
		L.setElement(2, 4, 0.4 / (2 * a * a));

		L.setElement(3, 0, (primitives[1] - v_n * isurfaces1[j].nx) / isurfaces1[j].nz);
		L.setElement(3, 1, isurfaces1[j].nx * isurfaces1[j].nx - 1);
		L.setElement(3, 2, (isurfaces1[j].nx * isurfaces1[j].ny) / isurfaces1[j].nz);
		L.setElement(3, 3, isurfaces1[j].nx);
		L.setElement(3, 4, 0);

		L.setElement(4, 0, (primitives[2] - v_n * isurfaces1[j].ny) / isurfaces1[j].nz);
		L.setElement(4, 1, - isurfaces1[j].nx * isurfaces1[j].ny / isurfaces1[j].nz);
		L.setElement(4, 2, (1 -isurfaces1[j].ny * isurfaces1[j].ny) / isurfaces1[j].nz);
		L.setElement(4, 3, -isurfaces1[j].ny);
		L.setElement(4, 4, 0);
		return L;
	}
	
}

Matrix I_m(int j, const Matrix &U)
{
	// Консервативыне переменные элемента1, находим по ребру. Нормаль к ребру = нормаль к элементу1.
	Vector conservatives(5);
	conservatives = U.getElements()[isurfaces1[j].element1];
	Vector primitives = get_primitives(conservatives);
	double v_n = primitives[1] * isurfaces1[j].nx + primitives[2] * isurfaces1[j].ny + primitives[3] * isurfaces1[j].nz;
	double a = 1.4 * primitives[4] / primitives[0];
	double h = a / 0.4;
	double e_k = 0.5 * (primitives[1] * primitives[1] + primitives[2] * primitives[2] + primitives[3] * primitives[3]);
	double h_0 = h + e_k;
	// Диагональная матрица 
	Matrix I(5, 5);
	I.setElement(0, 0, v_n - a);
	I.setElement(1, 1, v_n);
	I.setElement(2, 2, v_n + a);
	I.setElement(3, 3, v_n);
	I.setElement(4, 4, v_n);
	return I;
}

Matrix m_plus(Matrix M)
{
	for (int i = 0; i < M.getColumns(); ++i)
	{
		M.setElement(i, i, 0.5 * (M.getElements()[i][i] + abs(M.getElements()[i][i])));
	}
	return M;
}

Matrix m_minus(Matrix M)
{
	for (int i = 0; i < M.getColumns(); ++i)
	{
		M.setElement(i, i, 0.5 * (M.getElements()[i][i] - abs(M.getElements()[i][i])));
	}
	return M;
}

// Якобиан целиком
/*
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
*/

int main()
{
	setlocale(LC_ALL, "Russian");
	mesh::LoadMesh("C:/Users/wchhi/source/repos/EulerProj/EulerProj/meshes/mymesh.txt");
	//mesh::LoadMesh("C:/Users/Asus/Documents/Visual Studio 2013/Projects/EulerProject/meshes/mymesh.txt");
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
		// Матрица системы - блочная
		Tensor4 A(count_elements1, count_elements1, 5, 5);
		// Вектор правой части СЛАУ - блочный
		Matrix B(count_elements1, 5);
		// Проход по внутренним элементам
		for (int j = 0; j < count_internal_surfaces1; ++j)
		{
			Matrix I(5, 5);
			for (int k = 0; k < 5; ++k)
			{
				I.setElement(k, k, 1);
			}
			
			A.addMatrix(isurfaces1[j].element1, isurfaces1[j].element1, I + (R(j, U, epsilon) * m_plus(I_m(j, U)) * L(j, U, epsilon)) * (tau * isurfaces1[j].area / hexahedrons1[isurfaces1[j].element1].volume));
			A.addMatrix(isurfaces1[j].element1, isurfaces1[j].element2, (R(j, U, epsilon) * m_minus(I_m(j, U)) * L(j, U, epsilon)) * (tau * isurfaces1[j].area / hexahedrons1[isurfaces1[j].element1].volume));

			B.addVector(isurfaces1[j].element1, (R(j, U, epsilon) * m_plus(I_m(j, U)) * L(j, U, epsilon) * U[isurfaces1[j].element1] + R(j, U, epsilon) * m_minus(I_m(j, U)) * L(j, U, epsilon) * U[isurfaces1[j].element2] - F_hll(j, epsilon, U))
				* (tau * isurfaces1[j].area / hexahedrons1[isurfaces1[j].element1].volume) + U[isurfaces1[j].element1]);
		}
		// Проход по внешним элементам
		for (int j = 0; j < count_boundary_surfaces1; ++j)
		{
			Matrix I(5, 5);
			for (int k = 0; k < 5; ++k)
			{
				I.setElement(k, k, 1);
			}
			A.addMatrix(bsurfaces1[j].element1, bsurfaces1[j].element1, I + (R(j, U, epsilon) * m_plus(I_m(j, U)) * L(j, U, epsilon)) * (tau * bsurfaces1[j].area / hexahedrons1[bsurfaces1[j].element1].volume));
			B.addVector(bsurfaces1[j].element1, (R(j, U, epsilon) * m_plus(I_m(j, U)) * L(j, U, epsilon) * U[bsurfaces1[j].element1] - F_hll(j, epsilon, U)) * (tau * bsurfaces1[j].area / hexahedrons1[bsurfaces1[j].element1].volume)
				+ U[bsurfaces1[j].element1]);
		}
		// Надо установить фиктивный элемент для граничных рёбер

	}
	

	// efvm.SaveSolutionInGMSHFile();
	
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


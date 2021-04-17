#include "stdafx.h"
#include <iostream>
using namespace std;
#include "main.h"
#include "mesh.h"
#include "Euler2.h"
#include "vector"
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

// Функция умножения вектора на скаляр
vector <double> skalVec(vector < double> a, double b)
{
	for (int i = 0; i < a.size(); ++i)
	{
		a[i] = a[i] * b;
	}
	return a;
}

// Функция деления вектора на скаляр
vector <double> delVec(vector < double> a, double b)
{
	for (int i = 0; i < a.size(); ++i)
	{
		a[i] = a[i] / b;
	}
	return a;
}

// Функция сложения векторов
vector <double> plusVec(vector < double> a, vector <double> b)
{
	vector<double> c(a.size());
	for (int i = 0; i < a.size(); ++i)
	{
		c[i] += a[i] + b[i];
	}
	return c;
}

// Функция вычитания векторов
vector <double> minusVec(vector < double> a, vector <double> b)
{
	vector<double> c(a.size());
	for (int i = 0; i < a.size(); ++i)
	{
		c[i] += a[i] - b[i];
	}
	return c;
}

// Функция транспонирования матрицы
vector < vector <double> > transp(vector < vector <double> > a, int n, int m)
{
	int i, j;
	vector<vector<double> > arr(m, vector<double>(n));
	for (i = 0; i< n; i++)
		for (j = 0; j< m; j++)
			arr[j][i] = a[i][j];
	return arr;
}

// Функция умножения матрицы на вектор
vector <double> multiVec(vector < vector <double> > a, vector <double> b)
{
	vector<double> c(a.size());
	for (int i = 0; i < a.size(); ++i)
	{
		for (int j = 0; j < b.size(); ++j)
		{
			c[i] += a[i][j] * b[j];
		}
	}
	return c;
}

int main()
{
	setlocale(LC_ALL, "Russian");
	mesh::LoadMesh("D:/Загры/mymesh.txt");
	cout << "Информация о внутренних рёбрах:" << endl;
	for (int i = 0; i < 10; ++i)
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
	cout << "Граничные условия:" << endl;
	vector<double> v_lim(3);
	v_lim[0] = 3;
	v_lim[1] = 4;
	v_lim[2] = 5;
	double ro_lim = 1;
	double p_lim = 1.5;
	cout << "Скорость звука, начальная:" << endl;
	vector <double> a(count_points1);
	for (int j = 0; j < count_points1; ++j)
	{
		a[j] = sqrt(1.4 * primitive[j][4] / primitive[j][0]);
	}
	// Поток F
	cout << "Поток F, начальный:" << endl;
	vector<vector<double> > F(count_points1, vector<double>(5));
	for (int j = 0; j < count_points1; ++j)
	{
		F[j][0] = primitive[j][0] * primitive[j][1];
		F[j][1] = primitive[j][0] * primitive[j][1] * primitive[j][1] + primitive[j][4];
		F[j][0] = primitive[j][0] * primitive[j][1] * primitive[j][2];
		F[j][0] = primitive[j][0] * primitive[j][1] * primitive[j][3];
		F[j][0] = primitive[j][1] * (primitive[j][4] + primitive[j][0] * 10 / 4 * primitive[j][4] + primitive[j][0] * 0.5*(primitive[j][1] * primitive[j][1] + primitive[j][2] * primitive[j][2] + primitive[j][3] * primitive[j][3]));

	}
	// Решение
	//vector<vector<vector<double>> > U(10, vector<vector<double>>(count_points1, vector<double>(5)));
	vector<vector<double> > U(count_points1, vector<double>(5));
	// Нахождение 
	// Время - 1 секунда, шаг - 0.1
	for (int i = 0; i < 10; ++i)
	{
		cout << endl << i << "-Й МОМЕНТ ВРЕМЕНИ" << endl;
		cout << "Переход к консервативным переменным" << endl;
		// 1. Переход к консервативным переменным
		for (int j = 0; j < count_points1; ++j)
		{
			U[j][0] = primitive[j][0];
			U[j][1] = primitive[j][1] * U[j][0];
			U[j][2] = primitive[j][2] * U[j][0];
			U[j][3] = primitive[j][3] * U[j][0];
			U[j][4] = 10 / 4 * primitive[j][4] + primitive[j][0] * 0.5*(primitive[j][1] * primitive[j][1] + primitive[j][2] * primitive[j][2] + primitive[j][3] * primitive[j][3]);
		}
		cout << "Проходимся по внутренним граням" << endl;
		// Проходимся по внутренним элементам
		for (int j = 0; j < count_internal_surfaces1; ++j)
		{
			vector<double> F_hll(5);
			vector<vector<double> > T(5, vector<double>(5));
			vector<double> F_res(5);
			if (isurfaces1[j].nz == 1)
			{
				T[0][0] = 1;
				T[1][3] = 1;
				T[2][2] = 1;
				T[4][4] = 1;
				T[3][1] = -1;
			}
			else if (isurfaces1[j].nz == -1)
			{
				T[0][0] = 1;
				T[1][3] = -1;
				T[2][2] = -1;
				T[4][4] = 1;
				T[3][1] = 1;
			}
			else if (isurfaces1[j].nz*isurfaces1[j].nz != 1)
			{
				T[0][0] = 1;
				T[1][1] = isurfaces1[j].nx;
				T[1][2] = isurfaces1[j].ny;
				T[1][3] = isurfaces1[j].nz;
				T[2][1] = -isurfaces1[j].ny / (sqrt(1 - isurfaces1[j].nz*isurfaces1[j].nz));
				T[2][2] = isurfaces1[j].nx / (sqrt(1 - isurfaces1[j].nz*isurfaces1[j].nz));
				T[3][1] = -isurfaces1[j].nx*isurfaces1[j].nz / (sqrt(1 - isurfaces1[j].nz*isurfaces1[j].nz));
				T[3][2] = -isurfaces1[j].ny*isurfaces1[j].nz / (sqrt(1 - isurfaces1[j].nz*isurfaces1[j].nz));
				T[3][3] = sqrt(1 - isurfaces1[j].nz*isurfaces1[j].nz);
				T[4][4] = 1;
			}
			T = transp(T, 5, 5);
			// ... Рассчитываем F
			double absVec_L = sqrt(primitive[isurfaces1[j].element1][1] * primitive[isurfaces1[j].element1][1] + \
				primitive[isurfaces1[j].element1][2] * primitive[isurfaces1[j].element1][2] + \
				primitive[isurfaces1[j].element1][3] * primitive[isurfaces1[j].element1][3]);
			double absVec_R = sqrt(primitive[isurfaces1[j].element2][1] * primitive[isurfaces1[j].element2][1] + \
				primitive[isurfaces1[j].element2][2] * primitive[isurfaces1[j].element2][2] + \
				primitive[isurfaces1[j].element2][3] * primitive[isurfaces1[j].element2][3]);
			double S_r = max(absVec_L + a[isurfaces1[j].element1], absVec_R + a[isurfaces1[j].element1]);
			double S_l = -S_r;
			if (S_l >= 0)
			{
				F_hll = F[isurfaces1[j].element1];
			}
			else if (S_l <= 0 && S_r >= 0)
			{
				F_hll = delVec(plusVec(minusVec(skalVec(F[isurfaces1[j].element1], S_r), skalVec(F[isurfaces1[j].element2], S_l)), \
					skalVec(minusVec(U[isurfaces1[j].element2], U[isurfaces1[j].element1]), S_l * S_r)), S_r - S_l);
			}
			else if (S_r <= 0)
			{
				F_hll = F[isurfaces1[j].element2];
			}
			F_res = multiVec(T, F_hll);
			U[isurfaces1[j].element1] = minusVec(U[isurfaces1[j].element1], skalVec(F_res, 0.1 * isurfaces1[j].area));
			U[isurfaces1[j].element2] = plusVec(U[isurfaces1[j].element1], skalVec(F_res, 0.1 * isurfaces1[j].area));
		}
		// Проходимся по граничным элементам
		cout << "Проходимся по граничным граням" << endl;
		for (int j = 0; j < count_boundary_surfaces1; ++j)
		{
			vector<double> F_hll(5);
			vector<vector<double> > T(5, vector<double>(5));
			vector<double> F_res(5);
			if (bsurfaces1[j].nz == 1)
			{
				T[0][0] = 1;
				T[1][3] = 1;
				T[2][2] = 1;
				T[4][4] = 1;
				T[3][1] = -1;
			}
			else if (bsurfaces1[j].nz == -1)
			{
				T[0][0] = 1;
				T[1][3] = -1;
				T[2][2] = -1;
				T[4][4] = 1;
				T[3][1] = 1;
			}
			else if (bsurfaces1[j].nz*bsurfaces1[j].nz != 1)
			{
				T[0][0] = 1;
				T[1][1] = bsurfaces1[j].nx;
				T[1][2] = bsurfaces1[j].ny;
				T[1][3] = bsurfaces1[j].nz;
				T[2][1] = -bsurfaces1[j].ny / (sqrt(1 - bsurfaces1[j].nz*bsurfaces1[j].nz));
				T[2][2] = bsurfaces1[j].nx / (sqrt(1 - bsurfaces1[j].nz*bsurfaces1[j].nz));
				T[3][1] = -bsurfaces1[j].nx*bsurfaces1[j].nz / (sqrt(1 - bsurfaces1[j].nz*bsurfaces1[j].nz));
				T[3][2] = -bsurfaces1[j].ny*bsurfaces1[j].nz / (sqrt(1 - bsurfaces1[j].nz*bsurfaces1[j].nz));
				T[3][3] = sqrt(1 - bsurfaces1[j].nz*bsurfaces1[j].nz);
				T[4][4] = 1;
			}
			T = transp(T, 5, 5);

			if (bsurfaces1[j].type == FreedomExit)
			{
				// Как вычислить в центре грани?
				cout << "сверхзвуковая граница выхода" << endl;
				double v_n = primitive[bsurfaces1[j].element1][1] * bsurfaces1[j].nx + primitive[bsurfaces1[j].element1][2] * bsurfaces1[j].ny + primitive[bsurfaces1[j].element1][3] * bsurfaces1[j].nz;
				F_hll[0] = primitive[bsurfaces1[j].element1][0] * v_n;
				F_hll[1] = primitive[bsurfaces1[j].element1][0] * primitive[bsurfaces1[j].element1][1] * v_n + primitive[bsurfaces1[j].element1][4] * bsurfaces1[j].nx;
				F_hll[2] = primitive[bsurfaces1[j].element1][0] * primitive[bsurfaces1[j].element1][2] * v_n + primitive[bsurfaces1[j].element1][4] * bsurfaces1[j].ny;
				F_hll[3] = primitive[bsurfaces1[j].element1][0] * primitive[bsurfaces1[j].element1][3] * v_n + primitive[bsurfaces1[j].element1][4] * bsurfaces1[j].nz;
				F_hll[4] = primitive[bsurfaces1[j].element1][0] * v_n * (U[bsurfaces1[j].element1][4] + primitive[bsurfaces1[j].element1][4]);
			}
			else if (bsurfaces1[j].type == HardBorder)
			{
				cout << "поверхность, представляющая собой твёрдую непроницаемую стенку" << endl;
				F_hll[0] = 0;
				F_hll[1] = primitive[bsurfaces1[j].element1][4] * bsurfaces1[j].nx;
				F_hll[2] = primitive[bsurfaces1[j].element1][4] * bsurfaces1[j].ny;
				F_hll[3] = primitive[bsurfaces1[j].element1][4] * bsurfaces1[j].nz;
				F_hll[4] = 0;
			}
			else if (bsurfaces1[j].type == Symmetry)
			{
				cout << "плоскость симметрии" << endl;
				F_hll[0] = 0;
				F_hll[1] = primitive[bsurfaces1[j].element1][4] * bsurfaces1[j].nx;
				F_hll[2] = primitive[bsurfaces1[j].element1][4] * bsurfaces1[j].ny;
				F_hll[3] = primitive[bsurfaces1[j].element1][4] * bsurfaces1[j].nz;
				F_hll[4] = 0;
			}
			else //(bsurfaces1[j].type == Custom)
			{
				cout << "Сверхзвуковая граница входа" << endl;
				double v_n_lim = v_lim[0] * bsurfaces1[j].nx + v_lim[1] * bsurfaces1[j].ny + v_lim[2] * bsurfaces1[j].nz;
				F_hll[0] = ro_lim * v_n_lim;
				F_hll[1] = ro_lim * v_lim[0] * v_n_lim + p_lim * bsurfaces1[j].nx;
				F_hll[2] = ro_lim * v_lim[1] * v_n_lim + p_lim * bsurfaces1[j].ny;
				F_hll[3] = ro_lim * v_lim[2] * v_n_lim + p_lim * bsurfaces1[j].nz;
				F_hll[4] = ro_lim * v_n_lim * (p_lim + 10 / 4 * p_lim + 0.5 * ro_lim * (v_lim[0] * v_lim[0] + v_lim[1] * v_lim[1] + v_lim[2] * v_lim[2]));
			}
			F_res = multiVec(T, F_hll);
			U[bsurfaces1[j].element1] = minusVec(U[bsurfaces1[j].element1], skalVec(F_res, 0.1 * bsurfaces1[j].area));
		}
		cout << "Переход к примитивным переменным" << endl;
		// Переход к примитивным элементам
		for (int j = 0; j < count_points1; ++j)
		{
			primitive[j][0] = U[j][0];
			primitive[j][1] = U[j][1] / primitive[j][0];
			primitive[j][2] = U[j][2] / primitive[j][0];
			primitive[j][3] = U[j][3] / primitive[j][0];
			primitive[j][4] = 0.4*(U[j][4] - primitive[j][0] * 0.5*(primitive[j][1] * primitive[j][1] + primitive[j][2] * primitive[j][2] + primitive[j][3] * primitive[j][3]));
		}
	}

	//Euler efvm;
	//efvm.SaveMeshInGMSHFile();
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


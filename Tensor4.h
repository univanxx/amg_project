#ifndef TENSOR4_H
#define TENSOR4_H

#include "Tensors.h"
#include <iostream>

class Tensor4
{
private:
	int axle1;
	int axle2;
	int rows;
	int columns;
	double**** elements;

public:
	Tensor4(int number_axis1, int number_axis2, int number_rows, int number_columns);  // Создаётся нулевой тензор заданных размеров
	Tensor4(const Tensor4& tensor);  // Копирование
	double**** getElements() const  // Получение элементов тензора
	{
		return elements;
	};
	int getAxle1() const  // Получение количества осей
	{
		return axle1;
	};
	int getAxle2() const  // Получение количества осей
	{
		return axle2;
	};
	int getRows() const  // Получение количества строк
	{
		return rows;
	};
	int getColumns() const  // Получение количества столбцов
	{
		return columns;
	};
	void setElement(int axle1, int axle2, int row, int column, double element);  // Передача элемента в заданную позицию
	void addElement(int axle1, int axle2, int row, int column, double element);  // Добавление элемента к заданному значению
	void setMatrix(int axle1, int axle2, const Matrix &M); // Постановка матрицы на позицию
	void addMatrix(int axle1, int axle2, const Matrix &M); // Добавление матрицы на позицию

	Tensor4& operator =(const Tensor4& tensor);  // Оператор присваивания
	Tensor3 operator [](const int index);

	friend std::ostream& operator<<(std::ostream&, const Tensor4& tensor);
	~Tensor4()  // Деструктор
	{
		for (int m = 0; m < axle1; ++m)
		{
			for (int i = 0; i < axle2; ++i)
			{
				for (int j = 0; j < rows; ++j)
				{
					delete[] elements[m][i][j];
				}
				delete elements[m][i];
			}
			delete[] elements[m];
		}
		delete[] elements;
	};
};


#endif

#ifndef TENSORS_H
#define TENSORS_H

#include "Matrices.h"
#include <iostream>

class Tensor3
{
private:
	int axis;
	int rows;
	int columns;
	double*** elements;

public:
	Tensor3(int number_axis, int number_rows, int number_columns);  // Создаётся нулевой тензор заданных размеров
	Tensor3(const Tensor3& tensor);  // Копирование
	double*** getElements() const  // Получение элементов тензора
	{
		return elements;
	};
	int getAxes() const  // Получение количества осей
	{
		return axis;
	};
	int getRows() const  // Получение количества строк
	{
		return rows;
	};
	int getColumns() const  // Получение количества столбцов
	{
		return columns;
	};
	void setElement(int axle, int row, int column, double element);  // Передача элемента в заданную позицию
	void addElement(int axle, int row, int column, double element);  // Добавление элемента к заданному значению

	Tensor3& operator =(const Tensor3& tensor);  // Оператор присваивания
	Matrix& operator [](const int index);

	friend std::ostream& operator<<(std::ostream&, const Tensor3& tensor);
	~Tensor3()  // Деструктор
	{
		for (int i = 0; i < axis; ++i)
		{
			for (int j = 0; j < rows; ++j)
			{
				delete[] elements[i][j];
			}
			delete elements[i];
		}
		delete[] elements;
	};
};


#endif
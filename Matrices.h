#ifndef MATRICES_H
#define MATRICES_H

#include "Vectors.h"
#include <iostream>

class Matrix
{
private:
	int rows;
	int columns;
	double** elements;

public:
	Matrix(int number_rows, int number_columns);  // Создаётся нулевая матрица заданных размеров
	Matrix(const Matrix& matrix);  // Копирование
	double** getElements() const  // Получение элементов матрицы
	{
		return elements;
	};
	int getRows() const  // Получение количества строк
	{
		return rows;
	};
	int getColumns() const  // Получение количества столбцов
	{
		return columns;
	};
	void setElement(int row, int column, double element);  // Передача элемента в заданную позицию
	void addElement(int row, int column, double element);  // Добавление элемента к заданному значению
	void substractElement(int row, int column, double element);  // Вычитание элемента из заданной позиции
	void addVector(int row, const Vector& vector);  // Добавление вектора к заданной строке
	Matrix transpose() const;  // Транспонирование матрицы
	Matrix& operator =(const Matrix& matrix);  // Оператор присваивания
	Matrix operator *(const Matrix& matrix);  // Оператор умножения на матрицу
	Matrix operator +(const Matrix& matrix);  // Оператор сложения матриц
	Vector operator *(const Vector& vector) const;  // Оператор умножения на вектор -> в результате получается вектор
	Vector operator [](const int index);
	Matrix operator *(const double& number);  // Оператор умножения на скаляр
	friend std::ostream& operator<<(std::ostream&, const Matrix&);
	~Matrix()  // Деструктор
	{
		for (int i = 0; i < rows; i++)
			delete[] elements[i];
		delete[] elements;
	};
};


#endif


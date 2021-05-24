#ifndef VECTORS_H
#define VECTORS_H

#include <math.h>
#include <iostream>

class Vector
{
private:
	int size;
	double* elements;

public:
	Vector(int dimension);  // Создаётся нулевой вектор заданной длины 
	Vector(const Vector& vector); // Копирование
	int getSize() const  // Размер вектора
	{
		return size;
	};
	double* getElements() const  // Получение элементов вектора
	{
		return elements;
	};
	void addElement(int position, double element);  // Добавляем скаляр к элементу с нужной позицией
	double findNorm() const;  // Норма L_{2}
	Vector& operator =(const Vector& vector);  // Оператор присваивания
	Vector& operator =(const double* vector);  // Ещё один ператор присваивания
	Vector operator *(const double& number) const;  // Оператор умножения на скаляр
	Vector operator /(const double& number) const;  // Оператор деления на скаляр
	double operator *(const Vector& vector) const;  // Оператор умножения на вектор
	Vector& operator +=(const  Vector& vector);  // Оператор сложения векторов
	Vector& operator -=(const  Vector& vector);  // Оператор вычитания векторов
	Vector operator +(const  Vector& vector) const;  // Ещё один оператор сложения векторов
	Vector operator -(const  Vector& vector) const;  // Ещё один оператор вычитанияя векторов
	double& operator [](int index) const  // Оператор взятия индекса
	{
		return elements[index];
	};
	friend std::ostream& operator<<(std::ostream&, const Vector&); // Функция вывода вектора
	~Vector()  // Деструктор
	{
		delete[] elements;
	};
};


#endif


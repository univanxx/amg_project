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
	Matrix(int number_rows, int number_columns);  // �������� ������� ������� �������� ��������
	Matrix(const Matrix& matrix);  // �����������
	double** getElements() const  // ��������� ��������� �������
	{
		return elements;
	};
	int getRows() const  // ��������� ���������� �����
	{
		return rows;
	};
	int getColumns() const  // ��������� ���������� ��������
	{
		return columns;
	};
	void setElement(int row, int column, double element);  // �������� �������� � �������� �������
	void addElement(int row, int column, double element);  // ���������� �������� � ��������� ��������
	void substractElement(int row, int column, double element);  // ��������� �������� �� �������� �������
	void addVector(int row, const Vector& vector);  // ���������� ������� � �������� ������
	Matrix transpose() const;  // ���������������� �������
	Matrix& operator =(const Matrix& matrix);  // �������� ������������
	Matrix operator *(const Matrix& matrix);  // �������� ��������� �� �������
	Matrix operator +(const Matrix& matrix);  // �������� �������� ������
	Vector operator *(const Vector& vector) const;  // �������� ��������� �� ������ -> � ���������� ���������� ������
	Vector operator [](const int index);
	Matrix operator *(const double& number);  // �������� ��������� �� ������
	friend std::ostream& operator<<(std::ostream&, const Matrix&);
	~Matrix()  // ����������
	{
		for (int i = 0; i < rows; i++)
			delete[] elements[i];
		delete[] elements;
	};
};


#endif


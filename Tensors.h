#ifndef TENSORS_H
#define TENSORS_H

#include "Matrices.h"
#include <iostream>

class Tensor
{
private:
	int axis;
	int rows;
	int columns;
	double*** elements;

public:
	Tensor(int number_axis, int number_rows, int number_columns);  // �������� ������� ������ �������� ��������
	Tensor(const Tensor& tensor);  // �����������
	double*** getElements() const  // ��������� ��������� �������
	{
		return elements;
	};
	int getAxes() const  // ��������� ���������� ����
	{
		return axis;
	};
	int getRows() const  // ��������� ���������� �����
	{
		return rows;
	};
	int getColumns() const  // ��������� ���������� ��������
	{
		return columns;
	};
	void setElement(int axle, int row, int column, double element);  // �������� �������� � �������� �������
	void addElement(int axle, int row, int column, double element);  // ���������� �������� � ��������� ��������

	Tensor& operator =(const Tensor& tensor);  // �������� ������������
	Matrix operator [](const int index);
	//Matrix operator *(const Matrix& matrix) const;  // �������� ��������� �� ������� -> � ���������� ���������� ������

	friend std::ostream& operator<<(std::ostream&, const Tensor& tensor);
	~Tensor()  // ����������
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
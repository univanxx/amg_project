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
	Tensor3(int number_axis, int number_rows, int number_columns);  // �������� ������� ������ �������� ��������
	Tensor3(const Tensor3& tensor);  // �����������
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

	Tensor3& operator =(const Tensor3& tensor);  // �������� ������������
	Matrix& operator [](const int index);

	friend std::ostream& operator<<(std::ostream&, const Tensor3& tensor);
	~Tensor3()  // ����������
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
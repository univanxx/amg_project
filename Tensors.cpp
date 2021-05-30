//#include "stdafx.h"
#include "Tensors.h"
#include "Matrices.h"
#include "Vectors.h"

Tensor3::Tensor3(int number_axis, int number_rows, int number_columns) : axis(number_axis), rows(number_rows), columns(number_columns)
{
	elements = new double** [axis];
	for (int i = 0; i < axis; ++i)
	{
		elements[i] = new double* [rows];
		for (int j = 0; j < rows; ++j)
		{
			elements[i][j] = new double [columns];
			for (int k = 0; k < columns; ++k)
			{
				elements[i][j][k] = 0.;
			}
		}
	}
}


Tensor3::Tensor3(const Tensor3& tensor) : axis(tensor.axis), rows(tensor.rows), columns(tensor.columns)
{
	elements = new double** [axis];
	for (int i = 0; i < axis; ++i)
	{
		elements[i] = new double* [rows];
		for (int j = 0; j < rows; ++j)
		{
			elements[i][j] = new double[columns];
			for (int k = 0; k < columns; ++k)
			{
				elements[i][j][k] = tensor.getElements()[i][j][k];
			}
		}
	}
}


void Tensor3::setElement(int axle, int row, int column, double element)
{
	(*this).elements[axle][row][column] = element;
}

void Tensor3::addElement(int axle, int row, int column, double element)
{
	(*this).elements[axle][row][column] += element;
}


Tensor3& Tensor3::operator =(const Tensor3& tensor)
{
	axis = tensor.getAxes();
	rows = tensor.getRows();
	columns = tensor.getColumns();
	for (int i = 0; i < axis; ++i)
	{
		for (int j = 0; j < rows; ++j)
		{
			for (int k = 0; k < columns; ++k)
			{
				elements[i][j][k] = tensor.getElements()[i][j][k];
			}
		}
	}
	return (*this);
}

Matrix& Tensor3::operator [](const int index)
{
	int len_column = (*this).getColumns();
	int len_row = (*this).getRows();
	Matrix res(len_row, len_column);
	for (int k = 0; k < len_row; ++k)
	{
		for (int j = 0; j < len_column; ++j)
		{
			res.setElement(k, j, (*this).getElements()[index][k][j]);
		}
	}
	return res;
}



std::ostream& operator<<(std::ostream& out_stream, const Tensor3& tensor)
{
	for (int i = 0; i < tensor.axis; ++i)
	{
		out_stream << "Axle number " << i << "\n";
		for (int j = 0; j < tensor.rows; ++j)
		{
			for (int k = 0; k < tensor.columns; ++k)
			{
				out_stream << tensor.getElements()[i][j][k] << "  ";
			}
			out_stream << "\n";
		}
		out_stream << "\n";
	}
	return out_stream;
}
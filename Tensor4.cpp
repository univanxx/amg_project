//#include "stdafx.h"
#include "Tensor4.h"
#include "Tensors.h"
#include "Matrices.h"
#include "Vectors.h"

Tensor4::Tensor4(int number_axis1, int number_axis2, int number_rows, int number_columns) : axle1(number_axis1), axle2(number_axis2), rows(number_rows), columns(number_columns)
{
	elements = new double***[axle1];
	for (int i = 0; i < axle1; ++i)
	{
		elements[i] = new double** [axle2];
		for (int j = 0; j < axle2; ++j)
		{
			elements[i][j] = new double* [rows];
			for (int k = 0; k < rows; ++k)
			{
				elements[i][j][k] = new double[columns];
				for (int l = 0; l < columns; ++l)
				{
					elements[i][j][k][l] = 0.;
				}
			}
		}
	}
}


Tensor4::Tensor4(const Tensor4& tensor) : axle1(tensor.axle1), axle2(tensor.axle2), rows(tensor.rows), columns(tensor.columns)
{
	elements = new double*** [axle1];
	for (int i = 0; i < axle1; ++i)
	{
		elements[i] = new double** [axle2];
		for (int j = 0; j < axle2; ++j)
		{
			elements[i][j] = new double* [rows];
			for (int k = 0; k < rows; ++k)
			{
				elements[i][j][k] = new double[columns];
				for (int l = 0; l < columns; ++l)
				{
					elements[i][j][k][l] = tensor.getElements()[i][j][k][l];
				}
			}
		}
	}
}


void Tensor4::setElement(int axle1, int axle2, int row, int column, double element)
{
	(*this).elements[axle1][axle2][row][column] = element;
}

void Tensor4::addElement(int axle1, int axle2, int row, int column, double element)
{
	(*this).elements[axle1][axle2][row][column] += element;
}

void Tensor4::setMatrix(int axle1, int axle2, const Matrix &M)
{
	for (int i = 0; i < M.getRows(); ++i)
	{
		for (int j = 0; j < M.getColumns(); ++j)
		{
			(*this).setElement(axle1, axle2, i, j, M.getElements()[i][j]);
		}
	}
}

void Tensor4::addMatrix(int axle1, int axle2, const Matrix &M)
{
	for (int i = 0; i < M.getRows(); ++i)
	{
		for (int j = 0; j < M.getColumns(); ++j)
		{
			(*this).addElement(axle1, axle2, i, j, M.getElements()[i][j]);
		}
	}
}


Tensor4& Tensor4::operator =(const Tensor4& tensor)
{
	axle1 = tensor.getAxle1();
	axle2 = tensor.getAxle2();
	rows = tensor.getRows();
	columns = tensor.getColumns();
	{
		elements = new double*** [axle1];
		for (int i = 0; i < axle1; ++i)
		{
			elements[i] = new double** [axle2];
			for (int j = 0; j < axle2; ++j)
			{
				elements[i][j] = new double* [rows];
				for (int k = 0; k < rows; ++k)
				{
					elements[i][j][k] = new double[columns];
					for (int l = 0; l < columns; ++l)
					{
						elements[i][j][k][l] = tensor.getElements()[i][j][k][l];
					}
				}
			}
		}
	}
	return (*this);
}

Tensor3 Tensor4::operator [](const int index)
{
	int axle = (*this).getAxle2();
	int len_column = (*this).getColumns();
	int len_row = (*this).getRows();
	Tensor3 res(axle, len_row, len_column);
	for (int i = 0; i < axle; ++i)
	{
		for (int k = 0; k < len_row; ++k)
		{
			for (int j = 0; j < len_column; ++j)
			{
				res.setElement(i, k, j, (*this).getElements()[index][i][k][j]);
			}
		}
	}
	return res;
}



std::ostream& operator<<(std::ostream& out_stream, const Tensor4& tensor)
{
	for (int m = 0; m < tensor.axle1; ++m)
	{
		out_stream << "Axle1 number " << m << "\n";
		for (int i = 0; i < tensor.axle2; ++i)
		{
			out_stream << "Axle2 number " << i << "\n";
			for (int j = 0; j < tensor.rows; ++j)
			{
				for (int k = 0; k < tensor.columns; ++k)
				{
					out_stream << tensor.getElements()[m][i][j][k] << "  ";
				}
				out_stream << "\n";
			}
			out_stream << "\n";
		}
	}
	return out_stream;
}
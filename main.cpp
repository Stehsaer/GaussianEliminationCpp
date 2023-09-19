// main.cpp
// ���ߣ�SYSU �����ѧԺ ���Ž�
// ���ڣ�2023-9-19

#include <iostream>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <vector>
#include <tuple>

// ����
struct vector
{
private:
	size_t dimension;
	double* data;

public:

	// ���������
	vector()
	{
		dimension = 0;
		data = nullptr;
	}

	// ����
	vector(size_t dim)
	{
		dimension = dim;
		data = new double[dim];

		// ��ʼ��Ϊ0
		for (size_t i = 0; i < dim; i++)
		{
			data[i] = 0.0;
		}
	}

	// ����
	vector(std::initializer_list<double>&& val)
	{
		dimension = val.size();
		data = new double[dimension];

		for (size_t i = 0; i < dimension; i++)
			data[i] = *(val.begin() + i);
	}

	// ���ƹ���
	vector(const vector& val)
	{
		dimension = val.dimension;
		data = new double[dimension];
		memcpy(data, val.data, sizeof(double) * dimension); // ���ƹ��죺��������
	}

	// �ƶ�����
	vector(vector&& val) noexcept
	{
		dimension = val.dimension;
		data = val.data;
		val.data = nullptr;
	}

	// ��������
	~vector()
	{
		if(data)
			delete data;
	}

	//== ��Ԫ���������� ==

	vector operator+(const vector& val)
	{
		if (dimension != val.dimension) throw std::runtime_error("ά�Ȳ�ͬ���޷�����");

		vector r(dimension);

		for (size_t i = 0; i < dimension; i++)
		{
			r.data[i] = val.data[i] + data[i];
		}

		return r;
	}

	// ԭλ�ӷ�������
	void operator+=(const vector& val)
	{
		if (dimension != val.dimension) throw std::runtime_error("ά�Ȳ�ͬ���޷�����");

		for (size_t i = 0; i < dimension; i++)
		{
			data[i] += val.data[i];
		}
	}

	vector operator-(const vector& val)
	{
		if (dimension != val.dimension) throw std::runtime_error("ά�Ȳ�ͬ���޷�����");

		vector r(dimension);

		for (size_t i = 0; i < dimension; i++)
		{
			r.data[i] = val.data[i] + data[i];
		}

		return r;
	}

	// ԭλ����������
	void operator-=(const vector& val)
	{
		if (dimension != val.dimension) throw std::runtime_error("ά�Ȳ�ͬ���޷�����");

		for (size_t i = 0; i < dimension; i++)
		{
			data[i] -= val.data[i];
		}
	}

	vector operator*(double val)
	{
		vector r(dimension);

		for (size_t i = 0; i < dimension; i++)
		{
			r.data[i] = val * data[i];
		}

		return r;
	}

	// ԭλ�˷�������
	void operator*=(double val)
	{
		for (size_t i = 0; i < dimension; i++)
		{
			data[i] *= val;
		}
	}

	vector operator/(double val)
	{
		vector r(dimension);

		for (size_t i = 0; i < dimension; i++)
		{
			r.data[i] = val * data[i];
		}

		return r;
	}

	// ԭλ����������
	void operator/=(double val)
	{
		for (size_t i = 0; i < dimension; i++)
		{
			data[i] /= val;
		}
	}

	// ��ֵ
	void operator=(const vector& val)
	{
		if (data) delete data;
		
		dimension = val.dimension;
		data = new double[dimension];
		memcpy(data, val.data, dimension * sizeof(double));
	}

	// ��ȡĳ��
	double& operator[](size_t index) const
	{
		if (index >= dimension) throw std::runtime_error("��������");

		return data[index];
	}

	// ת�ַ���
	std::string str() const
	{
		std::stringstream stream;
		stream << "{";
		
		for (size_t i = 0; i < dimension; i++)
		{
			stream << data[i];

			if(i < dimension - 1)
				stream << ",";
		}

		stream << "}";

		return stream.str();
	}

	// ��ȡ����
	size_t dim() const
	{
		return dimension;
	}

	// ����������
	static void swap(vector& left, vector& right)
	{
		std::swap(left.dimension, right.dimension);
		std::swap(left.data, right.data);
	}
};

// ���д洢�ľ���
struct matrix
{
private:
	size_t width, height;
	vector* vectors;

public:

	// ֱ��ͨ����߹������
	matrix(size_t w, size_t h)
	{
		width = w; height = h;
		vectors = new vector[h];

		for (size_t i = 0; i < h; i++)
		{
			vectors[i] = vector(w);
		}
	}

	// ͨ����ʼ���б������
	matrix(std::initializer_list<vector>&& val)
	{
		width = val.begin()->dim();
		height = val.size();

		vectors = new vector[height];

		for (size_t i = 0; i < height; i++)
		{
			auto& v = *(val.begin() + i);

			if (v.dim() != width) throw std::runtime_error("Consistent dimension is required");

			vectors[i] = v;
		}
	}

	// ��������
	~matrix()
	{
		if (vectors)
			delete[] vectors;
	}

	// ȡĳ��
	vector& operator[](size_t index)
	{
		if (index > height) throw std::runtime_error("��������");

		return vectors[index];
	}

	// ��ȡ�߶�
	inline size_t get_width() const
	{
		return width;
	}

	// ��ȡ���
	inline size_t get_height() const
	{
		return height;
	}

	// ������
	void swap_row(size_t row1, size_t row2)
	{
		vector::swap((*this)[row1], (*this)[row2]);
	}

	// Ѱ����Ԫ֮�µ�������ֵ
	size_t get_largest_abs_row(size_t x, size_t pivot)
	{
		size_t largest_row = 0;
		double largest_value = 0;

		for (size_t i = pivot; i < height; i++)
		{
			if (abs(vectors[i][x]) > largest_value)
			{
				largest_value = abs(vectors[i][x]);
				largest_row = i;
			}
		}

		// ���ֵΪ0������SIZE_MAX
		if (largest_value == 0)
			return SIZE_MAX;

		return largest_row;
	}

	// print������̨
	void print()
	{
		std::cout << "{" << std::endl;
		for (size_t i = 0; i < height; i++)
		{
			std::cout << "    " << vectors[i].str() << std::endl;
		}
		std::cout << "}" << std::endl;
	}
};

int main()
{
	using pivot_pos = std::tuple<size_t, size_t>; // ��¼��Ԫλ��Ԫ�飬��ʽΪ[�У���]

	matrix m = {
		{0, 3, -6, 6, 4, -5},
		{3, -7, 8, -5, 8, 9},
		{3, -9, 12, -9, 6, 15}
	};

	std::cout << "����ǰ" << std::endl;
	m.print();

	// ǰ����Ԫ
	// ����Ԫ�н�����Ԫ

	size_t pivot_row = 0;

	// ��¼��Ԫ�У����ں�������
	std::vector<pivot_pos> pivot_positions;

	for (size_t column = 0; column < m.get_width(); column++)
	{
		// ��ȡ��Ԫ����������ֵ��
		size_t largest_row = m.get_largest_abs_row(column, pivot_row);

		// ������Ԫ���¾�Ϊ0
		if (largest_row == SIZE_MAX) continue;

		// ��¼��Ԫλ��
		pivot_positions.push_back(pivot_pos(column, pivot_row));

		// ������л�����Ԫ��
		if(pivot_row != largest_row)
			m.swap_row(pivot_row, largest_row);

		// ����Ԫ֮��ȫ����Ϊ0
		for (size_t row = pivot_row + 1; row < m.get_height(); row++)
		{
			m[row] -= m[pivot_row] * (m[row][column] / m[pivot_row][column]);
		}

		// ��ǰ����Ԫ��Ϊ1
		m[pivot_row] /= m[pivot_row][column]; 

		pivot_row++;

		// ��Ԫ�任���
		if (pivot_row == m.get_height()) break;		
	}

	// ������Ԫ
	for (auto pos : pivot_positions)
	{
		// �ṹ���󶨣�ע�⣺��Ҫ����ΪC++17����߰汾��
		auto& [p_column, p_row] = pos;

		// ����Ԫ֮��ȫ����Ϊ0
		for (size_t i = 0; i < p_row; i++)
		{
			m[i] -= m[p_row] * (m[i][p_column] / m[p_row][p_column]);
		}
	}

	std::cout << "�����" << std::endl;
	m.print();

	return 0;
}
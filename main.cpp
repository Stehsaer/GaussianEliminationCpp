// main.cpp
// 作者：SYSU 计算机学院 刘信杰
// 日期：2023-9-19

#include <iostream>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <vector>
#include <tuple>

// 向量
struct vector
{
private:
	size_t dimension;
	double* data;

public:

	// 构造空向量
	vector()
	{
		dimension = 0;
		data = nullptr;
	}

	// 构造
	vector(size_t dim)
	{
		dimension = dim;
		data = new double[dim];

		// 初始化为0
		for (size_t i = 0; i < dim; i++)
		{
			data[i] = 0.0;
		}
	}

	// 构造
	vector(std::initializer_list<double>&& val)
	{
		dimension = val.size();
		data = new double[dimension];

		for (size_t i = 0; i < dimension; i++)
			data[i] = *(val.begin() + i);
	}

	// 复制构造
	vector(const vector& val)
	{
		dimension = val.dimension;
		data = new double[dimension];
		memcpy(data, val.data, sizeof(double) * dimension); // 复制构造：拷贝数据
	}

	// 移动构造
	vector(vector&& val) noexcept
	{
		dimension = val.dimension;
		data = val.data;
		val.data = nullptr;
	}

	// 析构向量
	~vector()
	{
		if(data)
			delete data;
	}

	//== 逐元素四则运算 ==

	vector operator+(const vector& val)
	{
		if (dimension != val.dimension) throw std::runtime_error("维度不同，无法操作");

		vector r(dimension);

		for (size_t i = 0; i < dimension; i++)
		{
			r.data[i] = val.data[i] + data[i];
		}

		return r;
	}

	// 原位加法，更快
	void operator+=(const vector& val)
	{
		if (dimension != val.dimension) throw std::runtime_error("维度不同，无法操作");

		for (size_t i = 0; i < dimension; i++)
		{
			data[i] += val.data[i];
		}
	}

	vector operator-(const vector& val)
	{
		if (dimension != val.dimension) throw std::runtime_error("维度不同，无法操作");

		vector r(dimension);

		for (size_t i = 0; i < dimension; i++)
		{
			r.data[i] = val.data[i] + data[i];
		}

		return r;
	}

	// 原位减法，更快
	void operator-=(const vector& val)
	{
		if (dimension != val.dimension) throw std::runtime_error("维度不同，无法操作");

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

	// 原位乘法，更快
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

	// 原位除法，更快
	void operator/=(double val)
	{
		for (size_t i = 0; i < dimension; i++)
		{
			data[i] /= val;
		}
	}

	// 赋值
	void operator=(const vector& val)
	{
		if (data) delete data;
		
		dimension = val.dimension;
		data = new double[dimension];
		memcpy(data, val.data, dimension * sizeof(double));
	}

	// 获取某项
	double& operator[](size_t index) const
	{
		if (index >= dimension) throw std::runtime_error("索引超限");

		return data[index];
	}

	// 转字符串
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

	// 获取项数
	size_t dim() const
	{
		return dimension;
	}

	// 交换两向量
	static void swap(vector& left, vector& right)
	{
		std::swap(left.dimension, right.dimension);
		std::swap(left.data, right.data);
	}
};

// 按行存储的矩阵
struct matrix
{
private:
	size_t width, height;
	vector* vectors;

public:

	// 直接通过宽高构造矩阵
	matrix(size_t w, size_t h)
	{
		width = w; height = h;
		vectors = new vector[h];

		for (size_t i = 0; i < h; i++)
		{
			vectors[i] = vector(w);
		}
	}

	// 通过初始化列表构造矩阵
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

	// 析构矩阵
	~matrix()
	{
		if (vectors)
			delete[] vectors;
	}

	// 取某行
	vector& operator[](size_t index)
	{
		if (index > height) throw std::runtime_error("索引超限");

		return vectors[index];
	}

	// 获取高度
	inline size_t get_width() const
	{
		return width;
	}

	// 获取宽度
	inline size_t get_height() const
	{
		return height;
	}

	// 交换行
	void swap_row(size_t row1, size_t row2)
	{
		vector::swap((*this)[row1], (*this)[row2]);
	}

	// 寻找主元之下的最大绝对值
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

		// 最大值为0，返回SIZE_MAX
		if (largest_value == 0)
			return SIZE_MAX;

		return largest_row;
	}

	// print到控制台
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
	using pivot_pos = std::tuple<size_t, size_t>; // 记录主元位置元组，格式为[行，列]

	matrix m = {
		{0, 3, -6, 6, 4, -5},
		{3, -7, 8, -5, 8, 9},
		{3, -9, 12, -9, 6, 15}
	};

	std::cout << "化简前" << std::endl;
	m.print();

	// 前向消元
	// 逐主元列进行消元

	size_t pivot_row = 0;

	// 记录主元列，便于后续化简
	std::vector<pivot_pos> pivot_positions;

	for (size_t column = 0; column < m.get_width(); column++)
	{
		// 获取主元行下最大绝对值项
		size_t largest_row = m.get_largest_abs_row(column, pivot_row);

		// 该列主元行下均为0
		if (largest_row == SIZE_MAX) continue;

		// 记录主元位置
		pivot_positions.push_back(pivot_pos(column, pivot_row));

		// 将最大行换到主元处
		if(pivot_row != largest_row)
			m.swap_row(pivot_row, largest_row);

		// 将主元之下全部化为0
		for (size_t row = pivot_row + 1; row < m.get_height(); row++)
		{
			m[row] -= m[pivot_row] * (m[row][column] / m[pivot_row][column]);
		}

		// 提前把主元化为1
		m[pivot_row] /= m[pivot_row][column]; 

		pivot_row++;

		// 主元变换完成
		if (pivot_row == m.get_height()) break;		
	}

	// 后向消元
	for (auto pos : pivot_positions)
	{
		// 结构化绑定；注意：需要设置为C++17或更高版本！
		auto& [p_column, p_row] = pos;

		// 将主元之上全部化为0
		for (size_t i = 0; i < p_row; i++)
		{
			m[i] -= m[p_row] * (m[i][p_column] / m[p_row][p_column]);
		}
	}

	std::cout << "化简后" << std::endl;
	m.print();

	return 0;
}
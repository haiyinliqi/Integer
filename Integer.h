#ifndef _INTEGER_H_
#define _INTEGER_H_
#include<vector>
#include<string>
#include<cmath>
namespace NUMBER {
	using num_t = unsigned int;
	inline const num_t SIZE = 8;
	inline const num_t MOD = std::pow(10, SIZE);
	class Integer {//整数类
	protected:
		std::vector<num_t> _number;//整数
		bool _sign;//符号
		//绝对值判小于
		bool absolute_less(const Integer& other) const noexcept;
		//绝对值判大于
		bool absolute_​​greater(const Integer& other) const noexcept;
		//无符号加法
		Integer unsigned_add(const Integer& other) const noexcept;
		//无符号减法（this必须大于other）
		Integer unsigned_sub(const Integer& other) const noexcept;
		//普通无符号乘法运算
		Integer base_mult(const Integer& other) const noexcept;
		//Karatsuba算法（无符号）
		Integer Karatsuba(const Integer& other) const noexcept;
		//FFT算法（无符号）
		Integer FFT_mult(const Integer& other) const noexcept;
	public:
		Integer() noexcept;
		//通过char*构造
		Integer(const char* number) noexcept;
		//通过std::string构造
		Integer(const std::string& number) noexcept;
		//通过int构造
		Integer(const int& number) noexcept;
		//通过std::vector<num_t>构造
		Integer(const std::vector<num_t> number, const bool sign = false) noexcept;
		//转化为std::string
		std::string str() const noexcept;
		//运算重载
		bool operator==(const Integer& other) const noexcept;
		bool operator!=(const Integer& other) const noexcept;
		bool operator<(const Integer& other) const noexcept;
		bool operator>(const Integer& other) const noexcept;
		bool operator<=(const Integer& other) const noexcept;
		bool operator>=(const Integer& other) const noexcept;
		Integer operator-() const noexcept;
		Integer operator+(const Integer& other) const noexcept;
		void operator+=(const Integer& other) noexcept;
		Integer operator-(const Integer& other) const noexcept;
		void operator-=(const Integer& other) noexcept;
		Integer operator*(const Integer& other) const noexcept;
		void operator*=(const Integer& other) noexcept;
		Integer operator%(const Integer& other) const;
		void operator%=(const Integer& other);
		Integer operator/(const Integer& other) const;
		void operator/=(const Integer& other);
		//乘方运算
		Integer pow(const Integer& other) const noexcept;
		//带模数的乘方运算
		Integer pow(const Integer& other, const Integer& mod) const;
	};
	inline const Integer _0 = 0, _1 = 1, _2 = 2;
}
#endif


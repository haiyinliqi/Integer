#include"Integer.h"
#include<sstream>
#include<deque>
#include<cmath>
#include<stdexcept>
#include<iostream>
#include<functional>
#include<complex>
#include<numbers>
namespace NUMBER {
	bool Integer::absolute_less(const Integer& other) const noexcept {
		if (other._number.size() != _number.size()) {

			return _number.size() < other._number.size();
		}
		for (int i = _number.size() - 1;i >= 0;i--) {
			if (_number[i] != other._number[i]) {

				return _number[i] < other._number[i];
			}
		}

		return false;
	}
	bool Integer::absolute_greater(const Integer& other) const noexcept {
		if (other._number.size() != _number.size()) {

			return _number.size() > other._number.size();
		}
		for (int i = _number.size() - 1;i >= 0;i--) {
			if (_number[i] != other._number[i]) {

				return _number[i] > other._number[i];
			}
		}

		return false;
	}
	Integer Integer::unsigned_add(const Integer& other) const noexcept {
		Integer ans;
		ans._number.resize(std::max(_number.size(), other._number.size()));
		int size = std::min(_number.size(), other._number.size());
		unsigned long long carry = 0;
		for (int i = 0;i < size;i++) {
			carry += static_cast<unsigned long long>(_number[i]) + other._number[i];
			ans._number[i] = carry % MOD;
			carry /= MOD;
		}
		for (int i = size;i < _number.size();i++) {
			carry += _number[i];
			ans._number[i] = carry % MOD;
			carry /= MOD;
		}
		for (int i = size;i < other._number.size();i++) {
			carry += other._number[i];
			ans._number[i] = carry % MOD;
			carry /= MOD;
		}
		if (carry) {
			ans._number.push_back(carry);
		}

		return ans;
	}
	Integer Integer::unsigned_sub(const Integer& other) const noexcept {
		Integer ans = *this;
		unsigned long long borrow = 0;
		for (int i = 0;i < other._number.size();i++) {
			borrow += other._number[i];
			if (ans._number[i] < borrow) {
				ans._number[i] += MOD - borrow;
				borrow = 1;
			}
			else {
				ans._number[i] -= borrow;
				borrow = 0;
			}
		}
		for (int i = other._number.size();borrow;i++) {
			if (ans._number[i] < borrow) {
				ans._number[i] += MOD - borrow;
			}
			else {
				ans._number[i] -= borrow;
				borrow = 0;
			}
		}
		while (ans._number.size() > 1 && *(ans._number.end() - 1) == 0) {
			ans._number.pop_back();
		}

		return ans;
	}
	Integer Integer::base_mult(const Integer& other) const noexcept {
		Integer ans;
		ans._number.resize(_number.size() + other._number.size() - 1, 0);
		for (std::size_t i = 0;i < _number.size();i++) {
			unsigned long long carry = 0;
			for (std::size_t j = 0;j < other._number.size();j++) {
				carry += static_cast<unsigned long long>(_number[i]) * other._number[j] + ans._number[i + j];
				ans._number[i + j] = carry % MOD;
				carry /= MOD;
			}
			if (carry) {
				if (i + 1 < _number.size()) {
					ans._number[i + other._number.size()] = carry;
				}
				else {
					ans._number.push_back(carry);
				}
			}
		}

		return ans;
	}
	Integer Integer::Karatsuba(const Integer& other) const noexcept {
		std::size_t size = std::max(_number.size(), other._number.size());
		size /= 2;
		Integer a(std::vector<num_t>(size > _number.size() ? _number.end() : _number.begin() + size, _number.end()), false);
		Integer b(std::vector<num_t>(_number.begin(), size > _number.size() ? _number.end() : _number.begin() + size), false);
		Integer c(std::vector<num_t>(size > other._number.size() ? other._number.end() : other._number.begin() + size, other._number.end()), false);
		Integer d(std::vector<num_t>(other._number.begin(), size > other._number.size() ? other._number.end() : other._number.begin() + size), false);
		bool flag = false;
		if (a._number.size() == 0) {
			a._number.push_back(0);
			flag = true;
		}
		if (c._number.size() == 0) {
			c._number.push_back(0);
			flag = true;
		}
		Integer p1 = a * c, p2 = b * d, p3 = (a + b) * (c + d);
		p3 -= p1 + p2;
		Integer p4, ans;
		p4._number.resize(size);
		p4._number.insert(p4._number.end(), p3._number.begin(), p3._number.end());
		if (flag) {
			ans = _0;
		}
		else {
			ans._number.resize(size * 2);
			ans._number.insert(ans._number.end(), p1._number.begin(), p1._number.end());
		}
		ans += p4 + p2;

		return ans;
	}
	Integer Integer::FFT_mult(const Integer& other) const noexcept {
		using plural = std::complex<double>;
		std::function<std::vector<plural>(std::vector<plural>&, std::vector<plural>&)> FFT = [](std::vector<plural>& func, std::vector<plural>& val) -> std::vector<plural> {
			std::size_t n = func.size();
			std::vector<plural> ordered_func = func;
			std::vector<int> rev(n, 0);
			std::size_t log_n = std::log2(n);
			for (std::size_t i = 0; i < n; ++i) {
				rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (log_n - 1));
			}
			for (std::size_t i = 0; i < n; ++i) {
				if (i < rev[i]) {
					std::swap(ordered_func[i], ordered_func[rev[i]]);
				}
			}
			for (std::size_t len = 2; len <= n; len <<= 1) {
				std::size_t half_len = len / 2;
				std::size_t step = n / len;
				for (std::size_t start = 0; start < n; start += len) {
					for (std::size_t j = 0; j < half_len; ++j) {
						plural w = val[j * step];
						plural u = ordered_func[start + j];
						plural v = w * ordered_func[start + j + half_len];
						ordered_func[start + j] = u + v;
						ordered_func[start + j + half_len] = u - v;
					}
				}
			}

			return ordered_func;
			};
		std::size_t(*upper_power)(std::size_t) = [](std::size_t num)->std::size_t {
			std::size_t ans = 1;
			while (ans < num) {
				ans <<= 1;
			}

			return ans;
			};
		int k = (_number.size() + other._number.size()) * 2;
		k = upper_power(k);
		std::vector<plural> a(k), b(k);
		for (std::size_t i = 0;i < _number.size();i++) {
			a[i * 2] = _number[i] % (num_t)1e4;
			a[i * 2 + 1] = _number[i] / (num_t)1e4;
		}
		for (std::size_t i = 0;i < other._number.size();i++) {
			b[i * 2] = other._number[i] % (num_t)1e4;
			b[i * 2 + 1] = other._number[i] / (num_t)1e4;
		}
		std::vector<plural> val(k);
		for (int i = 0;i < k;i++) {
			val[i].real(std::cos(3.1415926535897932384624l * 2 * i / k));
			val[i].imag(std::sin(3.1415926535897932384624l * 2 * i / k));
		}
		std::vector<plural> ans1 = FFT(a, val);
		std::vector<plural> ans2 = FFT(b, val);
		for (std::size_t i = 0;i < k;i++) {
			ans1[i] *= ans2[i];
			ans1[i] = std::conj(ans1[i]);
		}
		std::vector<plural> product(k);
		product = FFT(ans1, val);
		for (std::size_t i = 0;i < k;i++) {
			product[i] = std::conj(product[i]);
			product[i] /= k;
		}
		Integer ans;
		unsigned long long carry = 0;
		for (std::size_t i = 0;i < k;i++) {
			if (i % 2 == 1) {
				carry += ((unsigned long long)std::llround(product[i].real())) * 1e4;
				ans._number.push_back(carry % MOD);
				carry /= MOD;
			}
			else {
				carry += std::llround(product[i].real());
			}
		}
		while (carry) {
			ans._number.push_back(carry % MOD);
			carry /= MOD;
		}
		while (ans._number.back() == 0) {
			ans._number.pop_back();
		}

		return ans;
	}
	Integer::Integer() noexcept :
		_number(), _sign(false) {
	}
	Integer::Integer(const char* number) noexcept :
		Integer(std::string(number)){
	}
	Integer::Integer(const std::string& number) noexcept {
		std::string unsigned_number;
		if (*number.begin() == '-') {
			_sign = true;
			unsigned_number = number.substr(1);
		}
		else {
			_sign = false;
			unsigned_number = number;
		}
		_number.resize((unsigned_number.size() + SIZE - 1) / SIZE);
		for (int i = unsigned_number.size() - 1;i >= 0;i -= SIZE) {
			if (i < SIZE - 1) {
				std::stringstream sstream(unsigned_number.substr(0, i + 1));
				sstream >> _number[(unsigned_number.size() - 1 - i) / SIZE];
			}
			else {
				std::stringstream sstream(unsigned_number.substr(i - SIZE + 1, SIZE));
				sstream >> _number[(unsigned_number.size() - 1 - i) / SIZE];
			}
		}
	}
	Integer::Integer(const int& number) noexcept :
		_number({ static_cast<unsigned int>(std::abs(number)) }), _sign(number < 0) {
	}
	Integer::Integer(const std::vector<num_t> number, const bool sign) noexcept :
		_number(number), _sign(sign) {
	}
	std::string Integer::str() const noexcept {
		std::deque<std::string> stack;
		for (std::size_t i = 0;i < _number.size();i++) {
			std::string new_str = std::to_string(_number[i]);
			while (new_str.size() < SIZE && i + 1 < _number.size()) {
				new_str = "0" + new_str;
			}
			stack.push_front(new_str);
		}
		if (_sign) {
			stack.push_front("-");
		}
		std::string ans;
		while (!stack.empty()) {
			ans += stack.front();
			stack.pop_front();
		}

		return ans;
	}
	bool Integer::operator==(const Integer& other) const noexcept {

		return _number == other._number && _sign == other._sign;
	}
	bool Integer::operator!=(const Integer& other) const noexcept {

		return !(_number == other._number);
	}
	bool Integer::operator<(const Integer& other) const noexcept {
		if (_sign ^ other._sign) {

			return _sign;
		}

		return _sign ? absolute_greater(other) : absolute_less(other);
	}
	bool Integer::operator>(const Integer& other) const noexcept {
		if (_sign ^ other._sign) {

			return !_sign;
		}

		return _sign ? absolute_less(other) : absolute_greater(other);
	}

	bool Integer::operator<=(const Integer& other) const noexcept {

		return !(*this > other);
	}
	bool Integer::operator>=(const Integer& other) const noexcept {

		return !(*this < other);
	}
	Integer Integer::operator-() const noexcept {
		if (_number.size() == 0 || (_number.size() == 1 && _number[0] == 0)) {

			return *this;
		}
		Integer ans = *this;
		ans._sign = !_sign;

		return ans;
	}
	Integer Integer::operator+(const Integer& other) const noexcept {
		Integer ans;
		if (_sign ^ other._sign) {
			bool sign = absolute_greater(other);
			ans = sign ? unsigned_sub(other) : other.unsigned_sub(*this);
			ans._sign = sign ? _sign : other._sign;

			return ans;
		}
		ans = this->unsigned_add(other);
		ans._sign = _sign;

		return ans;
	}
	void Integer::operator+=(const Integer& other) noexcept {
		*this = *this + other;
	}
	Integer Integer::operator-(const Integer& other) const noexcept {

		return *this + -other;
	}
	void Integer::operator-=(const Integer& other) noexcept {
		*this = *this - other;
	}
	Integer Integer::operator*(const Integer& other) const noexcept {
		if (*this == _0 || other == _0) {

			return _0;
		}
		Integer ans;
		if (_number.size() + other._number.size() > 2000 && _number.size() > 100 && other._number.size() > 100) {
			ans = FFT_mult(other);
		}
		else if (_number.size() + other._number.size() > 1000 && _number.size() > 50 && other._number.size() > 50) {
			ans = Karatsuba(other);
		}
		else {
			ans = base_mult(other);
		}
		ans._sign = _sign ^ other._sign;

		return ans;
	}
	void Integer::operator*=(const Integer& other) noexcept {
		*this = *this * other;
	}
	Integer Integer::operator%(const Integer& other) const {
		if (other == _0) {

			throw std::domain_error("除数不能为0");
		}
		else if (other == _2) {

			return std::vector<num_t>({ _number[0] % 2 });
		}
		else if (other == 5) {

			return std::vector<num_t>({ _number[0] % 5 });
		}
		Integer ans = _sign ? -(*this) : *this,
			abs_other = other._sign ? -other : other;
		std::vector<Integer> product;
		product.push_back(abs_other);
		product.push_back(abs_other + abs_other);
		while (product.back() <= ans) {
			product.push_back(product.back() + product.back());
		}
		std::size_t ptr = product.size() - 2;
		while (ans >= abs_other) {
			while (product[ptr] > ans) {
				ptr--;
			}
			ans -= product[ptr];
		}
		ans._sign = _sign;

		return ans;
	}
	void Integer::operator%=(const Integer& other) {
		*this = *this % other;
	}
	Integer Integer::operator/(const Integer& other) const {
		if (other == _0) {

			throw std::domain_error("除数不能为0");
		}
		if (other == _2) {
			Integer ans = *this;
			bool add = false;
			for (auto i = ans._number.rbegin();i != ans._number.rend();i++) {
				num_t copy = *i;
				*i /= 2;
				if (add) {
					*i += MOD / 2;
				}
				add = (copy % 2 == 1);
			}
			while (ans._number.size() > 1 && ans._number.back() == 0) {
				ans._number.pop_back();
			}

			return ans;
		}
		Integer abs_this = _sign ? -(*this) : *this,
			abs_other = other._sign ? -other : other,
			ans = _0;
		std::vector<Integer> product, power;
		product.push_back(abs_other);
		power.push_back(_1);
		product.push_back(abs_other + abs_other);
		while (product.back() <= abs_this) {
			product.push_back(product.back() + product.back());
			power.push_back(power.back() + power.back());
		}
		std::size_t ptr = product.size() - 2;
		while (abs_this >= abs_other) {
			while (product[ptr] > abs_this) {
				ptr--;
			}
			abs_this -= product[ptr];
			ans += power[ptr];
		}
		ans._sign = _sign ^ other._sign;

		return ans;
	}
	void Integer::operator/=(const Integer& other) {
		*this = *this / other;
	}
	Integer Integer::pow(const Integer& other) const noexcept {
		Integer ans = _1, product = *this;
		ans._sign = _sign;
		for (Integer copy = other;copy > _0;copy /= _2) {
			if (copy % _2 == _1) {
				ans *= product;
			}
			if (copy != _1) {
				product *= product;
			}
		}
		if (other < _0) {

			return (Integer)_1 / ans;
		}

		return ans;
	}
	Integer Integer::pow(const Integer& other, const Integer& mod) const {
		Integer ans = _1, product = *this;
		ans._sign = _sign;
		for (Integer copy = other;copy > _0;copy /= _2) {
			if (copy % _2 == _1) {
				ans *= product;
				ans %= mod;
			}
			if (copy != _1) {
				product *= product;
				product %= mod;
			}
		}
		if (other < _0) {

			return (Integer)_1 / ans;
		}

		return ans;
	}
}


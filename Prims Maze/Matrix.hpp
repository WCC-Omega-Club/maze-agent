
#ifndef matrix_h
#define matrix_h

#include <stdexcept>
#include <vector>
#include <functional>

///<summary>
///undefine to disable range checking
///</summary>
#define RANGE_CHECK
namespace numeric_lib {
	class overdetermined : public std::domain_error
	{
	public:
		overdetermined()
			: std::domain_error("solution is over-determined")
		{}
	};

	class underdetermined : public std::domain_error
	{
	public:
		underdetermined()
			: std::domain_error("solution is under-determined")
		{}
	};


	/**
	* @brief Kahn Summation method
	*/
	template<class T>
	class kahn_sum {
	public:
		kahn_sum() : sum(0.0), cor(0.0) {}
		kahn_sum<T>& operator+=(const T& val) {
			T old_sum = sum;
			T next = val - cor;
			cor = ((sum += next) - old_sum) - next;
			return *this;
		}
		kahn_sum<T>& operator-=(const T& val) {
			T old_sum = sum;
			T next = val + cor;
			cor = ((sum -= val) - old_sum) + next;
			return *this;
		}
		operator T&() { return sum; }
	private:
		T sum;  // running sum
		T cor;  // correction term
	};

	template<class T>
	class matrix {
	private:
		std::vector<T> elements; // array of elements
	public:
	
		const unsigned rows = 0;  // number of rows
		const unsigned cols = 0;  // number of columns

	protected:
		// range check function for matrix access
		void range_check(unsigned i, unsigned j) const;

	public:
		T& operator()(unsigned int i, unsigned int j) {
#ifdef RANGE_CHECK
			range_check(i, j);
#endif
			return elements[i*cols + j];
		}
		const T& operator()(unsigned i, unsigned j) const {
#ifdef RANGE_CHECK
			range_check(i, j);
#endif
			return elements[i*cols + j];
		}
		const T& element(unsigned i, unsigned j) const {
#ifdef RANGE_CHECK
			range_check(i, j);
#endif
			return elements[i*cols + j];
		}
		T& element(unsigned i, unsigned j) {
#ifdef RANGE_CHECK
			range_check(i, j);
#endif
			return elements[i*cols + j];
		}
	public:
		// constructors
		matrix() {};
		matrix(unsigned rows, unsigned cols, const T* elements = 0);
		matrix(const matrix<T>&);
		// destructor
		~matrix();

		// assignment
		matrix<T>& operator=(const matrix<T>&);

		// comparison
		bool operator==(const matrix<T>&) const;
		bool iszero() const;
		bool operator!() const {
			return iszero();
		}

		// scalar multiplication/division
		matrix<T>& operator*=(const T& a);
		matrix<T> operator*(const T& a) const {
			return matrix<T>(*this).operator*=(a);
		}
		matrix<T>& operator/=(const T& a);
		matrix<T> operator/(const T& a) {
			return matrix<T>(*this).operator/=(a);
		}
		matrix<T> operator-() const;
		matrix<T> operator+() const;

		// addition/subtraction
		matrix<T>& operator+=(const matrix<T>&);
		matrix<T>& operator-=(const matrix<T>&);
		matrix<T> operator+(const matrix<T>& M) const {
			return matrix<T>(*this).operator+=(M);
		}
		matrix<T> operator-(const matrix<T>& M) const {
			return matrix<T>(*this).operator-=(M);
		}

		// matrix multiplication
		matrix<T> operator*(const matrix<T>&) const;
		matrix<T>& operator*=(const matrix<T>& M) {
			return *this = *this * M;
		}
		/**
		*
		* @brief boolean matrix operations
		*		 TODO: Add template specialization
		*/
		matrix<T> operator &&(const matrix<T>& B) const;
		matrix<T> boolean_product(const matrix<T>& B) const;
		matrix<T> operator ^(unsigned int power) const {
			matrix<T> M(rows, cols);
			for (unsigned int i = 0; i < power; ++i)
				M = boolean_product(*this);
			return *this;
		}
		matrix<T> operator ||(const matrix<T>& B) const;
		// matrix division
		int size() { return rows *cols; }
		matrix<T> rightdiv(const matrix<T>&) const;
		matrix<T> leftdiv(const matrix<T>& D) const {
			return transpose().rightdiv(D.transpose()).transpose();
		}
		matrix<T> operator/(const matrix<T>& D) const {
			return leftdiv(D);
		}
		matrix<T>& operator/=(const matrix<T>& M) {
			return *this = *this / M;
		}

		// determinants
		matrix<T> minor(unsigned i, unsigned j) const;
		T det() const;
		T minor_det(unsigned i, unsigned j) const;

		// these member functions are only valid for squares
		matrix<T> inverse() const;
		matrix<T> pow(int exp) const;
		matrix<T> identity() const;
		bool isidentity() const;

		// vector operations
		matrix<T> getrow(unsigned j) const;
		matrix<T> getcol(unsigned i) const;
		matrix<T>& setcol(unsigned j, const matrix<T>& C);
		matrix<T>& setrow(unsigned i, const matrix<T>& R);
		matrix<T> delrow(unsigned i) const;
		matrix<T> delcol(unsigned j) const;

		matrix<T> transpose() const;
		matrix<T> operator~() const {
			return transpose();
		}
	};

	template<class T>
	matrix<T>::matrix(unsigned rows, unsigned cols, const T* elements = 0)
		: rows(rows), cols(cols), elements(rows*cols, T(0.0))
	{
		if (rows == 0 || cols == 0)
			throw std::range_error("attempt to create a degenerate matrix");
		// initialze from array
		if (elements)
			for (unsigned i = 0;i < rows*cols;i++)
				this->elements[i] = elements[i];
	};

	template<class T>
	matrix<T>::matrix(const matrix<T>& cp)
		: rows(cp.rows), cols(cp.cols), elements(cp.elements)
	{
	}

	template<class T>
	matrix<T>::~matrix()
	{
	}

	template<class T>
	matrix<T>& matrix<T>::operator=(const matrix<T>& cp)
	{
		if (cp.rows != rows && cp.cols != cols)
			throw std::domain_error("matrix op= not of same order");
		for (unsigned i = 0;i < rows*cols;i++)
			elements[i] = cp.elements[i];
		return *this;
	}


	template<class T>
	void matrix<T>::range_check(unsigned i, unsigned j) const
	{
		if (rows <= i)
			throw std::range_error("matrix access row out of range");
		if (cols <= j)
			throw std::range_error("matrix access col out of range");
	}

	template<class T>
	bool matrix<T>::operator==(const matrix<T>& A) const
	{
		if (A.rows != rows && A.cols != cols)
			throw std::domain_error("matrix op== not of same order");
		for (unsigned i = 0;i < rows;i++)
			for (unsigned j = 0;j < cols;j++)
				if (element(i, j) != A(i, j))
					return false;
		return true;
	}


	/**
	*
	*returns true if this is an additive identity matrix
	*/
	template<class T>
	bool matrix<T>::iszero() const
	{
		const T T0(0.0);  // additive identity for element
		for (unsigned i = 0;i < rows;i++)
			for (unsigned j = 0;j < cols;j++)
				if (element(i, j) != T0)
					return false;
		return true;
	}

	template<class T>
	matrix<T>& matrix<T>::operator*=(const T& a)
	{
		for (unsigned i = 0;i < rows*cols;i++)
			elements[i] *= a;
		return *this;
	}

	template<class T>
	matrix<T>& matrix<T>::operator/=(const T& a)
	{
		for (unsigned i = 0;i < rows*cols;i++)
			elements[i] /= a;
		return *this;
	}

	template<class T>
	matrix<T> matrix<T>::operator-() const
		// unary -
	{
		matrix<T> R(*this);
		for (unsigned i = 0;i < rows*cols;i++)
			R.elements[i] = -R.elements[i];
		return R;
	}

	/**
	*
	*unary +
	*/
	template<class T>
	matrix<T> matrix<T>::operator+() const
	{
		matrix<T> R(*this);
		for (unsigned i = 0;i < rows*cols;i++)
			R.elements[i] = +elements[i];
		return R;
	}

	template<class T>
	inline
		matrix<T> operator*(const T& a, const matrix<T>& M)
	{
		return matrix<T>(M) *= a;
	}

	template<class T>
	matrix<T>& matrix<T>::operator+=(const matrix<T>& M)
	{
		if (cols != M.cols || rows != M.rows)
			throw std::domain_error("op+= matrices must be of same order");
		for (unsigned i = 0;i < rows*cols;i++)
			elements[i] += M.elements[i];
		return *this;
	}

	template<class T>
	matrix<T>& matrix<T>::operator-=(const matrix<T>& M)
	{
		if (cols != M.cols || rows != M.rows)
			throw std::domain_error("op-= matrices must be of same order");
		for (unsigned i = 0;i < rows*cols;i++)
			elements[i] -= M.elements[i];
		return *this;
	}
	/**
	* return the minor matrix of a given element
	* this is the matrix with the column and row of the element deleted
	*/
	template<class T>
	matrix<T> matrix<T>::minor(unsigned i, unsigned j) const
	{
#ifdef RANGE_CHECK
		range_check(i, j);
#endif
		matrix<T>M(rows - 1, cols - 1);
		unsigned i2, j2;
		for (unsigned i1 = 0, i2 = 0; i1 < rows; i1++) {
			if (i != i1) {
				for (unsigned j1 = 0, j2 = 0; j1 < cols; j1++) {
					if (j != j1) {
						M.element(i2, j2) = element(i1, j1);
						j2++;
					}
				}
				i2++;
			}
		}
		return M;
	}

	///<summary>
	/// return the minor matrix of a given element
	/// this is the matrix with the column and row of the element deleted
	///</summary>
	template<class T>
	T matrix<T>::minor_det(unsigned i, unsigned j) const
	{
#ifdef RANGE_CHECK
		range_check(i, j);
#endif
		matrix<T> M(rows - 1, cols - 1);
		unsigned i2, j2;
		for (unsigned i1 = 0, i2 = 0; i1 < rows; i1++) {
			if (i != i1) {
				for (unsigned j1 = 0, j2 = 0; j1 < cols; j1++) {
					if (j != j1) {
						M.element(i2, j2) = element(i1, j1);
						j2++;
					}
				}
				i2++;
			}
		}
		return M.det();
	}

	template<class T>
	T matrix<T>::det() const
	{
		if (cols != rows)
			throw std::domain_error("matrix::det is only valid for square matrices");
		switch (cols) {
		case 1:
			return element(0, 0);

			// use hard-coded calculations for order 2x2 and 3x3 matrices

		case 2: {
			kahn_sum<T> sum;
			sum += element(0, 0)*element(1, 1);
			sum -= element(0, 1)*element(1, 0);
			return sum;
		}

		case 3: {
			kahn_sum<T> sum;
			sum += element(0, 0)*element(1, 1)*element(2, 2);
			sum += element(1, 0)*element(2, 1)*element(0, 2);
			sum += element(2, 0)*element(0, 1)*element(1, 2);
			sum -= element(0, 2)*element(1, 1)*element(2, 0);
			sum -= element(1, 2)*element(2, 1)*element(0, 0);
			sum -= element(2, 2)*element(0, 1)*element(1, 0);
			return sum;
		}

		default: {
			const T T0(0.0);
			kahn_sum<T> sum;
			for (unsigned i = 0;i < cols;i++) {
				T c = element(i, 0);
				if (c != T0) {
					c *= minor_det(i, 0);  // calc det of minor
					if (i % 2)
						sum -= c;
					else
						sum += c;
				}
			}
			return sum;
		}
		}
	}

	template<class T>
	matrix<T> matrix<T>::operator*(const matrix<T>& B) const
	{
		const matrix<T>& A = *this;
		if (A.cols != B.rows)
			throw std::domain_error("matrix multiply: incompatible orders");
		matrix<T> P(A.rows, B.cols);  // product
		for (unsigned j = 0;j < P.cols;j++)
			for (unsigned i = 0;i < P.rows;i++) {
				kahn_sum<T> sum;
				for (unsigned p = 0;p < A.cols;p++)
					sum += A(i, p)*B(p, j);
				P(i, j) = sum;
			}
		return P;
	}

	template<typename T>
	matrix<T> matrix<T>::operator &&(const matrix<T>& B) const
	{
		if (cols != M.cols || rows != M.rows)
			throw std::domain_error("op+= matrices must be of same order");
		for (unsigned i = 0;i < rows*cols;i++)
			elements[i] = elements[i] && B.elements[i];
		return *this;

	}

	template<typename T>
	matrix<T> matrix<T>::operator ||(const matrix<T>& B) const
	{
		if (cols != M.cols || rows != M.rows)
			throw std::domain_error("op+= matrices must be of same order");
		for (unsigned i = 0;i < rows*cols;i++)
			elements[i] = elements[i] || B.elements[i];
		return *this;

	}

	template<class T>
	matrix<T> matrix<T>::boolean_product(const matrix<T>& B) const
	{
		const matrix<T>& A = *this;
		if (A.cols != B.rows)
			throw std::domain_error("matrix boolean product: incompatible orders");
		matrix<T> P(A.rows, B.cols);  // boolean product
		for (unsigned j = 0;j < P.cols;j++)
			for (unsigned i = 0;i < P.rows;i++) {
				kahn_sum<T> sum;
				for (unsigned p = 0;p < A.cols;p++)
					sum = sum || (A(i, p) && B(p, j));
				P(i, j) = sum;
			}
		return P;
	}

	template<class T>
	matrix<T> matrix<T>::rightdiv(const matrix<T>& D) const
	{
		const matrix<T>& N = *this;

		if (N.rows != D.rows)
			throw std::domain_error("matrix divide: incompatible orders");

		matrix<T> Q(D.cols, N.cols);  // quotient matrix

		if (D.rows < D.cols)
			throw underdetermined();

		if (D.rows > D.cols) {
			bool solution = false;
			for (unsigned i = 0;i < D.rows;i++) {
				matrix<T> D2 = D.delrow(i);  // delete a row from the matrix
				matrix<T> N2 = N.delrow(i);
				matrix<T> Q2(Q);
				try {
					Q2 = N2.rightdiv(D2);
				}
				catch (underdetermined x) {
					continue;  // try again with next row
				}
				if (!solution) {
					// this is our possible solution
					solution = true;
					Q = Q2;
				}
				else {
					// do the solutions agree?
					if (Q.rows != Q2.cols)
						throw overdetermined();
				}
			}
			if (!solution)
				throw underdetermined();
			return Q;
		}

		// D.rows == D.cols
		// use Kramer's Rule
		//
		const T T0(0.0); // additive identity

		if (D.cols <= 3) {
			T ddet = D.det();
			if (ddet == T0)
				throw underdetermined();

			for (unsigned i = 0;i < N.cols;i++)
				for (unsigned j = 0;j < D.cols;j++) {
					matrix<T> A(D); // make a copy of the D matrix
									// replace column with numerator vector
					A.setcol(j, N.getcol(i));
					Q(j, i) = A.det() / ddet;
				}
		}
		else {
			// this method optimizes the determinant calculations
			// by saving a minor_det calculations.

			kahn_sum<T> sum;
			std::vector<T> minordet(D.cols);
			for (unsigned j = 0;j < D.cols;j++) {
				T c = D.minor_det(0, j);
				minordet[j] = c;
				T a = D(0, j);
				if (a != T0) {
					a *= c;
					if (j % 2)
						sum -= a;
					else
						sum += a;
				}
			}
			T ddet = sum;
			if (ddet == T0)
				throw underdetermined();
			for (unsigned i = 0;i < N.cols;i++) {
				for (unsigned j = 0;j < D.cols;j++) {
					matrix<T> A(D);
					A.setcol(j, N.getcol(i));
					kahn_sum<T> ndet;
					for (unsigned k = 0;k < D.cols;k++) {
						T a = A(0, k);
						if (a != T0) {
							if (k == j)
								a *= minordet[k];  // use previous calculation
							else
								a *= A.minor_det(0, k); // calculate minor's determinant
							if (k % 2)
								ndet -= a;
							else
								ndet += a;
						}
					}
					Q(j, i) = T(ndet) / ddet;
				}
			}
		}
		return Q;
	}

	// returns the multiplicative inverse of a matrix
	template<class T>
	matrix<T> matrix<T>::inverse() const
	{
		if (cols != rows)
			throw std::domain_error("matrix::inverse is defined for only square matrices");
		return identity() / (*this);
	}

	template<class T>
	matrix<T> matrix<T>::getrow(unsigned i) const
		// returns a vector matrix for row of i of the matrix
	{
#ifdef RANGE_CHECK
		range_check(i, 0);
#endif
		matrix<T> R(1, cols);
		for (unsigned x = 0;x < cols;x++)
			R(0, x) = element(i, x);
		return R;
	}

	template<class T>
	matrix<T> matrix<T>::getcol(unsigned j) const
		// returns a vector matrix for column j
	{
#ifdef RANGE_CHECK
		range_check(0, j);
#endif
		matrix<T> R(rows, 1);
		for (unsigned i = 0;i < rows;i++)
			R(i, 0) = element(i, j);
		return R;
	}

	template<class T>
	matrix<T> matrix<T>::delcol(unsigned j) const
		// returns a copy of the matrix with col j removed
	{
#ifdef RANGE_ERROR
		range_check(0, j);
#endif
		matrix<T> R(rows, cols - 1);
		for (unsigned y = 0;y < rows;y++) {
			unsigned x2 = 0;
			for (unsigned x = 0;x < cols;x++)
				if (x != j) {
					R(y, x2) = element(y, x);
					x2++;
				}
		}
		return R;
	}

	template<class T>
	matrix<T> matrix<T>::delrow(unsigned i) const
		// returns a copy of the matrix with row i removed
	{
#ifdef RANGE_CHECK
		range_check(i, 0);
#endif
		matrix<T> R(rows - 1, cols);
		for (unsigned x = 0;x < cols;x++) {
			unsigned y2 = 0;
			for (unsigned y = 0;y < rows;y++)
				if (y != i) {
					R(y2, x) = element(y, x);
					y2++;
				}
		}
		return R;
	}

	template<class T>
	matrix<T>& matrix<T>::setcol(unsigned j, const matrix<T>& C)
		// set a column of a matrix using a vector matrix
	{
#ifdef RANGE_CHECK
		range_check(0, j);
		if (C.rows != rows || C.cols != 1)
			throw std::range_error("matrix::setcol incompatible matrices");
#endif
		for (unsigned i = 0;i < C.rows;i++)
			element(i, j) = C(i, 0);
		return *this;
	}

	template<class T>
	matrix<T>& matrix<T>::setrow(unsigned i, const matrix<T>& R)
		// set a row of matrix using a vector matrix
	{
#ifdef RANGE_CHECK
		range_check(i, 0);
		if (R.cols != cols || R.rows != 1)
			throw std::range_error("matrix::setrow incompatible matrices");
#endif
		for (unsigned j = 0;j < R.cols;j++)
			element(i, j) = R(0, j);
		return *this;
	}

	template<class T>
	matrix<T> matrix<T>::identity() const
		// creates an identity matrix of the same order as the original matrix
	{
		if (cols != rows)
			throw std::domain_error("matrix::identity is defined only for square matrices");
		const T T1(1.0);
		matrix<T> I(cols, rows);
		// set the diag. to the element multiplicative identity
		for (unsigned i = 0;i < cols;i++)
			I(i, i) = T1;
		return I;
	}

	template<class T>
	bool matrix<T>::isidentity() const
		// returns true if the matrix is the multiplicative identity
	{
		if (cols != rows)
			return false;
		const T T1(1.0);
		const T T0(0.0);
		for (unsigned i = 0;i < rows;i++)
			for (unsigned j = 0;j < cols;j++)
				if (i == j) {
					if (element(i, j) != T1)
						return false;
				}
				else {
					if (element(i, j) != T0)
						return false;
				}
				return true;
	}

	template<class T>
	matrix<T> matrix<T>::pow(int exp) const
		// returns a matrix raised to a power - only valid for square matrices
	{
		if (exp == 0)
			return identity();  // return multiplicative identity matrix
		if (exp == 1)
			return *this;       // return the same matrix
		if (exp == -1)            // return the mulitplicative inverse
			return inverse();
		if (exp < 0)
			return inverse().pow(-exp);
		// exp > 1
		// check for special cases
		// multiplicative or additive identity?
		if (isidentity() || iszero())
			return *this;
		matrix<T> B(*this);
		for (int n = 1;n < exp;n++)
			B *= *this;
		return B;
	}

	template<class T>
	inline
		matrix<T> pow(const matrix<T>& M, int exp)
	{
		return M.pow(exp);
	}

	template<class T>
	matrix<T> matrix<T>::transpose() const
	{
		matrix<T> R(cols, rows);
		for (unsigned i = 0;i < rows;i++)
			for (unsigned j = 0;j < cols;j++)
				R(j, i) = element(i, j);
		return R;
	}

	template<class T>
	std::ostream& operator<<(std::ostream& o, const matrix<T>& M)
	{
		o << endl;
		for (int i = 0;i < M.rows;i++) {
			o << "[ ";
			for (int j = 0;j < M.cols;j++)
				o << M(i, j) << " ";
			o << "]" << endl;
		}
		o << endl;
		return o;
	}


#endif
}
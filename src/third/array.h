#ifndef _ARRAY_H_
#define _ARRAY_H_

#include <cmath>
#include <vector>
#include <iostream>

#undef min
#undef max

class Array
{
private:
    std::vector<double> _data;
    
public:
    // Ctor
    Array(unsigned N = 1, double value = 0.0) : _data(N, value)
    { }
    
    Array(const Array& rhs) : _data(rhs._data)
    { }
    
    // Dtor
    virtual ~Array()
    { }
    
    // Access
    unsigned size() const
    { return _data.size(); }
    
    const double* copy() const
    { return &(_data[0]); }
    
    void copyFrom(unsigned N, const double* const x)
    {
        _data = std::vector<double>(N);
        for (unsigned i = 0; i < N; ++i)
            _data[i] = x[i];
    }
    
    void copyTo(double* const& data)
    {
        for (unsigned i = 0; i < size(); ++i)
            data[i] = _data[i];
    }
    
    const std::vector<double>& getValues() const
    { return _data; }
    
    void setValues(const std::vector<double>& values)
    { _data = values; }
    
    double getValue(unsigned i) const
    { return _data[i]; }
    
    void setValue(unsigned i, const double x)
    { _data[i] = x; }
    
    double& operator [] (unsigned i)
    { return _data[i]; }
    
    const double operator [] (unsigned i) const
    { return _data[i]; }
    
    // Operators
    Array& operator = (const Array& rhs) {
        _data = rhs._data;
        return *this;
    }
    
    Array& operator *= (const double rhs) {
        for (unsigned i=0; i<size(); ++i) _data[i] *= rhs;
        return *this;
    }
    
    Array& operator /= (const double rhs) {
        for (unsigned i=0; i<size(); ++i) _data[i] /= rhs;
        return *this;
    }
    
    Array& operator += (const Array& rhs) {
        for (unsigned i=0; i<size(); ++i) _data[i] += rhs[i];
        return *this;
    }
    
    Array& operator -= (const Array& rhs) {
        for (unsigned i=0; i<size(); ++i) _data[i] -= rhs[i];
        return *this;
    }
    
    // Measures
    double length2() const {
        double sum = 0;
        for (unsigned i=0; i<size(); ++i)
            sum += _data[i]*_data[i];
        return sum;
    }
    
    double length() const
    { return sqrt( length2() ); }
    
    double normalize() {
        double len = length();
        if (len != 0.0)
            *this /= len;
        return len;
    }
    
    double getMax() const {
        double value = _data[0];
        for (unsigned i=1; i<size(); ++i)
            value = std::max(value, _data[i]);
        return value;
    }
    
    double getMin() const {
        double value = _data[0];
        for (unsigned i=1; i<size(); ++i)
            value = std::min(value, _data[i]);
        return value;
    }
    
    double getMean() const {
        double value = 0.0;
        for (unsigned i=0; i < size(); ++i)
            value += _data[i];
        return value / size();
    }
    
    Array slice(unsigned a, unsigned b) const {
        Array x(b-a);
        for (unsigned i = a; i < b; ++i)
            x[i-a] = _data[i];
        return x;
    }
};

inline
bool operator == (const Array& a, const Array& b) {
    for (unsigned i=0; i<a.size(); ++i)
        if (a[i] != b[i]) return false;
    return true;
}

inline
bool operator != (const Array& a, const Array& b)
{ return !(a == b); }

inline
Array operator - (const Array& rhs) {
    Array neg(rhs.size());
    for (unsigned i=0; i<rhs.size(); ++i) neg[i] = - rhs[i];
    return neg;
}

inline
Array operator + (const Array& a, const Array& b) {
    Array rt(a);
    rt += b;
    return rt;
}

inline
Array operator - (const Array& a, const Array& b) {
    Array rt(a);
    rt -= b;
    return rt;
}

inline
Array operator * (const Array& a, const double b) {
    Array rt(a);
    rt *= b;
    return rt;
}

inline
Array operator * (const double a, const Array& b) {
    Array rt(b);
    rt *= a;
    return rt;
}

inline
Array operator / (const Array& a, double b) {
    Array rt(a);
    rt /= b;
    return rt;
}

inline
double dot(const Array& a, const Array& b) {
    double sum = 0.;
    for (unsigned i=0; i<a.size(); ++i) sum += a[i]*b[i];
    return sum;
}

inline
std::ostream& operator << (std::ostream& out, const Array& a) {
    for (unsigned i=0; i<a.size(); ++i) out << a[i] << " ";
    return out;
}

inline
std::istream& operator >> (std::istream& in, Array& a) {
    for (unsigned i=0; i<a.size(); ++i) in >> a[i];
    return in;
}

inline
Array concatenate(const Array& a, const Array& b) {
    Array c(a.size() + b.size());
    for (unsigned i = 0 ; i < a.size(); ++i)
        c[i] = a[i];
    for (unsigned i = 0 ; i < b.size(); ++i)
        c[a.size()+i] = b[i];
    return c;
}

#endif

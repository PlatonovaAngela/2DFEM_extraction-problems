#include "simple_structure.h"

bool point::operator==(point a)
{
    static const double eps = 1.e-6;
    if (fabs(this->x - a.x) < eps && fabs(this->y - a.y) < eps)
        return true;
    else
        return false;
}

point::point(double x, double y)
{
    this->x = x;
    this->y = y;
}

void point::init(double x, double y)
{
    this->x = x;
    this->y = y;
}

line::line()
{
    vert_line = true;
    lp = -1;
    rp = -1;
    nBox1 = -1;
    nBox2 = -1;
}

bool box::is_conductor()
{
    return (eps != eps) ? true : false;
}

int box::get_que()
{
    return que;
}

void box::set_que(int que)
{
    this->que = que;
}

int* box::get_points()
{
    return p;
}

double box::get_a()
{
    return a;
}

void box::set_a(double a)
{
    this->a = a;
}

void box::init(int* points, double eps, double a)
{
    for (int i = 0; i < 4; ++i)
        p[i] = points[i];

    this->eps = eps;

    this->a = a;
}

point operator*(const point& vec, const double& a)
{
    return { vec.x * a, vec.y * a };
}

point operator*(const double& a, const point& vec)
{
    return vec * a;
}

point operator/(const point& vec, const double& a)
{
    return {vec.x / a, vec.y / a};
}

point operator-(const point& vec1, const point& vec2)
{
    return {vec1.x - vec2.x, vec1.y - vec2.y};
}

double dot(const point& vec1, const point& vec2)
{
    return vec1.x * vec2.x + vec1.y * vec2.y;
}

double fabs(const point& vec1)
{
    return sqrt(vec1.x * vec1.x + vec1.y * vec1.y);
}

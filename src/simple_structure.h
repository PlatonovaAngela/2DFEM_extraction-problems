#pragma once
#include <cmath>
#include <vector>

struct point
{
    double x;
    double y;

    bool operator==(point a);

    point(double x = 0.0, double y = 0.0);

    // Для совместимости со старым кодом
    void init(double x, double y);
};

struct box
{
    int p[4];

    double eps;
    
    double a;

    bool is_conductor();

    int que;

    // Для совместимости со старым кодом

    int get_que();

    void set_que(int que);

    int* get_points();

    double get_a();

    void set_a(double a);

    void init(int* points, double eps, double a);
};

struct line
{
    bool vert_line;

    int lp;
    int rp;

    int nBox1;
    int nBox2;

    line();
};


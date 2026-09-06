#pragma once
#include <vector>
#include "simple_structure.h"
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>

class task
{
private:

    std::vector<point> verts; // Вершины в геометрии

    std::vector<box>   boxes; // Боксы в геометрии

    std::vector<line> vert_lines;

    std::vector<line> horz_lines;

public:

    //Чтение файла формата старой версии программы
    int read_project(const std::string& projName);

    // Необязатльно (объединение боксов)
    int fix_boxes();

    // Сохранение геометрии
    int save_boxes_geom(const std::string& geomName);

    // Поиск линий у боксов
    int classify_lines();

    // Тест вывода ребер боксов
    void debug_save_lines();

    // Для построение сетки со сгущением от стенок проводников
    int box_processing(int n = 20, int m = 20);

    box get_box(int id);

    int get_box_size();

    point get_point(int id);

    line get_vert_line(int id);

    line get_horz_line(int id);

    int get_vert_lines_size();

    int get_horz_lines_size();

    double get_max_x();

    double get_min_x();

    double get_max_y();

    double get_min_y();

    int get_nets_size();
};


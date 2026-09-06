#include <iostream>

#include "InitData.h"
#include "task.h"
#include "FEM_solver.h"


int main(int argc, char** argv)
{
    // Чтение аргументов командной строки
    init_data input;
    int err = input.parser(argc, argv);
    if (err != 0)
        return err;

    task cur_task;
    cur_task.read_project(input.task_name);
    
    cur_task.fix_boxes();
    cur_task.box_processing(input.n1, input.n2);

    // Сохранение геометрии (опционально)
    //  g "C:\Users\egora\Desktop\Aspa\Codes\General_progs\geometry.txt"
    if (input.geom_name != "")
    {
        int err = cur_task.save_boxes_geom(input.geom_name);
        if (err < 0)
            return err;
    }

    // Поиск линии боксов
    cur_task.classify_lines();
    //cur_task.debug_save_lines();

    FEM_solver solver;

    err =  solver.make_base_mesh_ver2(cur_task);
    if (err != 0)
       return err;
    if(input.save_mesh != "")
    {
        int err = solver.save_mesh(input.save_mesh);
        if (err < 0)
            return err;
    }
    
    if (input.non_aprox == 1)
        solver.classify_verts_ver1(cur_task);
    else
        solver.classify_verts_ver2(cur_task);

    //solver.debug_mesh_with_classification();
    //solver.debug_approx();

    solver.perform_new_numeration();
    solver.make_system_of_linear_equations();
    

    std::vector<double> diagonals;
    int Netsize = cur_task.get_nets_size();
    for (int i = 0; i < Netsize; ++i)
    {
        solver.make_right(i);
        solver.solve_system();
        diagonals.push_back(solver.get_C());
    }

    std::vector<double> GRD = diagonals;
    for (int i = 0; i < Netsize; ++i)
    {
        std::cout << "C[" << i + 1 << "][" << i + 1 << "] = " << diagonals[i] << "\n";
        for (int j = i + 1; j < Netsize; ++j)
        {
            solver.make_right(i, j);
            solver.solve_system();
            double C = diagonals[i] + diagonals[j] - solver.get_C();
            C *= 0.5;
            std::cout << "C[" << i + 1 << "][" << j + 1 << "] = " << -C << "\n";
            GRD[i] -= C;
            GRD[j] -= C;
        }
        std::cout << "C[" << i + 1 << "][GND] = " << GRD[i] << "\n";
    }

    return 0;
    
}
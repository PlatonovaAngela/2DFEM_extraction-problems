#include "InitData.h"

init_data::init_data()
{
    task_name = "";
    geom_name = "";
    save_mesh = "";
    n1 = 10; 
    n2 = 10;
}

int init_data::parser(std::string field, std::string arg)
{
    //Парсинг help'a
    if (field == full_fields[0] || field == short_fields[0])
    {
        std::cout << full_fields[1] << "\t" << short_fields[1] << "\t" << "Input  geometry file." << std::endl;
        std::cout << full_fields[2] << "\t" << short_fields[2] << "\t" << "Output  geometry file." << std::endl;
        std::cout << full_fields[3] << "\t" << short_fields[3] << "\t" << "Output mesh file." << std::endl;
        std::cout << full_fields[4] << "\t" << short_fields[4] << "\t" << "n1 - number of elements on net edge." << std::endl;
        std::cout << full_fields[5] << "\t" << short_fields[5] << "\t" << "n2 - number of layers in meshing." << std::endl;
        std::cout << full_fields[6] << "\t" << short_fields[6] << "\t" << "Don't approximate vertexes." << std::endl;
        return -2;
    }
    // Парсинг имени файла
    else if(field == full_fields[1] || field == short_fields[1])
    {
        task_name = arg;
        for (int i = 0; i < fields_num; ++i)
        {
            if (full_fields[i] == arg)
            {
                std::cerr << "Argument is key value." << std::endl;
                return -1;
            }

            if (short_fields[i] == arg)
            {
                std::cerr << "Argument is key value." << std::endl;
                return -1;
            }

            if (arg == "")
            {
                std::cerr << "Empty argument input." << std::endl;
                return -1;
            }

            return 2;
        }

    }
    // Парсинг имени файла-геометрии
    else if (field == full_fields[2] || field == short_fields[2])
    {
        geom_name = arg;
        for (int i = 0; i < fields_num; ++i)
        {
            if (full_fields[i] == arg)
            {
                std::cerr << "Argument is key value." << std::endl;
                return -1;
            }

            if (short_fields[i] == arg)
            {
                std::cerr << "Argument is key value." << std::endl;
                return -1;
            }

            if (arg == "")
            {
                std::cerr << "Empty argument input." << std::endl;
                return -1;
            }

            return 2;
        }

    }
    // Парсинг имени файла-сетки
    else if (field == full_fields[3] || field == short_fields[3])
    {
        save_mesh = arg;
        for (int i = 0; i < fields_num; ++i)
        {
            if (full_fields[i] == arg)
            {
                std::cerr << "Argument is key value." << std::endl;
                return -1;
            }

            if (short_fields[i] == arg)
            {
                std::cerr << "Argument is key value." << std::endl;
                return -1;
            }

            if (arg == "")
            {
                std::cerr << "Empty argument input." << std::endl;
                return -1;
            }

            return 2;
        }

    }
    // Парсинг n1
    else if (field == full_fields[4] || field == short_fields[4])
    {
        n1 = atoi(arg.c_str());
        return 2;
    }
    // Парсинг n2
    else if (field == full_fields[5] || field == short_fields[5])
    {
        n2 = atoi(arg.c_str());
        return 2;
    }
    // Парсинг n2
    else if (field == full_fields[6] || field == short_fields[6])
    {
        non_aprox = 1;
        return 1;
    }
    else
        std::cerr << "Unknown argument input." << std::endl;
    return -1;
}

int init_data::parser(int argc, char** argv)
{
    int curarg = 1;
    int shift;
    while (curarg != argc)
    {
        if (curarg + 1 == argc)
        {
            shift = parser(argv[curarg], " ");
        }
        else
        {
            shift = parser(argv[curarg], argv[curarg + 1]);
        }
        if (shift < 0)
            return shift;
        else
            curarg += shift;

    }
    return 0;
}

#pragma once

#include <iostream>
#include <string>

const int fields_num = 8;

const std::string full_fields[] = { "--help", "--task_name", "--geom_name", "--save_mesh", "--n1", "--n2", "--no_aproximate", "--vem_solver"};

const std::string short_fields[] = { "h", "tn", "g", "m", "n1", "n2", "na", "vem"};

struct init_data
{
	std::string task_name;  // Имя задачи
	std::string geom_name;  // Файл для сохранения геометрии
	std::string save_mesh;  // Файл для сохранения сетки
	// Для разбиения сетки (разбиение стороны нета, и количество слоев до бокса с меньшим a)
	int n1;
	int n2;
	// запуск варианта без аппроксимации // теперь номер решателя
	int non_aprox;

	// Конструктор. Начальная инициализация.
	init_data();

	// Парсер для аргумента командной строки
	int parser(std::string field, std::string arg = "");

	// Парсер для всех аргументов 
	int parser(int argc, char** argv);
};


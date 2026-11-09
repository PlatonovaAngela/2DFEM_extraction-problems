#pragma once
#include "simple_structure.h"
#include "task.h"
#include <algorithm>
#include <Eigen/Sparse>
#include <Eigen/Dense>

struct FEM_point: public point
{
	FEM_point(double x = 0.0, double y = 0.0);
	
	bool left;
	bool down;
	bool right;
	bool up;

	bool approx;
	// Вершины в случае аппроксимации
	int vert1;
	int vert2;

	// Принадлежность границу -1 ~ du/dn = 0; -2 ~ u = 0 (Заземление)
	int boundary_ind;

	bool operator==(FEM_point a);

	int local_num;

	bool is_corner();
};

struct FEM_cell
{
	int p[4];

	int group;

	double eps;

	bool left;
	bool down;
	bool right;
	bool up;

	std::vector<int> left_points;
	std::vector<int> down_points;
	std::vector<int> right_points;
	std::vector<int> up_points;

	FEM_cell();
};

class FEM_solver
{
	private:
		std::vector<FEM_cell> cells;

		std::vector<FEM_point> verts;

		int number_of_groups;

		int number_of_independent;

		bool find_local_verts(std::vector<int>& local_verts, int icell);

		void calculate_local_matrix(std::vector<int>& local_verts, int icell, std::vector<std::vector<double>>& local_matrix);
		
		void calculate_local_matrix_ver2(std::vector<int>& local_verts, int icell, std::vector<std::vector<double>>& local_matrix);

		// i к какой степени свободы привязана функция
		double get_dNdx(int i, double a, double b, double x, double y);

		double get_dNdy(int i, double a, double b, double x, double y);

		void get_coeffs_Laplace_interp(double x1, double x2, double xc, double& coef1, double& coef2);

		// Элементы СЛАУ
		std::vector<Eigen::Triplet<double>> LS_Triplets;

		Eigen::VectorXd right;

		std::vector<double> answer;
	public:

		int make_base_mesh(task& problem);

		int add_get_index_point(FEM_point& point);

		int make_base_mesh_ver2(task& problem);

		int save_mesh(std::string filename);

		// Прописать элементы с аппроксимацией только внутри элемента
		int classify_verts_ver1(task& problem);

		// Прописать эелементы с аппроксимацией на границах боксов
		int classify_verts_ver2(task& problem);

		void debug_mesh_with_classification();

		// Создаёт нумерацию без аппроксимированных вершин
		int perform_new_numeration();

		// Создание СЛАУ
		int make_system_of_linear_equations();

		int make_right(int net1, int net2 = -4);

		int solve_system();

		int debug_save_vtk(std::string filename);

		void debug_approx();

		double get_C();

		void debug_save_portrait(std::vector<Eigen::Triplet<double>>& triplets);

		void debug_save_mesh_for_VEM(std::string filename);
		
		void fix_orientation();

		// Решатель VEM
		int vem_make_SLAE();
};

bool point_in_line(point p, task& problem, line curline);
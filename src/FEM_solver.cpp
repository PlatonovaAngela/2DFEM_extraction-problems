#include "FEM_solver.h"

bool FEM_solver::find_local_verts(std::vector<int>& local_verts, int icell)
{
    bool simple_cell = true;
    local_verts.clear();
    for (int i = 0; i < 4; ++i)
    {
        int vert1 = -1;
        int vert2 = -1;

        if (verts[cells[icell].p[i]].approx)
        {
            vert1 = verts[cells[icell].p[i]].vert1;
            vert2 = verts[cells[icell].p[i]].vert2;
            simple_cell = false;
        }
        else
        {
            vert1 = cells[icell].p[i];
        }

        auto it = std::find(local_verts.begin(), local_verts.end(), vert1);
        if (it == local_verts.end()) local_verts.push_back(vert1);

        if (vert2 != -1)
        {
            auto it = std::find(local_verts.begin(), local_verts.end(), vert2);
            if (it == local_verts.end()) local_verts.push_back(vert2);
        }

        switch (i)
        {
            case(0):
                for (int ivert : cells[icell].down_points)
                {
                    auto it = std::find(local_verts.begin(), local_verts.end(), ivert);
                    if (it == local_verts.end()) local_verts.push_back(ivert);
                    simple_cell = false;
                }
                break;
            case(1):
                for (int ivert : cells[icell].right_points)
                {
                    auto it = std::find(local_verts.begin(), local_verts.end(), ivert);
                    if (it == local_verts.end()) local_verts.push_back(ivert);
                    simple_cell = false;
                }
                break;
            case(2):
                for (int ivert : cells[icell].up_points)
                {
                    auto it = std::find(local_verts.begin(), local_verts.end(), ivert);
                    if (it == local_verts.end()) local_verts.push_back(ivert);
                    simple_cell = false;
                }
                break;
            case(3):
                for (int ivert : cells[icell].left_points)
                {
                    auto it = std::find(local_verts.begin(), local_verts.end(), ivert);
                    if (it == local_verts.end()) local_verts.push_back(ivert);
                    simple_cell = false;
                }
                break;
        }
    }

    /*for (int ivert : cells[icell].down_points)
    {
        local_verts.push_back(ivert);
        simple_cell = false;
    }

    for (int ivert : cells[icell].right_points)
    {
        local_verts.push_back(ivert);
        simple_cell = false;
    }

    for (int ivert : cells[icell].up_points)
    {
        local_verts.push_back(ivert);
        simple_cell = false;
    }

    for (int ivert : cells[icell].left_points)
    {
        local_verts.push_back(ivert);
        simple_cell = false;
    }*/
        

    return simple_cell;
}

void FEM_solver::calculate_local_matrix(std::vector<int>& local_verts, int icell, std::vector<std::vector<double>>& local_matrix)
{
    std::vector<double> qxval = { -1.0 / sqrt(3.0), -1.0 / sqrt(3.0),  1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
    std::vector<double> qyval = { -1.0 / sqrt(3.0), 1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
    std::vector<double> cur_qxvals;
    std::vector<double> cur_qyvals;

    cur_qxvals.resize(4);
    cur_qyvals.resize(4);

    int Nx = 2;
    int Ny = 2;

    bool ref_vert_left;
    bool ref_horz_up;

    if (cells[icell].left_points.size() != 0)
    {
        ref_vert_left = true;
        Ny += cells[icell].left_points.size();
    }
    else
    {
        ref_vert_left = false;
        Ny += cells[icell].right_points.size();
    }

    if (cells[icell].up_points.size() != 0)
    {
        ref_horz_up = true;
        Nx += cells[icell].up_points.size();
    }
    else
    {
        ref_horz_up = false;
        Nx += cells[icell].down_points.size();
    }

    std::vector<std::vector<int>> geom_elem_matrix;
    geom_elem_matrix.resize(Ny);
    for (int i = 0; i < Ny; ++i)
    {
        geom_elem_matrix[i].resize(Nx, -1);
    }

    geom_elem_matrix[0][0] = cells[icell].p[3];
    for (int i = 1; i < cells[icell].up_points.size() + 1; ++i)
        geom_elem_matrix[0][i] = cells[icell].up_points[i - 1];
    geom_elem_matrix[0][Nx-1] = cells[icell].p[2];
    for (int i = 1; i < cells[icell].right_points.size() + 1; ++i)
        geom_elem_matrix[Ny - 1 - i][Nx-1] = cells[icell].right_points[Ny - 2  - i];
    geom_elem_matrix[Ny - 1][Nx - 1] = cells[icell].p[1];
    for (int i = 1; i < cells[icell].down_points.size() + 1; ++i)
        geom_elem_matrix[Ny - 1][i] = cells[icell].down_points[i - 1];
    geom_elem_matrix[Ny - 1][0] = cells[icell].p[0];
    for (int i = 1; i < cells[icell].left_points.size() + 1; ++i)
        geom_elem_matrix[Ny - 1 - i][0] = cells[icell].left_points[Ny - 2 - i];


    local_matrix.clear();
    local_matrix.resize(local_verts.size());
    for (int i = 0; i < local_verts.size(); ++i)
        local_matrix[i].resize(local_verts.size(), 0.0);

    for (int i = 0; i < Ny - 1; ++i)
    {
        int vert_up = std::max(geom_elem_matrix[i][0], geom_elem_matrix[i][Nx - 1]);
        int vert_down = std::max(geom_elem_matrix[i + 1][0], geom_elem_matrix[i + 1][Nx - 1]);;
        double b = verts[vert_up].y - verts[vert_down].y;
        b *= 0.5;
        for (int j = 0; j < Nx - 1; ++j)
        {
            int vert_left = std::max(geom_elem_matrix[0][j], geom_elem_matrix[Ny - 1][j]);
            int vert_right = std::max(geom_elem_matrix[0][j + 1], geom_elem_matrix[Ny - 1][j + 1]);;
            double a = verts[vert_right].x - verts[vert_left].x;
            a *= 0.5;

            double S = 4 * a * b;

            for (int k = 0; k < 4; ++k)
            {
                cur_qxvals[k] = qxval[k] * a;
                cur_qyvals[k] = qyval[k] * b;
            }
        
            for (int p1 = 0; p1 < 4; p1++)
            {
                double dxmass1[4];
                double dymass1[4];
                for (int k = 0; k < 4; ++k)
                {
                    dxmass1[k] = get_dNdx(p1, a, b, cur_qxvals[k], cur_qyvals[k]);
                    dymass1[k] = get_dNdy(p1, a, b, cur_qxvals[k], cur_qyvals[k]);
                }

                int i1;
                int j1;

                switch (p1)
                {
                    case(0):
                        i1 = i + 1;
                        j1 = j;
                        break;
                    case(1):
                        i1 = i + 1;
                        j1 = j + 1;
                        break;
                    case(2):
                        i1 = i;
                        j1 = j + 1;
                        break;
                    case(3):
                        i1 = i;
                        j1 = j;
                        break;
                    default:
                        break;
                }

                std::vector<int> dof1;
                std::vector<double> coeffs1;
                
                double x1;
                double y1;
                if (geom_elem_matrix[i1][j1] == -1)
                {
                    if (ref_vert_left)
                        y1 = verts[geom_elem_matrix[i1][0]].y;
                    else
                        y1 = verts[geom_elem_matrix[i1][Nx - 1]].y;

                    if (ref_horz_up)
                        x1 = verts[geom_elem_matrix[0][j1]].x;
                    else
                        x1 = verts[geom_elem_matrix[Ny - 1][j1]].x;
                }
                else
                {
                    x1 = verts[geom_elem_matrix[i1][j1]].x;
                    y1 = verts[geom_elem_matrix[i1][j1]].y;
                }

                if(geom_elem_matrix[i1][j1] == -1)
                {
                    dof1.resize(6);
                    coeffs1.resize(6);
                    if (ref_vert_left)
                    {
                        dof1[0] = geom_elem_matrix[i1][0];
                        dof1[1] = geom_elem_matrix[0][Nx - 1];
                        dof1[2] = geom_elem_matrix[Ny - 1][Nx - 1];
                    }
                    else
                    {
                        dof1[0] = geom_elem_matrix[i1][Nx - 1];
                        dof1[1] = geom_elem_matrix[0][0];
                        dof1[2] = geom_elem_matrix[Ny - 1][0];
                    }

                    double vcoeff1, vcoeff2;
                    get_coeffs_Laplace_interp(verts[dof1[1]].y, verts[dof1[2]].y, y1, vcoeff1, vcoeff2);
                    double hcoeff1, hcoeff2;
                    get_coeffs_Laplace_interp(verts[dof1[0]].x, verts[dof1[2]].x, x1, hcoeff1, hcoeff2);
                    coeffs1[0] = 0.5 * hcoeff1;
                    coeffs1[1] = 0.5 * hcoeff2 * vcoeff1;
                    coeffs1[2] = 0.5 * hcoeff2 * vcoeff2;

                    if (ref_horz_up)
                    {
                        dof1[3] = geom_elem_matrix[0][j1];
                        dof1[4] = geom_elem_matrix[Ny - 1][0];
                        dof1[5] = geom_elem_matrix[Ny - 1][Nx - 1];
                    }
                    else
                    {
                        dof1[3] = geom_elem_matrix[Ny - 1][j1];
                        dof1[4] = geom_elem_matrix[0][0];
                        dof1[5] = geom_elem_matrix[0][Nx - 1];
                    }

                    get_coeffs_Laplace_interp(verts[dof1[3]].y, verts[dof1[4]].y, y1, vcoeff1, vcoeff2);
                    get_coeffs_Laplace_interp(verts[dof1[4]].x, verts[dof1[5]].x, x1, hcoeff1, hcoeff2);
                    coeffs1[3] = 0.5 * vcoeff1;
                    coeffs1[4] = 0.5 * vcoeff2 * hcoeff1;
                    coeffs1[5] = 0.5 * vcoeff2 * hcoeff2;
                }
                else
                {
                    dof1.push_back(geom_elem_matrix[i1][j1]);
                    coeffs1.push_back(1.0);
                }

                std::vector<int> approx_dof1;
                std::vector<double> approx_coeffs1;

                for (int k = 0; k < dof1.size(); ++k)
                {
                    if (verts[dof1[k]].approx)
                    {
                        double coef1, coef2;
                        if (fabs(verts[verts[dof1[k]].vert1].x - verts[verts[dof1[k]].vert2].x) < 1.e-6)
                            get_coeffs_Laplace_interp(verts[verts[dof1[k]].vert1].y, verts[verts[dof1[k]].vert2].y, verts[dof1[k]].y, coef1, coef2);
                        else
                            get_coeffs_Laplace_interp(verts[verts[dof1[k]].vert1].x, verts[verts[dof1[k]].vert2].x, verts[dof1[k]].x, coef1, coef2);
                        approx_dof1.push_back(verts[dof1[k]].vert1);
                        approx_dof1.push_back(verts[dof1[k]].vert2);
                        approx_coeffs1.push_back(coeffs1[k] * coef1);
                        approx_coeffs1.push_back(coeffs1[k] * coef2);
                    }
                    else
                    {
                        approx_dof1.push_back(dof1[k]);
                        approx_coeffs1.push_back(coeffs1[k]);
                    }
                }

                for (int p2 = 0; p2 < 4; p2++)
                {
                    double dxmass2[4];
                    double dymass2[4];
                    for (int k = 0; k < 4; ++k)
                    {
                        dxmass2[k] = get_dNdx(p2, a, b, cur_qxvals[k], cur_qyvals[k]);
                        dymass2[k] = get_dNdy(p2, a, b, cur_qxvals[k], cur_qyvals[k]);
                    }

                    int i2;
                    int j2;

                    switch (p2)
                    {
                    case(0):
                        i2 = i + 1;
                        j2 = j;
                        break;
                    case(1):
                        i2 = i + 1;
                        j2 = j + 1;
                        break;
                    case(2):
                        i2 = i;
                        j2 = j + 1;
                        break;
                    case(3):
                        i2 = i;
                        j2 = j;
                        break;
                    default:
                        break;
                    }

                    std::vector<int> dof2;
                    std::vector<double> coeffs2;

                    double x2;
                    double y2;
                    if (geom_elem_matrix[i2][j2] == -1)
                    {
                        if (ref_vert_left)
                            y2 = verts[geom_elem_matrix[i2][0]].y;
                        else
                            y2 = verts[geom_elem_matrix[i2][Nx - 1]].y;

                        if (ref_horz_up)
                            x2 = verts[geom_elem_matrix[0][j2]].x;
                        else
                            x2 = verts[geom_elem_matrix[Ny - 1][j2]].x;
                    }
                    else
                    {
                        x2 = verts[geom_elem_matrix[i2][j2]].x;
                        y2 = verts[geom_elem_matrix[i2][j2]].y;
                    }

                    if (geom_elem_matrix[i2][j2] == -1)
                    {
                        dof2.resize(6);
                        coeffs2.resize(6);
                        if (ref_vert_left)
                        {
                            dof2[0] = geom_elem_matrix[i2][0];
                            dof2[1] = geom_elem_matrix[0][Nx - 1];
                            dof2[2] = geom_elem_matrix[Ny - 1][Nx - 1];
                        }
                        else
                        {
                            dof2[0] = geom_elem_matrix[i2][Nx - 1];
                            dof2[1] = geom_elem_matrix[0][0];
                            dof2[2] = geom_elem_matrix[Ny - 1][0];
                        }

                        double vcoeff1, vcoeff2;
                        get_coeffs_Laplace_interp(verts[dof2[1]].y, verts[dof2[2]].y, y2, vcoeff1, vcoeff2);
                        double hcoeff1, hcoeff2;
                        get_coeffs_Laplace_interp(verts[dof2[0]].x, verts[dof2[2]].x, x2, hcoeff1, hcoeff2);
                        coeffs2[0] = 0.5 * hcoeff1;
                        coeffs2[1] = 0.5 * hcoeff2 * vcoeff1;
                        coeffs2[2] = 0.5 * hcoeff2 * vcoeff2;

                        if (ref_horz_up)
                        {
                            dof2[3] = geom_elem_matrix[0][j2];
                            dof2[4] = geom_elem_matrix[Ny - 1][0];
                            dof2[5] = geom_elem_matrix[Ny - 1][Nx - 1];
                        }
                        else
                        {
                            dof2[3] = geom_elem_matrix[Ny - 1][j2];
                            dof2[4] = geom_elem_matrix[0][0];
                            dof2[5] = geom_elem_matrix[0][Nx - 1];
                        }

                        get_coeffs_Laplace_interp(verts[dof2[3]].y, verts[dof2[4]].y, y2, vcoeff1, vcoeff2);
                        get_coeffs_Laplace_interp(verts[dof2[4]].x, verts[dof2[5]].x, x2, hcoeff1, hcoeff2);
                        coeffs2[3] = 0.5 * vcoeff1;
                        coeffs2[4] = 0.5 * vcoeff2 * hcoeff1;
                        coeffs2[5] = 0.5 * vcoeff2 * hcoeff2;
                    }
                    else
                    {
                        dof2.push_back(geom_elem_matrix[i2][j2]);
                        coeffs2.push_back(1.0);
                    }

                    std::vector<int> approx_dof2;
                    std::vector<double> approx_coeffs2;

                    for (int k = 0; k < dof2.size(); ++k)
                    {
                        if (verts[dof2[k]].approx)
                        {
                            double coef1, coef2;
                            if (fabs(verts[verts[dof2[k]].vert1].x - verts[verts[dof2[k]].vert2].x) < 1.e-6)
                                get_coeffs_Laplace_interp(verts[verts[dof2[k]].vert1].y, verts[verts[dof2[k]].vert2].y, verts[dof2[k]].y, coef1, coef2);
                            else
                                get_coeffs_Laplace_interp(verts[verts[dof2[k]].vert1].x, verts[verts[dof2[k]].vert2].x, verts[dof2[k]].x, coef1, coef2);
                            approx_dof2.push_back(verts[dof2[k]].vert1);
                            approx_dof2.push_back(verts[dof2[k]].vert2);
                            approx_coeffs2.push_back(coeffs2[k] * coef1);
                            approx_coeffs2.push_back(coeffs2[k] * coef2);
                        }
                        else
                        {
                            approx_dof2.push_back(dof2[k]);
                            approx_coeffs2.push_back(coeffs2[k]);
                        }
                    }
                
                    for (int idof1 = 0; idof1 < approx_dof1.size(); idof1++)
                    {
                        auto it1 = std::find(local_verts.begin(), local_verts.end(), approx_dof1[idof1]);
                        int ind1 = it1 - local_verts.begin();
                        for (int idof2 = 0; idof2 < approx_dof2.size(); idof2++)
                        {
                            auto it2 = std::find(local_verts.begin(), local_verts.end(), approx_dof2[idof2]);
                            int ind2 = it2 - local_verts.begin();
                            for (int k = 0; k < 4; ++k)
                                local_matrix[ind1][ind2] += approx_coeffs1[idof1] * approx_coeffs2[idof2] * (dxmass1[k] * dxmass2[k] + dymass1[k] * dymass2[k]) * S * 0.25;
                        }
                    }
                }
            }
        }
    }

    for (int i = 0; i < local_verts.size(); ++i)
    {
        for (int j = 0; j < local_verts.size(); ++j)
        {
            local_matrix[i][j] *= cells[icell].eps;
        }
    }
}

void FEM_solver::calculate_local_matrix_ver2(std::vector<int>& local_verts, int icell, std::vector<std::vector<double>>& local_matrix)
{
    std::vector<double> qxval = { -1.0 / sqrt(3.0), -1.0 / sqrt(3.0),  1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
    std::vector<double> qyval = { -1.0 / sqrt(3.0), 1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
    std::vector<double> cur_qxvals;
    std::vector<double> cur_qyvals;

    cur_qxvals.resize(4);
    cur_qyvals.resize(4);

    int Nx = 2;
    int Ny = 2;

    bool ref_vert_left;
    bool ref_horz_up;

    if (cells[icell].left_points.size() != 0)
    {
        ref_vert_left = true;
        Ny += cells[icell].left_points.size();
    }
    else
    {
        ref_vert_left = false;
        Ny += cells[icell].right_points.size();
    }

    if (cells[icell].up_points.size() != 0)
    {
        ref_horz_up = true;
        Nx += cells[icell].up_points.size();
    }
    else
    {
        ref_horz_up = false;
        Nx += cells[icell].down_points.size();
    }

    std::vector<std::vector<int>> geom_elem_matrix;
    geom_elem_matrix.resize(Ny);
    for (int i = 0; i < Ny; ++i)
    {
        geom_elem_matrix[i].resize(Nx, -1);
    }

    geom_elem_matrix[0][0] = cells[icell].p[3];
    for (int i = 1; i < cells[icell].up_points.size() + 1; ++i)
        geom_elem_matrix[0][i] = cells[icell].up_points[i - 1];
    geom_elem_matrix[0][Nx - 1] = cells[icell].p[2];
    for (int i = 1; i < cells[icell].right_points.size() + 1; ++i)
        geom_elem_matrix[Ny - 1 - i][Nx - 1] = cells[icell].right_points[Ny - 2 - i];
    geom_elem_matrix[Ny - 1][Nx - 1] = cells[icell].p[1];
    for (int i = 1; i < cells[icell].down_points.size() + 1; ++i)
        geom_elem_matrix[Ny - 1][i] = cells[icell].down_points[i - 1];
    geom_elem_matrix[Ny - 1][0] = cells[icell].p[0];
    for (int i = 1; i < cells[icell].left_points.size() + 1; ++i)
        geom_elem_matrix[Ny - 1 - i][0] = cells[icell].left_points[Ny - 2 - i];


    local_matrix.clear();
    local_matrix.resize(local_verts.size());
    for (int i = 0; i < local_verts.size(); ++i)
        local_matrix[i].resize(local_verts.size(), 0.0);

    // Рассмотр внутренних квадратов
    for (int i = 0; i < Ny - 1; ++i)
    {
        int vert_up = std::max(geom_elem_matrix[i][0], geom_elem_matrix[i][Nx - 1]);
        int vert_down = std::max(geom_elem_matrix[i + 1][0], geom_elem_matrix[i + 1][Nx - 1]);;
        double b = verts[vert_up].y - verts[vert_down].y;
        double yc = (verts[vert_up].y + verts[vert_down].y) * 0.5;
        b *= 0.5;
        for (int j = 0; j < Nx - 1; ++j)
        {
            int vert_left = std::max(geom_elem_matrix[0][j], geom_elem_matrix[Ny - 1][j]);
            int vert_right = std::max(geom_elem_matrix[0][j + 1], geom_elem_matrix[Ny - 1][j + 1]);;
            double a = verts[vert_right].x - verts[vert_left].x;
            a *= 0.5;
            double xc = (verts[vert_right].x + verts[vert_left].x) * 0.5;
            double S = 4 * a * b;
            
            for (int qind = 0; qind < 4; ++qind)
            {
                std::vector<int> indexes;
                std::vector<double> dx;
                std::vector<double> dy;

                for (int p1 = 0; p1 < 4; p1++)
                {
                    int i1;
                    int j1;

                    switch (p1)
                    {
                    case(0):
                        i1 = i + 1;
                        j1 = j;
                        break;
                    case(1):
                        i1 = i + 1;
                        j1 = j + 1;
                        break;
                    case(2):
                        i1 = i;
                        j1 = j + 1;
                        break;
                    case(3):
                        i1 = i;
                        j1 = j;
                        break;
                    default:
                        break;
                    }
                
                    double x1;
                    double y1;
                    if (geom_elem_matrix[i1][j1] == -1)
                    {
                        if (ref_vert_left)
                            y1 = verts[geom_elem_matrix[i1][0]].y;
                        else
                            y1 = verts[geom_elem_matrix[i1][Nx - 1]].y;

                        if (ref_horz_up)
                            x1 = verts[geom_elem_matrix[0][j1]].x;
                        else
                            x1 = verts[geom_elem_matrix[Ny - 1][j1]].x;
                    }
                    else
                    {
                        x1 = verts[geom_elem_matrix[i1][j1]].x;
                        y1 = verts[geom_elem_matrix[i1][j1]].y;
                    }
                
                    double cdx = get_dNdx(p1, a, b, a * qxval[qind], b * qyval[qind]);
                    double cdy = get_dNdy(p1, a, b, a * qxval[qind], b * qyval[qind]);

                    std::vector<int> dof1;
                    std::vector<double> coeffs1;


                    if (geom_elem_matrix[i1][j1] == -1)
                    {
                        double vcoeff1, vcoeff2;
                        double hcoeff1, hcoeff2;
                        dof1.resize(6);
                        coeffs1.resize(6);

                        if (ref_horz_up)
                        {
                            dof1[2] = geom_elem_matrix[0][j1];
                            dof1[0] = geom_elem_matrix[Ny - 1][0];
                            dof1[1] = geom_elem_matrix[Ny - 1][Nx - 1];
                        }
                        else
                        {
                            dof1[2] = geom_elem_matrix[Ny - 1][j1];
                            dof1[0] = geom_elem_matrix[0][0];
                            dof1[1] = geom_elem_matrix[0][Nx - 1];
                        }

                        get_coeffs_Laplace_interp(verts[dof1[2]].y, verts[dof1[0]].y, y1, vcoeff1, vcoeff2);
                        get_coeffs_Laplace_interp(verts[dof1[0]].x, verts[dof1[1]].x, x1, hcoeff1, hcoeff2);
                        coeffs1[2] = 0.5 * vcoeff1;
                        coeffs1[0] = 0.5 * vcoeff2 * hcoeff1;
                        coeffs1[1] = 0.5 * vcoeff2 * hcoeff2;

                        if (ref_vert_left)
                        {
                            dof1[5] = geom_elem_matrix[i1][0];
                            dof1[3] = geom_elem_matrix[0][Nx - 1];
                            dof1[4] = geom_elem_matrix[Ny - 1][Nx - 1];
                        }
                        else
                        {
                            dof1[5] = geom_elem_matrix[i1][Nx - 1];
                            dof1[3] = geom_elem_matrix[0][0];
                            dof1[4] = geom_elem_matrix[Ny - 1][0];
                        }

                        
                        get_coeffs_Laplace_interp(verts[dof1[3]].y, verts[dof1[4]].y, y1, vcoeff1, vcoeff2);
                        get_coeffs_Laplace_interp(verts[dof1[5]].x, verts[dof1[3]].x, x1, hcoeff1, hcoeff2);
                        coeffs1[5] = 0.5 * hcoeff1;
                        coeffs1[3] = 0.5 * hcoeff2 * vcoeff1;
                        coeffs1[4] = 0.5 * hcoeff2 * vcoeff2;
                    }
                    else
                    {
                        dof1.push_back(geom_elem_matrix[i1][j1]);
                        coeffs1.push_back(1.0);
                    }
                
                    for (int idof = 0; idof < dof1.size(); ++idof)
                    {
                        indexes.push_back(dof1[idof]);
                        dx.push_back(coeffs1[idof] * cdx);
                        dy.push_back(coeffs1[idof] * cdy);
                    }
                }
                
                std::vector<double> Locdx;
                std::vector<double> Locdy;
                std::vector<int> Locindexes;

                for (int k = 0; k < indexes.size(); ++k)
                {
                    if (verts[indexes[k]].approx)
                    {
                        double coef1, coef2;
                        if (fabs(verts[verts[indexes[k]].vert1].x - verts[verts[indexes[k]].vert2].x) < 1.e-6)
                            get_coeffs_Laplace_interp(verts[verts[indexes[k]].vert1].y, verts[verts[indexes[k]].vert2].y, verts[indexes[k]].y, coef1, coef2);
                        else
                            get_coeffs_Laplace_interp(verts[verts[indexes[k]].vert1].x, verts[verts[indexes[k]].vert2].x, verts[indexes[k]].x, coef1, coef2);
                        auto it = std::find(local_verts.begin(), local_verts.end(), verts[indexes[k]].vert1);
                        Locindexes.push_back(it - local_verts.begin());
                        it = std::find(local_verts.begin(), local_verts.end(), verts[indexes[k]].vert2);
                        Locindexes.push_back(it - local_verts.begin());
                        Locdx.push_back(dx[k] * coef1);
                        Locdx.push_back(dx[k] * coef2);
                        Locdy.push_back(dy[k] * coef1);
                        Locdy.push_back(dy[k] * coef2);
                    }
                    else
                    {
                        auto it = std::find(local_verts.begin(), local_verts.end(), indexes[k]);
                        Locindexes.push_back(it - local_verts.begin());
                        Locdx.push_back(dx[k]);
                        Locdy.push_back(dy[k]);
                    }
                }
            
                for (int ival1 = 0; ival1 < Locindexes.size(); ++ival1)
                {
                    int ind1 = Locindexes[ival1];
                    for (int ival2 = 0; ival2 < Locindexes.size(); ++ival2)
                    {
                        int ind2 = Locindexes[ival2];
                        local_matrix[ind1][ind2] += (a * b) * (Locdx[ival1] * Locdx[ival2] + Locdy[ival1] * Locdy[ival2]);
                    }
                }
            }
        }
    }

    for (int i = 0; i < local_verts.size(); ++i)
    {
        for (int j = 0; j < local_verts.size(); ++j)
        {
            local_matrix[i][j] *= cells[icell].eps;
        }
    }

}

double FEM_solver::get_dNdx(int i, double a, double b, double x, double y)
{
    switch (i)
    {
        case(0):
            return -(b - y) / (4 * a * b);
        case(1):
            return  (b - y) / (4 * a * b);
        case(2):
            return  (b + y) / (4 * a * b);
        case(3):
            return -(b + y) / (4 * a * b);
        default:
            break;
    }
    return 0.0;
}

double FEM_solver::get_dNdy(int i, double a, double b, double x, double y)
{
    switch (i)
    {
        case(0):
            return -(a - x) / (4 * a * b);
        case(1):
            return -(a + x) / (4 * a * b);
        case(2):
            return  (a + x) / (4 * a * b);
        case(3):
            return  (a - x) / (4 * a * b);
        default:
            break;
    }
    return 0.0;
}

void FEM_solver::get_coeffs_Laplace_interp(double x1, double x2, double xc, double& coef1, double& coef2)
{
    coef1 = (xc - x2) / (x1 - x2);
    coef2 = (xc - x1) / (x2 - x1);
}

int FEM_solver::make_base_mesh(task& problem)
{
    number_of_groups = problem.get_box_size();

    std::vector<int> vertp;
    std::vector<int> horzp;

    auto func = [&](int ind, FEM_point vert)
    {
        return (verts[ind] == vert) ? true : false;
    };

    for (int ibox = 0; ibox < problem.get_box_size(); ++ibox)
    {
        box curbox = problem.get_box(ibox);
        if (curbox.is_conductor()) continue;

        double xmin, xmax;
        double ymin, ymax;

        xmin = problem.get_point(curbox.p[0]).x;
        xmax = problem.get_point(curbox.p[1]).x;
        ymin = problem.get_point(curbox.p[0]).y;
        ymax = problem.get_point(curbox.p[2]).y;

        int Nx = (xmax - xmin) / curbox.a;
        if (fabs(Nx * curbox.a - (xmax - xmin)) > 1.e-6)
            ++Nx;
        if (Nx == 1)
            Nx++;


        int Ny = (ymax - ymin) / curbox.a;
        if (fabs(Ny * curbox.a - (ymax - ymin)) > 1.e-6)
            ++Ny;
        if (Ny == 1)
            Ny++;

        double hx = (xmax - xmin) / Nx;
        double hy = (ymax - ymin) / Ny;
        for (int i = 0; i < Ny; ++i)
        {
            for (int j = 0; j < Nx; ++j)
            {
                FEM_cell curcell;


                int id;
                double x;
                double y;
                FEM_point curpoint;
                bool fl1 = false;
                bool fl2 = false;

                // Левая нижняя точка
                if (j != 0)
                {
                    id = cells.back().p[1];
                    curcell.p[0] = id;
                }
                else
                {
                    x = xmin + hx * j;
                    y = ymax - hy * (i + 1);
                    curpoint = FEM_point(x, y);
                    id = -1;
                    fl1 = false;
                    fl2 = false;
                    auto it = vertp.begin();
                    for (; it != vertp.end(); ++it)
                    {
                        if (func(*it, curpoint))
                        {
                            id = *it;
                            fl1 = true;
                            break;
                        }
                    }

                    if (i + 1 == Ny)
                    {
                        auto it = horzp.begin();
                        for (; it != horzp.end(); ++it)
                        {
                            if (func(*it, curpoint))
                            {
                                id = *it;
                                fl2 = true;
                                break;
                            }
                        }
                    }
                    if (id == -1)
                    {
                        id = verts.size();
                        verts.push_back(curpoint);
                    }
                    if (j == 0 && !fl1)
                        vertp.push_back(id);
                    if (i + 1 == Ny && !fl2)
                        horzp.push_back(id);
                    curcell.p[0] = id;
                    if (j == 0)
                    {
                        verts[id].left = true;
                        curcell.left = true;
                    }
                    if (i + 1 == Ny)
                    {
                        verts[id].down = true;
                        curcell.down = true;
                    }
                }
                    
                // Правая нижняя точка
                x = xmin + hx * (j + 1);
                y = ymax - hy * (i + 1);
                curpoint = FEM_point(x,y);
                id = -1;
                fl1 = false;
                fl2 = false;
                if (j + 1 == Nx)
                {
                    auto it = vertp.begin();
                    for (; it != vertp.end(); ++it)
                    {
                        if (func(*it, curpoint))
                        {
                            id = *it;
                            fl1 = true;
                            break;
                        }
                    }

                }
                if (i + 1 == Ny)
                {
                    auto it = horzp.begin();
                    for (; it != horzp.end(); ++it)
                    {
                        if (func(*it, curpoint))
                        {
                            id = *it;
                            fl2 = true;
                            break;
                        }
                    }
                }
                if (id == -1)
                {
                    id = verts.size();
                    verts.push_back(curpoint);
                }
                if (j + 1 == Nx && !fl1)
                    vertp.push_back(id);
                if (i + 1 == Ny && !fl2)
                    horzp.push_back(id);
                curcell.p[1] = id;
                if (j + 1 == Nx)
                {
                    verts[id].right = true;
                    curcell.right = true;
                }
                if (i + 1 == Ny)
                {
                    verts[id].down = true;
                    curcell.down = true;
                }

                // Правая верхняя точка
                if (i == 0)
                {
                    x = xmin + hx * (j + 1);
                    y = ymax - hy * i;
                    curpoint = FEM_point(x, y);
                    id = -1;
                    fl1 = false;
                    fl2 = false;
                    if (j + 1 == Nx)
                    {
                        auto it = vertp.begin();
                        for (; it != vertp.end(); ++it)
                        {
                            if (func(*it, curpoint))
                            {
                                id = *it;
                                fl1 = true;
                                break;
                            }
                        }

                    }
                    auto it = horzp.begin();
                    for (; it != horzp.end(); ++it)
                    {
                            if (func(*it, curpoint))
                            {
                                id = *it;
                                fl2 = true;
                                break;
                            }
                    }
                    
                    if (id == -1)
                    {
                        id = verts.size();
                        verts.push_back(curpoint);
                    }
                    if (j + 1 == Nx && !fl1)
                        vertp.push_back(id);
                    if (i == 0 && !fl2)
                        horzp.push_back(id);
                    curcell.p[2] = id;
                    if (j + 1 == Nx)
                    {
                        verts[id].right = true;
                        curcell.right = true;
                    }
                    if (i == 0)
                    {
                        verts[id].up = true;
                        curcell.up = true;
                    }
                }
                else
                {
                    id = cells[cells.size() - Nx].p[1];
                    curcell.p[2] = id;
                }

                // Левая верхняя точка
                if (i == 0)
                {   
                    x = xmin + hx * j;
                    y = ymax - hy * i;
                    curpoint = FEM_point(x, y);
                    id = -1;
                    fl1 = false;
                    fl2 = false;
                    if (j == 0)
                    {
                        auto it = vertp.begin();
                        for (; it != vertp.end(); ++it)
                        {
                            if (func(*it, curpoint))
                            {
                                id = *it;
                                fl1 = true;
                                break;
                            }
                        }

                    }
                    if (i == 0)
                    {
                        auto it = horzp.begin();
                        for (; it != horzp.end(); ++it)
                        {
                            if (func(*it, curpoint))
                            {
                                id = *it;
                                fl2 = true;
                                break;
                            }
                        }
                    }
                    if (id == -1)
                    {
                        id = verts.size();
                        verts.push_back(curpoint);
                    }
                    if (j == 0 && !fl1)
                        vertp.push_back(id);
                    if (i == 0 && !fl2)
                        horzp.push_back(id);
                    curcell.p[3] = id;
                    if (j == 0)
                    {
                        verts[id].left = true;
                        curcell.left = true;
                    }
                    if (i == 0)
                    {
                        verts[id].up = true;
                        curcell.up = true;
                    }
                }
                else
                {
                    id = cells[cells.size() - Nx].p[0];
                    curcell.p[3] = id;
                }
                    

                curcell.group = ibox;
                curcell.eps = problem.get_box(ibox).eps;
                cells.push_back(curcell);
            }
        }
    }
    return 0;
}

int FEM_solver::add_get_index_point(FEM_point& point)
{
    auto it = std::find(verts.begin(), verts.end(), point);
    if (it == verts.end())
    {
        verts.push_back(point);
        return verts.size() - 1;
    }
    else
    {
        return it - verts.begin();
    }
}

int FEM_solver::make_base_mesh_ver2(task& problem)
{
    number_of_groups = problem.get_box_size();

    for (int icell = 0; icell < problem.get_box_size(); ++icell)
    {
        box area = problem.get_box(icell);
        if (area.is_conductor()) continue;

        double xmin = problem.get_point(area.p[0]).x;
        double xmax = problem.get_point(area.p[1]).x;
        double ymin = problem.get_point(area.p[0]).y;
        double ymax = problem.get_point(area.p[2]).y;

        int Nx = int((xmax - xmin) / area.get_a());
        int Ny = int((ymax - ymin) / area.get_a());
        if (fabs(Nx * area.get_a() - (xmax - xmin)) > 1.e-6)
            ++Nx;
        if (fabs(Ny * area.get_a() - (ymax - ymin)) > 1.e-6)
            ++Ny;

        if (Nx == 1) Nx++;
        if (Ny == 1) Ny++;

        double* xcoords = new double[Nx + 1];
        double* ycoords = new double[Ny + 1];

        double hx = (xmax - xmin) / Nx;
        ++Nx;
        xcoords[0] = xmin;
        for (int ivert = 1; ivert < Nx; ++ivert)
        {
            xcoords[ivert] = xcoords[ivert - 1] + hx;
        }
        if (xcoords[Nx - 1] > xmax)
            xcoords[Nx - 1] = xmax;

        double hy = (ymax - ymin) / Ny;
        ++Ny;
        ycoords[0] = ymax;
        for (int ivert = 1; ivert < Ny; ++ivert)
        {
            ycoords[ivert] = ycoords[ivert - 1] - hy;
        }
        if (ycoords[Ny - 1] < ymin)
            ycoords[Ny - 1] = ymin;

        --Nx;
        --Ny;

        for (int iyvert = 0; iyvert < Ny; ++iyvert)
        {
            for (int ixvert = 0; ixvert < Nx; ++ixvert)
            {

                FEM_cell area;

                // Левая нижняя точка
                FEM_point v0;
                v0.x = xcoords[ixvert];
                v0.y = ycoords[iyvert + 1];
                int ind = add_get_index_point(v0);
                if (iyvert + 1 == Ny)
                    verts[ind].down = true;
                if (ixvert == 0)
                    verts[ind].left = true;
                area.p[0] = ind;

                // Правая нижняя точка
                v0 = FEM_point();
                v0.x = xcoords[ixvert + 1];
                v0.y = ycoords[iyvert + 1];
                ind = add_get_index_point(v0);
                if (iyvert + 1 == Ny)
                    verts[ind].down = true;
                if (ixvert + 1 == Nx)
                    verts[ind].right = true;
                area.p[1] = ind;

                // Правая верхняя точка
                v0 = FEM_point();
                v0.x = xcoords[ixvert + 1];
                v0.y = ycoords[iyvert];
                ind = add_get_index_point(v0);
                if (iyvert == 0)
                    verts[ind].up = true;
                if (ixvert + 1 == Nx)
                    verts[ind].right = true;
                area.p[2] = ind;

                //левая верхняя точка
                v0 = FEM_point();
                v0.x = xcoords[ixvert];
                v0.y = ycoords[iyvert];
                ind = add_get_index_point(v0);
                if (iyvert == 0)
                    verts[ind].up = true;
                if (ixvert == 0)
                    verts[ind].left = true;
                area.p[3] = ind;

                if (iyvert == 0)
                    area.up = true;
                if (iyvert + 1 == Ny)
                    area.down = true;
                if (ixvert == 0)
                    area.left = true;
                if (ixvert + 1 == Nx)
                    area.right = true;

                area.group = icell;
                area.eps = problem.get_box(icell).eps;
                cells.push_back(area);
            }
        }

        delete[] xcoords;
        delete[] ycoords;
    }

    return 0;
}

int FEM_solver::save_mesh(std::string filename)
{
    std::ofstream file(filename);
    if (!file.is_open())
        return 1;
    file << cells.size() << '\n';
    for (int icell = 0; icell < cells.size(); ++icell) 
        file << cells[icell].p[0] << " " << cells[icell].p[1] << " " << cells[icell].p[2] << " " << cells[icell].p[3] << "\n";
    file << verts.size() << '\n';
    for (int ivert = 0; ivert < verts.size(); ++ivert)
        file << verts[ivert].x << " " << verts[ivert].y << "\n";
    file.close();

    return 0;
}

int FEM_solver::classify_verts_ver1(task& problem)
{   
    int* left_cells_sizes = new int[number_of_groups + 1];
    int* down_cells_sizes = new int[number_of_groups + 1];
    int* right_cells_sizes = new int[number_of_groups + 1];
    int* up_cells_sizes = new int[number_of_groups + 1];

    for (int i = 0; i < number_of_groups + 1; ++i)
    {
        left_cells_sizes[i] = 0;
        down_cells_sizes[i] = 0;
        right_cells_sizes[i] = 0;
        up_cells_sizes[i] = 0;
    }

    for (int icell = 0; icell < cells.size(); ++icell)
    {
        if(cells[icell].left)
            left_cells_sizes[cells[icell].group + 1]++;
        if (cells[icell].down)
            down_cells_sizes[cells[icell].group + 1]++;
        if (cells[icell].right)
            right_cells_sizes[cells[icell].group + 1]++;
        if (cells[icell].up)
            up_cells_sizes[cells[icell].group + 1]++;
    }

    for (int i = 0; i < number_of_groups; ++i)
    {
        left_cells_sizes[i + 1] += left_cells_sizes[i];
        down_cells_sizes[i + 1] += down_cells_sizes[i];
        right_cells_sizes[i + 1] += right_cells_sizes[i];
        up_cells_sizes[i + 1] += up_cells_sizes[i];
    }

    int cur_left = 0;
    int cur_down = 0;
    int cur_right = 0;
    int cur_up = 0;
    int* left_cells = new int[left_cells_sizes[number_of_groups]];
    int* down_cells = new int[down_cells_sizes[number_of_groups]];
    int* right_cells = new int[right_cells_sizes[number_of_groups]];
    int* up_cells = new int[up_cells_sizes[number_of_groups]];

    for (int icell = 0; icell < cells.size(); ++icell)
    {
        if (cells[icell].left)
            left_cells[cur_left++] = icell;
        if (cells[icell].down)
            down_cells[cur_down++] = icell;
        if (cells[icell].right)
            right_cells[cur_right++] = icell;
        if (cells[icell].up)
            up_cells[cur_up++] = icell;
    }

    //Рассмотрение вертикальных линий
    auto vert_sort_func = [&](int i, int j)
        {
            if (verts[i].y > verts[j].y)
                return true;
            return false;
        };
    for (int iline = 0; iline < problem.get_vert_lines_size(); ++iline)
    {
        line cur_line = problem.get_vert_line(iline);
        if (cur_line.nBox2 == -1)
            continue;

        if (problem.get_box(cur_line.nBox1).is_conductor() || problem.get_box(cur_line.nBox2).is_conductor())
            continue;

        std::vector<int> line_points;
        
        int id = right_cells_sizes[cur_line.nBox1];
        if (right_cells_sizes[cur_line.nBox1] != right_cells_sizes[cur_line.nBox1 + 1])
        {
            if (point_in_line(verts[cells[right_cells[id]].p[2]], problem, cur_line))
                line_points.push_back(cells[right_cells[id]].p[2]);
        }
        for (int i = right_cells_sizes[cur_line.nBox1]; i < right_cells_sizes[cur_line.nBox1 + 1]; ++i)
        {
            if (point_in_line(verts[cells[right_cells[i]].p[1]], problem, cur_line))
                line_points.push_back(cells[right_cells[i]].p[1]);
        }

        id = left_cells_sizes[cur_line.nBox2];
        if (left_cells_sizes[cur_line.nBox2] != left_cells_sizes[cur_line.nBox2 + 1])
        {
            if (point_in_line(verts[cells[left_cells[id]].p[3]], problem, cur_line))
            {
                auto it = std::find(line_points.begin(), line_points.end(), cells[left_cells[id]].p[3]);
                if (it == line_points.end())
                    line_points.push_back(cells[left_cells[id]].p[3]);
            }
        }
        for (int i = left_cells_sizes[cur_line.nBox2]; i < left_cells_sizes[cur_line.nBox2 + 1]; ++i)
        {
            if (point_in_line(verts[cells[left_cells[i]].p[0]], problem, cur_line))
            {
               auto it = std::find(line_points.begin(), line_points.end(), cells[left_cells[i]].p[0]);
               if (it == line_points.end())
                 line_points.push_back(cells[left_cells[i]].p[0]);
            }
        }

        std::sort(line_points.begin(), line_points.end(), vert_sort_func);

        auto it = line_points.begin();
        for (int i = right_cells_sizes[cur_line.nBox1]; i < right_cells_sizes[cur_line.nBox1 + 1]; ++i)
        {
            for (; it != line_points.end(); ++it)
            {
                if (verts[cells[right_cells[i]].p[1]].y < verts[*it].y && verts[cells[right_cells[i]].p[2]].y > verts[*it].y)
                {
                     if(std::find(cells[right_cells[i]].right_points.begin(), cells[right_cells[i]].right_points.end(), *it) == 
                         cells[right_cells[i]].right_points.end())
                         cells[right_cells[i]].right_points.push_back(*it);
                }                  
                if (verts[cells[right_cells[i]].p[1]].y > verts[*it].y)
                    break;
            }
        }

        it = line_points.begin();
        for (int i = left_cells_sizes[cur_line.nBox2]; i < left_cells_sizes[cur_line.nBox2 + 1]; ++i)
        {
            for (; it != line_points.end(); ++it)
            {
                if (verts[cells[left_cells[i]].p[0]].y < verts[*it].y && verts[cells[left_cells[i]].p[3]].y > verts[*it].y)
                {
                    if (std::find(cells[left_cells[i]].left_points.begin(), cells[left_cells[i]].left_points.end(), *it) ==
                        cells[left_cells[i]].left_points.end())
                        cells[left_cells[i]].left_points.push_back(*it);
                }
                if (verts[cells[left_cells[i]].p[1]].y > verts[*it].y)
                    break;
            }
        }
    }

    //Рассмотрение горизонтальных линий
    auto horz_sort_func = [&](int i, int j)
        {
            if (verts[i].x < verts[j].x)
                return true;
            return false;
        };
    for (int iline = 0; iline < problem.get_horz_lines_size(); ++iline)
    {
        line cur_line = problem.get_horz_line(iline);

        if (cur_line.nBox2 == -1)
            continue;

        if (problem.get_box(cur_line.nBox1).is_conductor() || problem.get_box(cur_line.nBox2).is_conductor())
            continue;

        std::vector<int> line_points;

        int id = down_cells_sizes[cur_line.nBox1];
        if (down_cells_sizes[cur_line.nBox1] != down_cells_sizes[cur_line.nBox1 + 1])
        {
            if (point_in_line(verts[cells[down_cells[id]].p[0]], problem, cur_line))
                line_points.push_back(cells[down_cells[id]].p[0]);
        }
        for (int i = down_cells_sizes[cur_line.nBox1]; i < down_cells_sizes[cur_line.nBox1 + 1]; ++i)
        {
            if (point_in_line(verts[cells[down_cells[i]].p[1]], problem, cur_line))
                line_points.push_back(cells[down_cells[i]].p[1]);
        }

        id = up_cells_sizes[cur_line.nBox2];
        if (up_cells_sizes[cur_line.nBox2] != up_cells_sizes[cur_line.nBox2 + 1])
        {
            if (point_in_line(verts[cells[up_cells[id]].p[3]], problem, cur_line))
            {
                auto it = std::find(line_points.begin(), line_points.end(), cells[up_cells[id]].p[3]);
                if (it == line_points.end())
                    line_points.push_back(cells[up_cells[id]].p[3]);
            }
        }
        for (int i = up_cells_sizes[cur_line.nBox2]; i < up_cells_sizes[cur_line.nBox2 + 1]; ++i)
        {
            if (point_in_line(verts[cells[up_cells[i]].p[2]], problem, cur_line))
            {
                auto it = std::find(line_points.begin(), line_points.end(), cells[up_cells[i]].p[2]);
                if (it == line_points.end())
                    line_points.push_back(cells[up_cells[i]].p[2]);
            }
        }

        std::sort(line_points.begin(), line_points.end(), horz_sort_func);

        auto it = line_points.begin();
        for (int i = down_cells_sizes[cur_line.nBox1]; i < down_cells_sizes[cur_line.nBox1 + 1]; ++i)
        {
            for (; it != line_points.end(); ++it)
            {
                if (verts[cells[down_cells[i]].p[0]].x < verts[*it].x && verts[cells[down_cells[i]].p[1]].x > verts[*it].x)
                {
                    if (std::find(cells[down_cells[i]].down_points.begin(), cells[down_cells[i]].down_points.end(), *it) ==
                        cells[down_cells[i]].down_points.end())
                        cells[down_cells[i]].down_points.push_back(*it);
                }
                if (verts[cells[down_cells[i]].p[1]].x < verts[*it].x)
                    break;
            }
        }

        it = line_points.begin();
        for (int i = up_cells_sizes[cur_line.nBox2]; i < up_cells_sizes[cur_line.nBox2 + 1]; ++i)
        {
            for (; it != line_points.end(); ++it)
            {
                if (verts[cells[up_cells[i]].p[3]].x < verts[*it].x && verts[cells[up_cells[i]].p[2]].x > verts[*it].x)
                {
                    if (std::find(cells[up_cells[i]].up_points.begin(), cells[up_cells[i]].up_points.end(), *it) ==
                        cells[up_cells[i]].up_points.end())
                        cells[up_cells[i]].up_points.push_back(*it);
                }
                if (verts[cells[up_cells[i]].p[2]].x < verts[*it].x)
                    break;
            }
        }
    }

    // Простановка индексов для боксов
    std::vector<int> G2L_nets;
    G2L_nets.resize(problem.get_box_size(), -1);
    int net_size = 0;
    for (int ibox = 0; ibox < problem.get_box_size(); ++ibox)
    {
        if (problem.get_box(ibox).is_conductor())
        {
            G2L_nets[ibox] = net_size++;
        }
    }

    double ymin = problem.get_min_y();

    // Обработка вертикальных границ
    for (int ivline = 0; ivline < problem.get_vert_lines_size(); ivline++)
    {
        line cur_line = problem.get_vert_line(ivline);

        if (cur_line.nBox2 == -1)
        {
            bool left_line = false;
            if(point_in_line(problem.get_point(problem.get_box(cur_line.nBox1).p[0]), problem, cur_line))
                left_line = true;
            if (left_line)
            {
                for (int icell = left_cells_sizes[cur_line.nBox1]; icell < left_cells_sizes[cur_line.nBox1 + 1]; ++icell)
                    verts[cells[left_cells[icell]].p[0]].boundary_ind = -1;
                verts[cells[left_cells[left_cells_sizes[cur_line.nBox1]]].p[3]].boundary_ind = -1;
            }
            else
            {
                for (int icell = right_cells_sizes[cur_line.nBox1]; icell < right_cells_sizes[cur_line.nBox1 + 1]; ++icell)
                    verts[cells[right_cells[icell]].p[1]].boundary_ind = -1;
                verts[cells[right_cells[right_cells_sizes[cur_line.nBox1]]].p[2]].boundary_ind = -1;
            }
            continue;
        }

        int condname = -1;
        int dielname = -1;
        if (problem.get_box(cur_line.nBox1).is_conductor() && !problem.get_box(cur_line.nBox2).is_conductor())
        {
            condname = cur_line.nBox1;
            dielname = cur_line.nBox2;
        }

        if (!problem.get_box(cur_line.nBox1).is_conductor() && problem.get_box(cur_line.nBox2).is_conductor())
        {
            condname = cur_line.nBox2;
            dielname = cur_line.nBox1;
        }

        if (condname == -1 && dielname == -1)
            continue;

        bool left_line = false;
        if (point_in_line(problem.get_point(problem.get_box(dielname).p[0]), problem, cur_line))
            left_line = true;
        if (left_line)
        {
            for (int icell = left_cells_sizes[dielname]; icell < left_cells_sizes[dielname + 1]; ++icell)
                verts[cells[left_cells[icell]].p[0]].boundary_ind = G2L_nets[condname];
            verts[cells[left_cells[left_cells_sizes[dielname]]].p[3]].boundary_ind = G2L_nets[condname];
        }
        else
        {
            for (int icell = right_cells_sizes[dielname]; icell < right_cells_sizes[dielname + 1]; ++icell)
                verts[cells[right_cells[icell]].p[1]].boundary_ind = G2L_nets[condname];
            verts[cells[right_cells[right_cells_sizes[dielname]]].p[2]].boundary_ind = G2L_nets[condname];
        }
    }

    // Обработка горизонтальных границ
    for (int ihline = 0; ihline < problem.get_horz_lines_size(); ihline++)
    {
        line cur_line = problem.get_horz_line(ihline);
        int free_ind = -1;
        if (fabs(problem.get_point(cur_line.lp).y - ymin) < 1.e-6)
            free_ind = -2;


        if (cur_line.nBox2 == -1)
        {
            bool down_line = false;
            if (point_in_line(problem.get_point(problem.get_box(cur_line.nBox1).p[0]), problem, cur_line))
                down_line = true;
            if (down_line)
            {
                for (int icell = down_cells_sizes[cur_line.nBox1]; icell < down_cells_sizes[cur_line.nBox1 + 1]; ++icell)
                    verts[cells[down_cells[icell]].p[1]].boundary_ind = free_ind;
                verts[cells[down_cells[down_cells_sizes[cur_line.nBox1]]].p[0]].boundary_ind = free_ind;
            }
            else
            {
                for (int icell = up_cells_sizes[cur_line.nBox1]; icell < up_cells_sizes[cur_line.nBox1 + 1]; ++icell)
                    verts[cells[up_cells[icell]].p[2]].boundary_ind = free_ind;
                verts[cells[up_cells[up_cells_sizes[cur_line.nBox1]]].p[3]].boundary_ind = free_ind;
            }
            continue;
        }

        int condname = -1;
        int dielname = -1;
        if (problem.get_box(cur_line.nBox1).is_conductor() && !problem.get_box(cur_line.nBox2).is_conductor())
        {
            condname = cur_line.nBox1;
            dielname = cur_line.nBox2;
        }

        if (!problem.get_box(cur_line.nBox1).is_conductor() && problem.get_box(cur_line.nBox2).is_conductor())
        {
            condname = cur_line.nBox2;
            dielname = cur_line.nBox1;
        }

        if (condname == -1 && dielname == -1)
            continue;

        bool down_line = false;
        if (point_in_line(problem.get_point(problem.get_box(dielname).p[0]), problem, cur_line))
            down_line = true;
        if (down_line)
        {
            for (int icell = down_cells_sizes[dielname]; icell < down_cells_sizes[dielname + 1]; ++icell)
                verts[cells[down_cells[icell]].p[1]].boundary_ind = G2L_nets[condname];
            verts[cells[down_cells[down_cells_sizes[dielname]]].p[0]].boundary_ind = G2L_nets[condname];
        }
        else
        {
            for (int icell = up_cells_sizes[dielname]; icell < up_cells_sizes[dielname + 1]; ++icell)
                verts[cells[up_cells[icell]].p[2]].boundary_ind = G2L_nets[condname];
            verts[cells[up_cells[up_cells_sizes[dielname]]].p[3]].boundary_ind = G2L_nets[condname];
        }
    }


    delete[] left_cells_sizes;
    delete[] down_cells_sizes;
    delete[] right_cells_sizes;
    delete[] up_cells_sizes;
    delete[] left_cells;
    delete[] down_cells;
    delete[] right_cells;
    delete[] up_cells;
    return 0;
}

int FEM_solver::classify_verts_ver2(task& problem)
{
    int* left_cells_sizes = new int[number_of_groups + 1];
    int* down_cells_sizes = new int[number_of_groups + 1];
    int* right_cells_sizes = new int[number_of_groups + 1];
    int* up_cells_sizes = new int[number_of_groups + 1];

    for (int i = 0; i < number_of_groups + 1; ++i)
    {
        left_cells_sizes[i] = 0;
        down_cells_sizes[i] = 0;
        right_cells_sizes[i] = 0;
        up_cells_sizes[i] = 0;
    }

    for (int icell = 0; icell < cells.size(); ++icell)
    {
        if (cells[icell].left)
            left_cells_sizes[cells[icell].group + 1]++;
        if (cells[icell].down)
            down_cells_sizes[cells[icell].group + 1]++;
        if (cells[icell].right)
            right_cells_sizes[cells[icell].group + 1]++;
        if (cells[icell].up)
            up_cells_sizes[cells[icell].group + 1]++;
    }

    for (int i = 0; i < number_of_groups; ++i)
    {
        left_cells_sizes[i + 1] += left_cells_sizes[i];
        down_cells_sizes[i + 1] += down_cells_sizes[i];
        right_cells_sizes[i + 1] += right_cells_sizes[i];
        up_cells_sizes[i + 1] += up_cells_sizes[i];
    }

    int cur_left = 0;
    int cur_down = 0;
    int cur_right = 0;
    int cur_up = 0;
    int* left_cells = new int[left_cells_sizes[number_of_groups]];
    int* down_cells = new int[down_cells_sizes[number_of_groups]];
    int* right_cells = new int[right_cells_sizes[number_of_groups]];
    int* up_cells = new int[up_cells_sizes[number_of_groups]];

    for (int icell = 0; icell < cells.size(); ++icell)
    {
        if (cells[icell].left)
            left_cells[cur_left++] = icell;
        if (cells[icell].down)
            down_cells[cur_down++] = icell;
        if (cells[icell].right)
            right_cells[cur_right++] = icell;
        if (cells[icell].up)
            up_cells[cur_up++] = icell;
    }

    //Рассмотрение вертикальных линий
    bool master_nbox1; // Флаг того что мастер-бокс - nbox1
    auto vert_sort_func = [&](int i, int j)
        {
            if (verts[i].y > verts[j].y)
                return true;
            return false;
        };
    for (int iline = 0; iline < problem.get_vert_lines_size(); ++iline)
    {
        line cur_line = problem.get_vert_line(iline);
        if (cur_line.nBox2 == -1)
            continue;

        if (problem.get_box(cur_line.nBox1).is_conductor() || problem.get_box(cur_line.nBox2).is_conductor())
            continue;

        master_nbox1 = (problem.get_box(cur_line.nBox1).a > problem.get_box(cur_line.nBox2).a) ? true : false;

        std::vector<int> line_points;

        int id = right_cells_sizes[cur_line.nBox1];
        if (right_cells_sizes[cur_line.nBox1] != right_cells_sizes[cur_line.nBox1 + 1])
        {
            if (point_in_line(verts[cells[right_cells[id]].p[2]], problem, cur_line))
                line_points.push_back(cells[right_cells[id]].p[2]);
        }
        for (int i = right_cells_sizes[cur_line.nBox1]; i < right_cells_sizes[cur_line.nBox1 + 1]; ++i)
        {
            if (point_in_line(verts[cells[right_cells[i]].p[1]], problem, cur_line))
                line_points.push_back(cells[right_cells[i]].p[1]);
        }

        id = left_cells_sizes[cur_line.nBox2];
        if (left_cells_sizes[cur_line.nBox2] != left_cells_sizes[cur_line.nBox2 + 1])
        {
            if (point_in_line(verts[cells[left_cells[id]].p[3]], problem, cur_line))
            {
                auto it = std::find(line_points.begin(), line_points.end(), cells[left_cells[id]].p[3]);
                if (it == line_points.end())
                    line_points.push_back(cells[left_cells[id]].p[3]);
            }
        }
        for (int i = left_cells_sizes[cur_line.nBox2]; i < left_cells_sizes[cur_line.nBox2 + 1]; ++i)
        {
            if (point_in_line(verts[cells[left_cells[i]].p[0]], problem, cur_line))
            {
                auto it = std::find(line_points.begin(), line_points.end(), cells[left_cells[i]].p[0]);
                if (it == line_points.end())
                    line_points.push_back(cells[left_cells[i]].p[0]);
            }
        }

        std::sort(line_points.begin(), line_points.end(), vert_sort_func);

        for (int ivert = 1; ivert < line_points.size() - 1; ivert++)
        {
            if (verts[line_points[ivert]].left && !verts[line_points[ivert]].right && master_nbox1 && !verts[line_points[ivert]].is_corner())
                verts[line_points[ivert]].approx = true;
            if (!verts[line_points[ivert]].left && verts[line_points[ivert]].right && !master_nbox1 && !verts[line_points[ivert]].is_corner())
                verts[line_points[ivert]].approx = true;
        }

        
        for (int ivert = 1; ivert < line_points.size() - 1; ivert++)
        {
            if (verts[line_points[ivert]].approx)
            {
                int cur_vert = ivert - 1;
                while (verts[line_points[cur_vert]].approx)
                    cur_vert--;
                verts[line_points[ivert]].vert1 = line_points[cur_vert];
                cur_vert = ivert + 1;
                while (verts[line_points[cur_vert]].approx)
                    cur_vert++;
                verts[line_points[ivert]].vert2 = line_points[cur_vert];
            }
        }

        auto it = line_points.begin();
        for (int i = right_cells_sizes[cur_line.nBox1]; i < right_cells_sizes[cur_line.nBox1 + 1]; ++i)
        {
            for (; it != line_points.end(); ++it)
            {
                if (verts[cells[right_cells[i]].p[1]].y < verts[*it].y && verts[cells[right_cells[i]].p[2]].y > verts[*it].y)
                {
                    if (std::find(cells[right_cells[i]].right_points.begin(), cells[right_cells[i]].right_points.end(), *it) ==
                        cells[right_cells[i]].right_points.end() && verts[*it].is_corner())
                        cells[right_cells[i]].right_points.push_back(*it);
                }
                if (verts[cells[right_cells[i]].p[1]].y > verts[*it].y)
                    break;
            }
        }

        it = line_points.begin();
        for (int i = left_cells_sizes[cur_line.nBox2]; i < left_cells_sizes[cur_line.nBox2 + 1]; ++i)
        {
            for (; it != line_points.end(); ++it)
            {
                if (verts[cells[left_cells[i]].p[0]].y < verts[*it].y && verts[cells[left_cells[i]].p[3]].y > verts[*it].y)
                {
                    if (std::find(cells[left_cells[i]].left_points.begin(), cells[left_cells[i]].left_points.end(), *it) ==
                        cells[left_cells[i]].left_points.end() && verts[*it].is_corner())
                        cells[left_cells[i]].left_points.push_back(*it);
                }
                if (verts[cells[left_cells[i]].p[1]].y > verts[*it].y)
                    break;
            }
        }
    }

    //Рассмотрение горизонтальных линий
    auto horz_sort_func = [&](int i, int j)
        {
            if (verts[i].x < verts[j].x)
                return true;
            return false;
        };
    for (int iline = 0; iline < problem.get_horz_lines_size(); ++iline)
    {
        line cur_line = problem.get_horz_line(iline);
        if (cur_line.nBox2 == -1)
            continue;

        if (problem.get_box(cur_line.nBox1).is_conductor() || problem.get_box(cur_line.nBox2).is_conductor())
            continue;

        std::vector<int> line_points;

        master_nbox1 = (problem.get_box(cur_line.nBox1).a > problem.get_box(cur_line.nBox2).a) ? true : false;

        int id = down_cells_sizes[cur_line.nBox1];
        if (down_cells_sizes[cur_line.nBox1] != down_cells_sizes[cur_line.nBox1 + 1])
        {
            if (point_in_line(verts[cells[down_cells[id]].p[0]], problem, cur_line))
                line_points.push_back(cells[down_cells[id]].p[0]);
        }
        for (int i = down_cells_sizes[cur_line.nBox1]; i < down_cells_sizes[cur_line.nBox1 + 1]; ++i)
        {
            if (point_in_line(verts[cells[down_cells[i]].p[1]], problem, cur_line))
                line_points.push_back(cells[down_cells[i]].p[1]);
        }

        id = up_cells_sizes[cur_line.nBox2];
        if (up_cells_sizes[cur_line.nBox2] != up_cells_sizes[cur_line.nBox2 + 1])
        {
            if (point_in_line(verts[cells[up_cells[id]].p[3]], problem, cur_line))
            {
                auto it = std::find(line_points.begin(), line_points.end(), cells[up_cells[id]].p[3]);
                if (it == line_points.end())
                    line_points.push_back(cells[up_cells[id]].p[3]);
            }
        }
        for (int i = up_cells_sizes[cur_line.nBox2]; i < up_cells_sizes[cur_line.nBox2 + 1]; ++i)
        {
            if (point_in_line(verts[cells[up_cells[i]].p[2]], problem, cur_line))
            {
                auto it = std::find(line_points.begin(), line_points.end(), cells[up_cells[i]].p[2]);
                if (it == line_points.end())
                    line_points.push_back(cells[up_cells[i]].p[2]);
            }
        }

        std::sort(line_points.begin(), line_points.end(), horz_sort_func);

        for (int ivert = 1; ivert < line_points.size() - 1; ivert++)
        {
            if (verts[line_points[ivert]].up && !verts[line_points[ivert]].down && master_nbox1 && !verts[line_points[ivert]].is_corner())
                verts[line_points[ivert]].approx = true;
            if (!verts[line_points[ivert]].up && verts[line_points[ivert]].down && !master_nbox1 && !verts[line_points[ivert]].is_corner())
                verts[line_points[ivert]].approx = true;
        }


        for (int ivert = 1; ivert < line_points.size() - 1; ivert++)
        {
            if (verts[line_points[ivert]].approx)
            {
                int cur_vert = ivert - 1;
                while (verts[line_points[cur_vert]].approx)
                    cur_vert--;
                verts[line_points[ivert]].vert1 = line_points[cur_vert];
                cur_vert = ivert + 1;
                while (verts[line_points[cur_vert]].approx)
                    cur_vert++;
                verts[line_points[ivert]].vert2 = line_points[cur_vert];
            }
        }


        auto it = line_points.begin();
        for (int i = down_cells_sizes[cur_line.nBox1]; i < down_cells_sizes[cur_line.nBox1 + 1]; ++i)
        {
            for (; it != line_points.end(); ++it)
            {
                if (verts[cells[down_cells[i]].p[0]].x < verts[*it].x && verts[cells[down_cells[i]].p[1]].x > verts[*it].x)
                {
                    if (std::find(cells[down_cells[i]].down_points.begin(), cells[down_cells[i]].down_points.end(), *it) ==
                        cells[down_cells[i]].down_points.end() && verts[*it].is_corner())
                        cells[down_cells[i]].down_points.push_back(*it);
                }
                if (verts[cells[down_cells[i]].p[1]].x < verts[*it].x)
                    break;
            }
        }

        it = line_points.begin();
        for (int i = up_cells_sizes[cur_line.nBox2]; i < up_cells_sizes[cur_line.nBox2 + 1]; ++i)
        {
            for (; it != line_points.end(); ++it)
            {
                if (verts[cells[up_cells[i]].p[3]].x < verts[*it].x && verts[cells[up_cells[i]].p[2]].x > verts[*it].x)
                {
                    if (std::find(cells[up_cells[i]].up_points.begin(), cells[up_cells[i]].up_points.end(), *it) ==
                        cells[up_cells[i]].up_points.end() && verts[*it].is_corner())
                        cells[up_cells[i]].up_points.push_back(*it);
                }
                if (verts[cells[up_cells[i]].p[2]].x < verts[*it].x)
                    break;
            }
        }
    }

    // Простановка индексов для боксов
    std::vector<int> G2L_nets;
    G2L_nets.resize(problem.get_box_size(), -1);
    int net_size = 0;
    for (int ibox = 0; ibox < problem.get_box_size(); ++ibox)
    {
        if (problem.get_box(ibox).is_conductor())
            G2L_nets[ibox] = net_size++;
    }

    double ymin = problem.get_min_y();

    // Обработка вертикальных границ
    for (int ivline = 0; ivline < problem.get_vert_lines_size(); ivline++)
    {
        line cur_line = problem.get_vert_line(ivline);

        if (cur_line.nBox2 == -1)
        {
            bool left_line = false;
            if (point_in_line(problem.get_point(problem.get_box(cur_line.nBox1).p[0]), problem, cur_line))
                left_line = true;
            if (left_line)
            {
                for (int icell = left_cells_sizes[cur_line.nBox1]; icell < left_cells_sizes[cur_line.nBox1 + 1]; ++icell)
                    verts[cells[left_cells[icell]].p[0]].boundary_ind = -1;
                verts[cells[left_cells[left_cells_sizes[cur_line.nBox1]]].p[3]].boundary_ind = -1;
            }
            else
            {
                for (int icell = right_cells_sizes[cur_line.nBox1]; icell < right_cells_sizes[cur_line.nBox1 + 1]; ++icell)
                    verts[cells[right_cells[icell]].p[1]].boundary_ind = -1;
                verts[cells[right_cells[right_cells_sizes[cur_line.nBox1]]].p[2]].boundary_ind = -1;
            }
            continue;
        }

        int condname = -1;
        int dielname = -1;
        if (problem.get_box(cur_line.nBox1).is_conductor() && !problem.get_box(cur_line.nBox2).is_conductor())
        {
            condname = cur_line.nBox1;
            dielname = cur_line.nBox2;
        }

        if (!problem.get_box(cur_line.nBox1).is_conductor() && problem.get_box(cur_line.nBox2).is_conductor())
        {
            condname = cur_line.nBox2;
            dielname = cur_line.nBox1;
        }

        if (condname == -1 && dielname == -1)
            continue;

        if (verts[34].boundary_ind == 0)
            std::cout << 1;

        bool left_line = false;
        if (point_in_line(problem.get_point(problem.get_box(dielname).p[0]), problem, cur_line))
            left_line = true;
        if (left_line)
        {
            for (int icell = left_cells_sizes[dielname]; icell < left_cells_sizes[dielname + 1]; ++icell)
                verts[cells[left_cells[icell]].p[0]].boundary_ind = G2L_nets[condname];
            verts[cells[left_cells[left_cells_sizes[dielname]]].p[3]].boundary_ind = G2L_nets[condname];
        }
        else
        {
            for (int icell = right_cells_sizes[dielname]; icell < right_cells_sizes[dielname + 1]; ++icell)
                verts[cells[right_cells[icell]].p[1]].boundary_ind = G2L_nets[condname];
                
            verts[cells[right_cells[right_cells_sizes[dielname]]].p[2]].boundary_ind = G2L_nets[condname];
        }
    }

    // Обработка горизонтальных границ
    for (int ihline = 0; ihline < problem.get_horz_lines_size(); ihline++)
    {
        line cur_line = problem.get_horz_line(ihline);
        int free_ind = -1;
        if (fabs(problem.get_point(cur_line.lp).y - ymin) < 1.e-6)
            free_ind = -2;


        if (cur_line.nBox2 == -1)
        {
            bool down_line = false;
            if (point_in_line(problem.get_point(problem.get_box(cur_line.nBox1).p[0]), problem, cur_line))
                down_line = true;
            if (down_line)
            {
                for (int icell = down_cells_sizes[cur_line.nBox1]; icell < down_cells_sizes[cur_line.nBox1 + 1]; ++icell)
                    verts[cells[down_cells[icell]].p[1]].boundary_ind = free_ind;
                verts[cells[down_cells[down_cells_sizes[cur_line.nBox1]]].p[0]].boundary_ind = free_ind;
            }
            else
            {
                for (int icell = up_cells_sizes[cur_line.nBox1]; icell < up_cells_sizes[cur_line.nBox1 + 1]; ++icell)
                    verts[cells[up_cells[icell]].p[2]].boundary_ind = free_ind;
                verts[cells[up_cells[up_cells_sizes[cur_line.nBox1]]].p[3]].boundary_ind = free_ind;
            }
            continue;

        }

        int condname = -1;
        int dielname = -1;
        if (problem.get_box(cur_line.nBox1).is_conductor() && !problem.get_box(cur_line.nBox2).is_conductor())
        {
            condname = cur_line.nBox1;
            dielname = cur_line.nBox2;
        }

        if (!problem.get_box(cur_line.nBox1).is_conductor() && problem.get_box(cur_line.nBox2).is_conductor())
        {
            condname = cur_line.nBox2;
            dielname = cur_line.nBox1;
        }

        if (condname == -1 && dielname == -1)
            continue;

        bool down_line = false;
        if (point_in_line(problem.get_point(problem.get_box(dielname).p[0]), problem, cur_line))
            down_line = true;
        if (down_line)
        {
            for (int icell = down_cells_sizes[dielname]; icell < down_cells_sizes[dielname + 1]; ++icell)
                verts[cells[down_cells[icell]].p[1]].boundary_ind = G2L_nets[condname];
            verts[cells[down_cells[down_cells_sizes[dielname]]].p[0]].boundary_ind = G2L_nets[condname];
        }
        else
        {
            for (int icell = up_cells_sizes[dielname]; icell < up_cells_sizes[dielname + 1]; ++icell)
                verts[cells[up_cells[icell]].p[2]].boundary_ind = G2L_nets[condname];
            verts[cells[up_cells[up_cells_sizes[dielname]]].p[3]].boundary_ind = G2L_nets[condname];
        }
    }


    delete[] left_cells_sizes;
    delete[] down_cells_sizes;
    delete[] right_cells_sizes;
    delete[] up_cells_sizes;
    delete[] left_cells;
    delete[] down_cells;
    delete[] right_cells;
    delete[] up_cells;
    return 0;
}

void FEM_solver::debug_mesh_with_classification()
{
    std::ofstream file("C:\\Users\\egora\\Desktop\\Aspa\\Codes\\General_progs\\debug_mesh.txt");
    file << cells.size() << '\n';
    for (FEM_cell cell : cells)
    {
        int vert_number = 0;
        for (int i = 0; i < 4; ++i)
        {
            if (!verts[cell.p[i]].approx)
                vert_number++;
        }
        vert_number += cell.down_points.size();
        vert_number += cell.right_points.size();
        vert_number += cell.up_points.size();
        vert_number += cell.left_points.size();
        file << cell.p[0] << ' ' << cell.p[1] << " " << cell.p[2] << " " << cell.p[3] << " " << vert_number << '\n';
    }

    file << verts.size() << '\n';
    for (FEM_point point : verts)
    {
        file << point.x << ' ' << point.y << '\n';
    }
    for (FEM_point point : verts)
    {
        file << ((point.approx)? 0 : 1) << '\n';
    }
    

    file.close();
}

int FEM_solver::perform_new_numeration()
{
    number_of_independent = 0;
    for (int ivert = 0; ivert < verts.size(); ++ivert)
    {
        if (!verts[ivert].approx)
            verts[ivert].local_num = number_of_independent++;
    }
    return 0;
}

int FEM_solver::make_system_of_linear_equations()
{
    std::vector<int> local_verts;

    int cur_group = -1;
    std::vector<std::vector<double>> local_matrix;
    std::vector<std::vector<double>> group_matrix;
    bool simple_matrix;

    for (int icell = 0; icell < cells.size(); ++icell)
    {
        simple_matrix = find_local_verts(local_verts, icell);
        if (cur_group == cells[icell].group && simple_matrix)
        {
            local_matrix = group_matrix;
        }
        else if(simple_matrix)
        {
            cur_group = cells[icell].group;
            calculate_local_matrix_ver2(local_verts, icell, local_matrix);
            group_matrix = local_matrix;
        }
        else
        {
            calculate_local_matrix_ver2(local_verts, icell, local_matrix);
        }

        for (int i = 0; i < local_verts.size(); ++i)
        {
            if (verts[local_verts[i]].boundary_ind == -2 || verts[local_verts[i]].boundary_ind >= 0)
                continue;
            for (int j = 0; j < local_verts.size(); ++j)
            {
                if (verts[local_verts[j]].boundary_ind == -2)
                    continue;
                LS_Triplets.push_back({ local_verts[i], local_verts[j],
                    local_matrix[i][j] });
            } 
        }
    }

    right = Eigen::VectorXd(number_of_independent);
    
    return 0;
}

int FEM_solver::make_right(int net1, int net2)
{
    right.setZero();
    for (int ivert = 0; ivert < verts.size(); ++ivert)
    {
        if (verts[ivert].boundary_ind == net1 ||
            verts[ivert].boundary_ind == net2)
        {
            right.coeffRef(verts[ivert].local_num) = 1.0;
        }
    }
    return 0;
}

int FEM_solver::solve_system()
{
    std::vector<Eigen::Triplet<double>> triplets;
    for (int itrip = 0; itrip < LS_Triplets.size(); itrip++)
    {
        if (verts[LS_Triplets[itrip].col()].boundary_ind >= 0)
        {
            right.coeffRef(verts[LS_Triplets[itrip].row()].local_num) -= LS_Triplets[itrip].value() * right.coeffRef(verts[LS_Triplets[itrip].col()].local_num);
            continue;
        }
        triplets.push_back({ verts[LS_Triplets[itrip].row()].local_num, verts[LS_Triplets[itrip].col()].local_num, LS_Triplets[itrip].value() });
    }

    for (int ivert = 0; ivert < verts.size(); ++ivert)
    {
        if (verts[ivert].boundary_ind >= 0 || verts[ivert].boundary_ind == -2)
            triplets.push_back({ verts[ivert].local_num, verts[ivert].local_num, 1.0 });
    }

    std::ofstream file("triplets.txt");
    for (auto tr : triplets)
        file << tr.row() << " " << tr.col() << " " << tr.value() << '\n';
    file.close();

    Eigen::SparseMatrix<double> SMatrix(number_of_independent, number_of_independent);
    
    SMatrix.setFromTriplets(triplets.begin(), triplets.end());
    //Eigen::SparseLU<Eigen::SparseMatrix<double>> MSolver;
    //MSolver.analyzePattern(SMatrix);
    //MSolver.factorize(SMatrix);
    //if (MSolver.info() != Eigen::Success)
    //{
    //    std::cerr << "Solver's error." << std::endl;
    //}

    //SMatrix.makeCompressed();
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, 1, Eigen::DiagonalPreconditioner<double>>  MSolver(SMatrix);
    
    Eigen::VectorXd ans(number_of_independent);
    ans = MSolver.solve(right);
    if (MSolver.info() != Eigen::Success)
    {
        std::cerr << "Solver's error." << std::endl;
    }

    answer.clear();
    answer.resize(verts.size());
    for (int ivert = 0; ivert < verts.size(); ++ivert)
    {
        if (verts[ivert].approx)
        {
            double coeff1;
            double coeff2;
            if (fabs(verts[verts[ivert].vert1].x - verts[verts[ivert].vert2].x) < 1.e-6)
            {
                get_coeffs_Laplace_interp(verts[verts[ivert].vert1].y, verts[verts[ivert].vert2].y, verts[ivert].y, coeff1, coeff2);
            }
            else
            {
                get_coeffs_Laplace_interp(verts[verts[ivert].vert1].x, verts[verts[ivert].vert2].x, verts[ivert].x, coeff1, coeff2);
            }
            answer[ivert] = coeff1 * ans.coeffRef(verts[verts[ivert].vert1].local_num) + coeff2 * ans.coeffRef(verts[verts[ivert].vert2].local_num);
        }
        else
        {
            answer[ivert] = ans.coeffRef(verts[ivert].local_num);
        }
    }

    //debug_save_portrait(triplets);
    //debug_save_vtk("answer");
    return 0;
}

int FEM_solver::debug_save_vtk(std::string filename)
{
    std::ofstream file;
    file.open(filename + ".pvd");
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type = \"Collection\" version = \"0.1\">\n";
    file << "\t<Collection>\n";
    file << "<DataSet timestep=\"0\" part=\"0\" file=\"" << filename + ".vtu" << "\" />\n";
    file << "\t</Collection>\n";
    file << "</VTKFile>";
    file.close();

    file.open(filename + ".vtu");
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\"  version=\"0.1\"  >\n";
    file << "<UnstructuredGrid>\n";
    file << "<Piece  NumberOfPoints=\"" << verts.size()
        << "\" NumberOfCells=\"" << cells.size() << "\">\n";
    file << "<Points>\n";
    file << "<DataArray  type = \"Float64\"  NumberOfComponents = \"3\"  format = \"ascii\">";
    for (int i = 0; i < verts.size(); ++i)
        file << verts[i].x << ' ' << verts[i].y << " 0  ";
    file << "</DataArray>\n";
    file << "</Points>\n";
    file << "<Cells>\n";
    file << "<DataArray  type = \"UInt32\"  Name = \"connectivity\"  format = \"ascii\">";
    for (FEM_cell el : cells)
        file << el.p[3] << ' ' << el.p[0] << " " << el.p[1] << "  " << el.p[2] << "  ";
    file << "</DataArray>\n";
    file << "<DataArray  type = \"UInt32\"  Name = \"offsets\"  format = \"ascii\">";
    for (int i = 4; i <= cells.size() * 4; i += 4)
        file << i << ' ';
    file << "</DataArray>\n";
    file << "<DataArray  type = \"UInt32\"  Name = \"types\"  format = \"ascii\">";
    for (int i = 0; i <= cells.size(); ++i)
        file << 9 << ' ';
    file << "</DataArray>\n";
    file << "</Cells>";
    file << "<PointData  Scalars = \"U\">\n";
    file << "<DataArray  type = \"Float64\"  Name = \"U\"  format = \"ascii\">";
    for (int i = 0; i < verts.size(); i++)
        file << answer[i] << ' ';
    file << "</DataArray>\n";
    file << "</PointData>\n";
    file << "</Piece>\n";
    file << "</UnstructuredGrid>\n";
    file << "</VTKFile>";
    file.close();

    return 0;
}

void FEM_solver::debug_approx()
{
    std::ofstream file("test_verts.txt");
    for (int ivert = 0; ivert < verts.size(); ++ivert)
    {
        if (verts[ivert].approx)
        {
            int p1 = std::max(verts[ivert].vert1, verts[ivert].vert2);
            int p2 = std::min(verts[ivert].vert1, verts[ivert].vert2);
            file << ivert <<" " << p1 << ' ' << p2 << '\n';
        }
        else
        {
            file << ivert << " " << -1 << ' ' << -1 << '\n';
        }
        //file << ivert << " " << verts[ivert].x << " " << verts[ivert].y << '\n';
    }
    /*for (int icell = 0; icell < cells.size(); ++icell)
        file << icell << " " << cells[icell].p[0] << " " << cells[icell].p[1] << " " << cells[icell].p[2] << " " << cells[icell].p[3] << '\n';*/
    file.close();
}

double FEM_solver::get_C()
{
    double Cij = 0.0;
    std::vector<double> qxval = { -1.0 / sqrt(3.0), -1.0 / sqrt(3.0),  1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
    std::vector<double> qyval = { -1.0 / sqrt(3.0), 1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 1.0 / sqrt(3.0) };
    std::vector<double> cur_qxvals;
    std::vector<double> cur_qyvals;
    cur_qxvals.resize(4);
    cur_qyvals.resize(4);

    for (int icell = 0; icell < cells.size(); ++icell)
    {
        double Er = cells[icell].eps;

        int Nx = 2;
        int Ny = 2;

        bool ref_vert_left;
        bool ref_horz_up;

        if (cells[icell].left_points.size() != 0)
        {
            ref_vert_left = true;
            Ny += cells[icell].left_points.size();
        }
        else
        {
            ref_vert_left = false;
            Ny += cells[icell].right_points.size();
        }

        if (cells[icell].up_points.size() != 0)
        {
            ref_horz_up = true;
            Nx += cells[icell].up_points.size();
        }
        else
        {
            ref_horz_up = false;
            Nx += cells[icell].down_points.size();
        }

        std::vector<std::vector<int>> geom_elem_matrix;
        geom_elem_matrix.resize(Ny);
        for (int i = 0; i < Ny; ++i)
        {
            geom_elem_matrix[i].resize(Nx, -1);
        }

        geom_elem_matrix[0][0] = cells[icell].p[3];
        for (int i = 1; i < cells[icell].up_points.size() + 1; ++i)
            geom_elem_matrix[0][i] = cells[icell].up_points[i - 1];
        geom_elem_matrix[0][Nx - 1] = cells[icell].p[2];
        for (int i = 1; i < cells[icell].right_points.size() + 1; ++i)
            geom_elem_matrix[Ny - 1 - i][Nx - 1] = cells[icell].right_points[Ny - 2 - i];
        geom_elem_matrix[Ny - 1][Nx - 1] = cells[icell].p[1];
        for (int i = 1; i < cells[icell].down_points.size() + 1; ++i)
            geom_elem_matrix[Ny - 1][i] = cells[icell].down_points[i - 1];
        geom_elem_matrix[Ny - 1][0] = cells[icell].p[0];
        for (int i = 1; i < cells[icell].left_points.size() + 1; ++i)
            geom_elem_matrix[Ny - 1 - i][0] = cells[icell].left_points[Ny - 2 - i];

        std::vector<std::vector<double>> coeffs;
        coeffs.resize(Ny);
        for (int i = 0; i < Ny; ++i)
            coeffs[i].resize(Nx, 0.0);

        for (int i = 0; i < Ny; i++)
        {
            for (int j = 0; j < Nx; j++)
            {
                if (geom_elem_matrix[i][j] == -1)
                {
                    std::vector<int> dof1;
                    std::vector<double> coeffs1;
                    double vcoeff1, vcoeff2;
                    double hcoeff1, hcoeff2;
                    dof1.resize(6);
                    coeffs1.resize(6);

                    double x1;
                    double y1;

                    if (ref_vert_left)
                        y1 = verts[geom_elem_matrix[i][0]].y;
                    else
                        y1 = verts[geom_elem_matrix[i][Nx - 1]].y;

                    if (ref_horz_up)
                        x1 = verts[geom_elem_matrix[0][j]].x;
                    else
                        x1 = verts[geom_elem_matrix[Ny - 1][j]].x;

                    if (ref_horz_up)
                    {
                        dof1[2] = geom_elem_matrix[0][j];
                        dof1[0] = geom_elem_matrix[Ny - 1][0];
                        dof1[1] = geom_elem_matrix[Ny - 1][Nx - 1];
                    }
                    else
                    {
                        dof1[2] = geom_elem_matrix[Ny - 1][j];
                        dof1[0] = geom_elem_matrix[0][0];
                        dof1[1] = geom_elem_matrix[0][Nx - 1];
                    }

                    get_coeffs_Laplace_interp(verts[dof1[2]].y, verts[dof1[0]].y, y1, vcoeff1, vcoeff2);
                    get_coeffs_Laplace_interp(verts[dof1[0]].x, verts[dof1[1]].x, x1, hcoeff1, hcoeff2);
                    coeffs1[2] = 0.5 * vcoeff1;
                    coeffs1[0] = 0.5 * vcoeff2 * hcoeff1;
                    coeffs1[1] = 0.5 * vcoeff2 * hcoeff2;

                    if (ref_vert_left)
                    {
                        dof1[5] = geom_elem_matrix[i][0];
                        dof1[3] = geom_elem_matrix[0][Nx - 1];
                        dof1[4] = geom_elem_matrix[Ny - 1][Nx - 1];
                    }
                    else
                    {
                        dof1[5] = geom_elem_matrix[i][Nx - 1];
                        dof1[3] = geom_elem_matrix[0][0];
                        dof1[4] = geom_elem_matrix[Ny - 1][0];
                    }


                    get_coeffs_Laplace_interp(verts[dof1[3]].y, verts[dof1[4]].y, y1, vcoeff1, vcoeff2);
                    get_coeffs_Laplace_interp(verts[dof1[5]].x, verts[dof1[3]].x, x1, hcoeff1, hcoeff2);
                    coeffs1[5] = 0.5 * hcoeff1;
                    coeffs1[3] = 0.5 * hcoeff2 * vcoeff1;
                    coeffs1[4] = 0.5 * hcoeff2 * vcoeff2;

                    for (int k = 0; k < 6; k++)
                        coeffs[i][j] += answer[dof1[k]] * coeffs1[k];
                }
                else
                {
                    coeffs[i][j] = answer[geom_elem_matrix[i][j]];
                }
            }
        }

        // Рассмотр внутренних квадратов
        for (int i = 0; i < Ny - 1; ++i)
        {
            int vert_up = std::max(geom_elem_matrix[i][0], geom_elem_matrix[i][Nx - 1]);
            int vert_down = std::max(geom_elem_matrix[i + 1][0], geom_elem_matrix[i + 1][Nx - 1]);;
            double b = verts[vert_up].y - verts[vert_down].y;
            double yc = (verts[vert_up].y + verts[vert_down].y) * 0.5;
            b *= 0.5;
            for (int j = 0; j < Nx - 1; ++j)
            {
                int vert_left = std::max(geom_elem_matrix[0][j], geom_elem_matrix[Ny - 1][j]);
                int vert_right = std::max(geom_elem_matrix[0][j + 1], geom_elem_matrix[Ny - 1][j + 1]);;
                double a = verts[vert_right].x - verts[vert_left].x;
                a *= 0.5;
                double xc = (verts[vert_right].x + verts[vert_left].x) * 0.5;
                double S = 4 * a * b;
                for (int qind = 0; qind < 4; ++qind)
                {
                    double dx[4];
                    double dy[4];
                    for (int p = 0; p < 4; ++p)
                    {
                        dx[p] = get_dNdx(p, a, b, a * qxval[qind], b * qyval[qind]);
                        dy[p] = get_dNdy(p, a, b, a * qxval[qind], b * qyval[qind]);
                    }

                    dx[3] *= coeffs[i][j];
                    dy[3] *= coeffs[i][j];
                    dx[0] *= coeffs[i+1][j];
                    dy[0] *= coeffs[i+1][j];
                    dx[1] *= coeffs[i+1][j+1];
                    dy[1] *= coeffs[i+1][j+1];
                    dx[2] *= coeffs[i][j + 1];
                    dy[2] *= coeffs[i][j + 1];

                    double sumdx = 0.0;
                    double sumdy = 0.0;
                    for (int p = 0; p < 4; ++p)
                    {
                        sumdx += dx[p];
                        sumdy += dy[p];
                    }

                    Cij += Er * a * b * (sumdx * sumdx + sumdy * sumdy);
                }
            }

        }
    }

    Cij *= 8.85e-3;
    return Cij;
}

void FEM_solver::debug_save_portrait(std::vector<Eigen::Triplet<double>>& triplets)
{
    std::ofstream file("portrait.txt");
    for (auto tr : triplets)
        file << tr.row() << " " << tr.col() << '\n';
    file.close();
}

FEM_point::FEM_point(double x, double y): point(x, y)
{
    left = false;
    down = false;
    right = false;
    up = false;
    approx = false;
    boundary_ind = -3;
    local_num = -1;
}

bool FEM_point::operator==(FEM_point a)
{
    return this->point::operator==(a) ;
}

bool FEM_point::is_corner()
{
    if (left && up) return true;
    if (left && down) return true;
    if (right && up) return true;
    if (right && down) return true;
    return false;
}

FEM_cell::FEM_cell()
{
    left = false;
    down = false;
    right = false;
    up = false;
}

bool point_in_line(point p, task& problem, line curline)
{
    if (fabs(p.x - problem.get_point(curline.lp).x) < 1.e-8 &&
        fabs(p.y - problem.get_point(curline.lp).y) < 1.e-8)
        return true;
    if (fabs(p.x - problem.get_point(curline.rp).x) < 1.e-8 &&
        fabs(p.y - problem.get_point(curline.rp).y) < 1.e-8)
        return true;
    if (fabs(p.x - problem.get_point(curline.lp).x) < 1.e-8 && curline.vert_line &&
        p.y > problem.get_point(curline.lp).y && p.y < problem.get_point(curline.rp).y)
        return true;
    if (fabs(p.y - problem.get_point(curline.lp).y) < 1.e-8 && !curline.vert_line &&
        p.x > problem.get_point(curline.lp).x && p.x < problem.get_point(curline.rp).x)
        return true;
    return false;
}

#include "task.h"

int task::read_project(const std::string& projName)
{
    std::ifstream file;
    std::string buffer;

    //Чтение файла с геометрией
    file.open(projName);
    if (!file.is_open())
    {
        std::cerr << "Error open file." << std::endl;
        return -7;
    }   
    // Чтение боксов
    double x1, y1, x2, y2, eps;
    int type;
    int pointbuf[4];
    double a;
    point vert;
    int ind;
    while (true)
    {
        file >> buffer;
        if (buffer != "box")
            break;

        boxes.push_back(box());

        file >> x1;
        file >> y1;
        file >> x2;
        file >> y2;
        file >> type;
        file >> eps;

        //Добавление точек
        vert.x = x1;
        vert.y = y1;
        auto it = std::find(verts.begin(), verts.end(), vert);
        if (it == verts.end())
        {
            ind = verts.size();
            verts.push_back(vert);
        }
        else
        {
            ind = it - verts.begin();
        }
        boxes[boxes.size()-1].p[0] = ind;

        vert = { x2, y1 };
        it = std::find(verts.begin(), verts.end(), vert);
        if (it == verts.end())
        {
            ind = verts.size();
            verts.push_back(vert);
        }
        else
        {
            ind = it - verts.begin();
        }
        boxes[boxes.size() - 1].p[1] = ind;

        vert = {x2, y2};
        it = std::find(verts.begin(), verts.end(), vert);
        if (it == verts.end())
        {
            ind = verts.size();
            verts.push_back(vert);
        }
        else
        {
            ind = it - verts.begin();
        }
        boxes[boxes.size() - 1].p[2] = ind;

        vert = { x1, y2 };
        it = std::find(verts.begin(), verts.end(), vert);
        if (it == verts.end())
        {
            ind = verts.size();
            verts.push_back(vert);
        }
        else
        {
            ind = it - verts.begin();
        }
        boxes[boxes.size() - 1].p[3] = ind;
        
        if (type > 0)
            boxes[boxes.size() - 1].eps = NAN;
        else
            boxes[boxes.size() - 1].eps = eps;


        if (eps != eps)
        {
            boxes[boxes.size() - 1].a = 0.0;
        }
        else
        {
            boxes[boxes.size() - 1].a = 10.0;
        }

        buffer = "";
    }
    file.close();
    return 0;
}

int task::fix_boxes()
{
    std::vector<bool> exist;
    exist.resize(boxes.size(), true);
    int changes = 1;
    while (changes != 0)
    {
        changes = 0;
        for (int ibox1 = 0; ibox1 < boxes.size(); ++ibox1)
        {
            if (!exist[ibox1]) continue;
            if (boxes[ibox1].eps != boxes[ibox1].eps) continue;
            for (int ibox2 = 0; ibox2 < boxes.size(); ++ibox2)
            {
                if (!exist[ibox2]) continue;
                if (boxes[ibox1].eps != boxes[ibox2].eps) continue;
                int* p1;
                int* p2;
                p1 = boxes[ibox1].p;
                p2 = boxes[ibox2].p;
                // Общая нижняя грань [для box'a 1]
                if (verts[p1[0]] == verts[p2[3]] && verts[p1[1]] == verts[p2[2]])
                {
                    ++changes;
                    p1[0] = p2[0];
                    p1[1] = p2[1];
                    exist[ibox2] = false;
                    continue;
                }
                // Общая правая грань [для box'a 1]
                if (verts[p1[1]] == verts[p2[0]] && verts[p1[2]] == verts[p2[3]])
                {
                    ++changes;
                    p1[1] = p2[1];
                    p1[2] = p2[2];
                    exist[ibox2] = false;
                    continue;
                }
                // Общая верхняя грань [для box'a 1]
                if (verts[p1[3]] == verts[p2[0]] && verts[p1[2]] == verts[p2[1]])
                {
                    ++changes;
                    p1[2] = p2[2];
                    p1[3] = p2[3];
                    exist[ibox2] = false;
                    continue;
                }
                // Общая левая грань [для box'a 1]
                if (verts[p1[0]] == verts[p2[1]] && verts[p1[3]] == verts[p2[2]])
                {
                    ++changes;
                    p1[0] = p2[0];
                    p1[3] = p2[3];
                    exist[ibox2] = false;
                    continue;
                }
            }
        }
    }

    std::vector<box> boxvrem;
    for (int i = 0; i < boxes.size(); ++i)
    {
        if (exist[i])
        {
            boxvrem.push_back(box());
            for (int j = 0; j < 4; ++j)
                boxvrem[boxvrem.size() - 1].p[j] = boxes[i].p[j];
            boxvrem[boxvrem.size() - 1].eps = boxes[i].eps;
            boxvrem[boxvrem.size() - 1].a = boxes[i].a;
        }
    }

    std::swap(boxes, boxvrem);

    // Чистка от лишних точек
    std::vector<int> existvert;
    existvert.resize(verts.size(), -1);
    int ivert = 0;
    for (int ibox = 0; ibox < boxes.size(); ++ibox)
    {
        int* p = boxes[ibox].p;
        for (int i = 0; i < 4; ++i)
        {
            if (existvert[p[i]] == -1)
            {
                existvert[p[i]] = ivert;
                p[i] = ivert++;
            }
            else
            {
                p[i] = existvert[p[i]];
            }
        }
    }
    std::vector<point> newverts;
    newverts.resize(ivert);
    for (int ivert = 0; ivert < verts.size(); ++ivert)
    {
        if (existvert[ivert] == -1) continue;
        newverts[existvert[ivert]] = verts[ivert];
    }
    std::swap(verts, newverts);
    return 0;
}

int task::save_boxes_geom(const std::string& geomName)
{
    std::ofstream file(geomName);
    if (file.is_open() != true)
    {
        std::cerr << "Can't save geometry" << std::endl;
        return 1;
    }

    auto func = [&](box b1, box b2)
        {
            if (b1.eps != b1.eps)
                return true;
            if (b2.eps != b2.eps)
                return false;
            if (b1.eps < b2.eps)
                return true;
            return false;
        };

    std::sort(boxes.begin(), boxes.end(), func);

    file << "Boxes:\n";
    file << boxes.size() << '\n';
    for (int i = 0; i < boxes.size(); ++i)
    {
        file << boxes[i].p[0] << ' ';
        file << boxes[i].p[1] << ' ';
        file << boxes[i].p[2] << ' ';
        file << boxes[i].p[3] << ' ';
        double eps = boxes[i].eps;
        if (eps != eps)
        {
            file << "None ";
        }
        else
        {
            file << eps << ' ';
        }
        file << boxes[i].a << '\n';
    }
    file << "Verticles:\n";
    file << verts.size() << '\n';
    for (int i = 0; i < verts.size(); ++i)
    {
        file << verts[i].x << ' ' << verts[i].y << '\n';
    }

    file.close();
    return 0;
}

int task::classify_lines()
{
    std::vector<line> right_lines;
    std::vector<line> down_lines;
    std::vector<line> left_lines;
    std::vector <line> up_lines;

    right_lines.resize(boxes.size());
    down_lines.resize(boxes.size());
    left_lines.resize(boxes.size());
    up_lines.resize(boxes.size());

    for (int ibox = 0; ibox < boxes.size(); ++ibox)
    {
        // Левая линия
        line line1;
        line1.vert_line = true;
        line1.lp = boxes[ibox].p[0];
        line1.rp = boxes[ibox].p[3];
        line1.nBox1 = ibox;
        left_lines[ibox] = line1;

        // Нижняя линия
        line line2;
        line2.vert_line = false;
        line2.lp = boxes[ibox].p[0];
        line2.rp = boxes[ibox].p[1];
        line2.nBox1 = ibox;
        down_lines[ibox] = line2;

        // правая линия
        line line3;
        line3.vert_line = true;
        line3.lp = boxes[ibox].p[1];
        line3.rp = boxes[ibox].p[2];
        line3.nBox1 = ibox;
        right_lines[ibox] = line3;

        // Нижняя линия
        line line4;
        line4.vert_line = false;
        line4.lp = boxes[ibox].p[3];
        line4.rp = boxes[ibox].p[2];
        line4.nBox1 = ibox;
        up_lines[ibox] = line4;
    }

    // Поиск правых среди левых
    while (right_lines.size() != 0)
    {
        auto it1 = right_lines.begin();
        bool flag = false;
        for (auto it2 = left_lines.begin(); it2 != left_lines.end(); ++it2)
        {
            if (fabs(verts[it1->rp].x - verts[it2->rp].x) > 1.e-6)
                continue;

            //Эквивалентность
            if(it1->lp == it2->lp && it1->rp == it2->rp)
            {
                it1->nBox2 = it2->nBox1;
                flag = true;
                left_lines.erase(it2);
                vert_lines.push_back(*it1);
                right_lines.erase(it1);
                break;
            }

            // Пересекаются
            if ((verts[it1->rp].y <= verts[it2->rp].y && verts[it1->rp].y >= verts[it2->lp].y) ||
                (verts[it1->lp].y <= verts[it2->rp].y && verts[it1->lp].y >= verts[it2->lp].y))
            {
                if (verts[it1->rp].y == verts[it2->lp].y || verts[it1->lp].y == verts[it2->rp].y)
                    continue;
                int lrp;
                int llp;
                line up;
                line down;
                line cur;
                if (fabs(verts[it1->rp].y - verts[it2->rp].y) > 1.e-6)
                {
                    up.vert_line = true;
                    up.rp = (verts[it1->rp].y > verts[it2->rp].y) ? it1->rp : it2->rp;
                    up.lp = (verts[it1->rp].y > verts[it2->rp].y) ? it2->rp : it1->rp;
                    up.nBox1 = (verts[it1->rp].y > verts[it2->rp].y) ? it1->nBox1 : it2->nBox1;
                    lrp = up.lp;
                }
                else
                {
                    lrp = it1->rp;
                }

                if (fabs(verts[it1->lp].y - verts[it2->lp].y) > 1.e-6)
                {
                    down.vert_line = true;
                    down.rp = (verts[it1->lp].y > verts[it2->lp].y) ? it1->lp : it2->lp;
                    down.lp = (verts[it1->lp].y > verts[it2->lp].y) ? it2->lp : it1-> lp;
                    down.nBox1 = (verts[it1->lp].y > verts[it2->lp].y) ? it2->nBox1 : it1->nBox1;
                    llp = down.rp;
                }
                else
                {
                    llp = it1->lp;
                }

                cur.lp = llp;
                cur.rp = lrp;
                cur.vert_line = true;
                cur.nBox1 = it1->nBox1;
                cur.nBox2 = it2->nBox1;
                vert_lines.push_back(cur);

                right_lines.erase(it1);
                left_lines.erase(it2);

                if (up.lp != up.rp)
                {
                    if (up.nBox1 == cur.nBox1)
                        right_lines.push_back(up);
                    else
                        left_lines.push_back(up);
                }

                if (down.lp != down.rp)
                {
                    if (down.nBox1 == cur.nBox1)
                        right_lines.push_back(down);
                    else
                        left_lines.push_back(down);
                }

                flag = true;
                break;
            }
        }

        // Значит граничное
        if(!flag)
        {
            vert_lines.push_back(*it1);
            right_lines.erase(it1);
        }
    }
    for (auto it2 = left_lines.begin(); it2 != left_lines.end(); ++it2)
    {
        vert_lines.push_back(*it2);
    }
    left_lines.clear();
    
    // Поиск нижних среди верхних
    while (down_lines.size() != 0)
    {
        auto it1 = down_lines.begin();
        bool flag = false;
        for (auto it2 = up_lines.begin(); it2 != up_lines.end(); ++it2)
        {
            if (fabs(verts[it1->rp].y - verts[it2->rp].y) > 1.e-6)
                continue;
            //Эквивалентность
            if (it1->lp == it2->lp && it1->rp == it2->rp)
            {
                it1->nBox2 = it2->nBox1;
                flag = true;
                up_lines.erase(it2);
                horz_lines.push_back(*it1);
                down_lines.erase(it1);
                break;
            }

            // Пересекаются
            if ((verts[it1->rp].x <= verts[it2->rp].x && verts[it1->rp].x >= verts[it2->lp].x) ||
                (verts[it1->lp].x <= verts[it2->rp].x && verts[it1->lp].x >= verts[it2->lp].x))
            {
                if (verts[it1->rp].x == verts[it2->lp].x || verts[it1->lp].x == verts[it2->rp].x)
                    continue;
                int lrp;
                int llp;
                line up;
                line down;
                line cur;
                if (fabs(verts[it1->rp].x - verts[it2->rp].x) > 1.e-6)
                {
                    up.vert_line = false;
                    up.rp = (verts[it1->rp].x > verts[it2->rp].x) ? it1->rp : it2->rp;
                    up.lp = (verts[it1->rp].x > verts[it2->rp].x) ? it2->rp : it1->rp;
                    up.nBox1 = (verts[it1->rp].x > verts[it2->rp].x) ? it1->nBox1 : it2->nBox1;
                    lrp = up.lp;
                }
                else
                {
                    lrp = it1->rp;
                }

                if (fabs(verts[it1->lp].x - verts[it2->lp].x) > 1.e-6)
                {
                    down.vert_line = false;
                    down.rp = (verts[it1->lp].x > verts[it2->lp].x) ? it1->lp : it2->lp;
                    down.lp = (verts[it1->lp].x > verts[it2->lp].x) ? it2->lp : it1->lp;
                    down.nBox1 = (verts[it1->lp].x > verts[it2->lp].x) ? it2->nBox1 : it1->nBox1;
                    llp = down.rp;
                }
                else
                {
                    llp = it1->lp;
                }

                cur.lp = llp;
                cur.rp = lrp;
                cur.vert_line = false;
                cur.nBox1 = it1->nBox1;
                cur.nBox2 = it2->nBox1;
                horz_lines.push_back(cur);

                down_lines.erase(it1);
                up_lines.erase(it2);

                if (up.lp != up.rp)
                {
                    if (up.nBox1 == cur.nBox1)
                        down_lines.push_back(up);
                    else
                        up_lines.push_back(up);
                }

                if (down.lp != down.rp)
                {
                    if (down.nBox1 == cur.nBox1)
                        down_lines.push_back(down);
                    else
                        up_lines.push_back(down);
                }

                flag = true;
                break;
            }
        }

        // Значит граничное
        if (!flag)
        {
            horz_lines.push_back(*it1);
            down_lines.erase(it1);
        }
    }
    for (auto it2 = up_lines.begin(); it2 != up_lines.end(); ++it2)
    {
        horz_lines.push_back(*it2);
    }
    up_lines.clear();

    return 0;
}

void task::debug_save_lines()
{
    std::ofstream file("C:\\Users\\egora\\Desktop\\Aspa\\Codes\\General_progs\\test_lines.txt");
    for (line a : vert_lines)
    {
        file << verts[a.lp].x << " " << verts[a.lp].y << " " << verts[a.rp].x << " " << verts[a.rp].y << "\n";
    }
    for (line a : horz_lines)
    {
        file << verts[a.lp].x << " " << verts[a.lp].y << " " << verts[a.rp].x << " " << verts[a.rp].y << "\n";
    }
    file.close();
}

int task::box_processing(int n, int m)
{
    // начальная инициализация
    const int conrazb = n;
    const int layers  = m;

    std::vector<box> newBoxes;
    std::vector<point> newVerts;

    auto addgetvert = [&](point  node)
        {
            auto it = std::find(verts.begin(), verts.end(), node);
            if (it == verts.end())
            {
                verts.push_back(node);
                return int(verts.size() - 1);
            }
            return int(it - verts.begin());
        };

    for (int ibox = 0; ibox < boxes.size(); ++ibox)
    {
        if (boxes[ibox].eps != boxes[ibox].eps)
        {
            int* p = boxes[ibox].p;
            double lx = verts[p[1]].x - verts[p[0]].x;
            double ly = verts[p[2]].y - verts[p[1]].y;
            double a = std::fmin(lx / conrazb, ly / conrazb);
            boxes[ibox].a = a;
            boxes[ibox].que = 0;
        }
        else
        {
            boxes[ibox].que = -1;
        }
    }

        double xmax, ymax;
        double xmin, ymin;
        xmax = get_max_x();
        xmin = get_min_x();
        ymax = get_max_y();
        ymin = get_min_y();

        int curque = 0;
        bool flag = false;
        std::vector<int> for_delete;
        while (true)
        {
            for (int ibox = 0; ibox < boxes.size(); ibox++)
            {
                if (boxes[ibox].get_que() == curque)
                {
                    double a = boxes[ibox].get_a() * 2.0;
                    // 1. Обработка нижней границе [т.е ищем среди верхних дургих элементов]
                    int a1 = boxes[ibox].get_points()[0];
                    int b1 = boxes[ibox].get_points()[1];
                    double y1 = verts[a1].y;
                    int boxessize = boxes.size();
                    if (fabs(y1 - ymin) > 1.e-6)
                    {
                        double y2 = y1 - a * 15;
                        for (int iboxdop = 0; iboxdop < boxessize; ++iboxdop)
                        {

                            if (iboxdop == ibox) continue;
                            if (abs(y1 - verts[boxes[iboxdop].get_points()[2]].y) > 1.e-6)
                                continue;
                            if (verts[boxes[iboxdop].get_points()[2]].x - 1.e-6 < verts[a1].x)
                                continue;
                            if (verts[boxes[iboxdop].get_points()[3]].x > verts[b1].x - 1.e-6)
                                continue;
                            if (a + 1.e-6 > boxes[iboxdop].get_a())
                                continue;
                            // Сюда дошло значит пересекаются
                            int a2 = boxes[iboxdop].get_points()[3];
                            int b2 = boxes[iboxdop].get_points()[2];
                            int per_case = 0;
                            // 1. вариант пересечения [a2 < a1  b2 < b1]
                            if (verts[a2].x < verts[a1].x && verts[b2].x < verts[b1].x)
                                per_case = 1;
                            if (verts[a1].x < verts[a2].x && verts[b1].x < verts[b2].x)
                                per_case = 2;
                            if (fabs(verts[a1].x - verts[a2].x) < 1.e-6 && abs(verts[b1].x - verts[b2].x) < 1.e-6)
                                per_case = 3;
                            if (verts[a2].x < verts[a1].x && verts[b1].x < verts[b2].x)
                                per_case = 4;
                            if (fabs(verts[a1].x - verts[a2].x) < 1.e-6 && verts[b2].x < verts[b1].x)
                                per_case = 5;
                            if (verts[a1].x < verts[a2].x && abs(verts[b1].x - verts[b2].x) < 1.e-6)
                                per_case = 6;
                            if (verts[a1].x < verts[a2].x && verts[b2].x < verts[b1].x)
                                per_case = 7;
                            if (fabs(verts[a1].x - verts[a2].x) < 1.e-6 && verts[b1].x < verts[b2].x)
                                per_case = 8;
                            if (verts[a2].x < verts[a1].x && abs(verts[b1].x - verts[b2].x) < 1.e-6)
                                per_case = 9;

                            bool fully = ((y2 - a) <= verts[boxes[iboxdop].get_points()[1]].y) ? true : false;

                            point v1;
                            v1.init(std::max(verts[a1].x, verts[a2].x), verts[boxes[iboxdop].get_points()[0]].y);
                            point v2;
                            v2.init(std::max(verts[a1].x, verts[a2].x), (fully) ? verts[boxes[iboxdop].get_points()[0]].y : y2);
                            point v3;
                            v3.init(verts[a2].x, (fully) ? verts[boxes[iboxdop].get_points()[0]].y : y2);
                            point v4;
                            v4.init(std::min(verts[b1].x, verts[b2].x), verts[boxes[iboxdop].get_points()[0]].y);
                            point v5;
                            v5.init(std::min(verts[b1].x, verts[b2].x), (fully) ? verts[boxes[iboxdop].get_points()[0]].y : y2);
                            point v6;
                            v6.init(verts[b2].x, (fully) ? verts[boxes[iboxdop].get_points()[0]].y : y2);
                            point v7;
                            v7.init(std::max(verts[a1].x, verts[a2].x), verts[boxes[iboxdop].get_points()[2]].y);
                            point v8;
                            v8.init(std::min(verts[b1].x, verts[b2].x), verts[boxes[iboxdop].get_points()[2]].y);

                            // Левые боковые
                            if (per_case == 1 || per_case == 4 || per_case == 9)
                            {
                                box curbox1;
                                int v[4];
                                v[0] = addgetvert(v3);
                                v[1] = addgetvert(v2);
                                v[2] = addgetvert(v7);
                                v[3] = boxes[iboxdop].get_points()[3];
                                curbox1.init(v, boxes[iboxdop].eps, boxes[iboxdop].get_a());
                                boxes.push_back(curbox1);
                                if (!fully)
                                {
                                    v[0] = boxes[iboxdop].get_points()[0];
                                    v[1] = addgetvert(v1);
                                    v[2] = addgetvert(v2);
                                    v[3] = addgetvert(v3);
                                    curbox1.init(v, boxes[iboxdop].eps, boxes[iboxdop].get_a());
                                    boxes.push_back(curbox1);
                                }
                            }

                            // Правые боковые
                            if (per_case == 2 || per_case == 4 || per_case == 8)
                            {
                                box curbox1;
                                int v[4];
                                v[0] = addgetvert(v5);
                                v[1] = addgetvert(v6);
                                v[2] = boxes[iboxdop].get_points()[2];
                                v[3] = addgetvert(v8);
                                curbox1.init(v, boxes[iboxdop].eps, boxes[iboxdop].get_a());
                                boxes.push_back(curbox1);
                                if (!fully)
                                {
                                    v[0] = addgetvert(v4);
                                    v[1] = boxes[iboxdop].get_points()[1];
                                    v[2] = addgetvert(v6);
                                    v[3] = addgetvert(v5);
                                    curbox1.init(v, boxes[iboxdop].eps, boxes[iboxdop].get_a());
                                    boxes.push_back(curbox1);
                                }
                            }

                            // самые важные
                            box curbox2;
                            int v[4];
                            v[0] = addgetvert(v2);
                            v[1] = addgetvert(v5);
                            v[2] = addgetvert(v8);
                            v[3] = addgetvert(v7);
                            curbox2.init(v, boxes[iboxdop].eps, a);
                            curbox2.set_que(curque + 1);
                            boxes.push_back(curbox2);
                            if (!fully)
                            {
                                v[0] = addgetvert(v1);
                                v[1] = addgetvert(v4);
                                v[2] = addgetvert(v5);
                                v[3] = addgetvert(v2);
                                curbox2.init(v, boxes[iboxdop].eps, boxes[iboxdop].get_a());
                                boxes.push_back(curbox2);
                            }

                            for_delete.push_back(iboxdop);

                        }
                    }

                    // 2. Обработка правой границы
                    a1 = boxes[ibox].get_points()[1];
                    b1 = boxes[ibox].get_points()[2];
                    double x1 = verts[a1].x;
                    if (fabs(x1 - xmax) > 1.e-6)
                    {
                        double x2 = x1 + 15 * a;
                        for (int iboxdop = 0; iboxdop < boxessize; ++iboxdop)
                        {
                            if (iboxdop == ibox) continue;
                            if (abs(x1 - verts[boxes[iboxdop].get_points()[0]].x) > 1.e-6)
                                continue;
                            if (verts[boxes[iboxdop].get_points()[3]].y - 1.e-6 < verts[a1].y)
                                continue;
                            if (verts[boxes[iboxdop].get_points()[0]].y > verts[b1].y - 1.e-6)
                                continue;
                            if (a + 1.e-6 > boxes[iboxdop].get_a())
                                continue;
                            // Сюда дошло значит пересекаются
                            int a2 = boxes[iboxdop].get_points()[0];
                            int b2 = boxes[iboxdop].get_points()[3];
                            int per_case = 0;
                            // 1. вариант пересечения [a2 < a1  b2 < b1]
                            if (verts[a2].y < verts[a1].y && verts[b2].y < verts[b1].y)
                                per_case = 1;
                            if (verts[a1].y < verts[a2].y && verts[b1].y < verts[b2].y)
                                per_case = 2;
                            if (fabs(verts[a1].y - verts[a2].y) < 1.e-6 && abs(verts[b1].y - verts[b2].y) < 1.e-6)
                                per_case = 3;
                            if (verts[a2].y < verts[a1].y && verts[b1].y < verts[b2].y)
                                per_case = 4;
                            if (fabs(verts[a1].y - verts[a2].y) < 1.e-6 && verts[b2].y < verts[b1].y)
                                per_case = 5;
                            if (verts[a1].y < verts[a2].y && abs(verts[b1].y - verts[b2].y) < 1.e-6)
                                per_case = 6;
                            if (verts[a1].y < verts[a2].y && verts[b2].y < verts[b1].y)
                                per_case = 7;
                            if (fabs(verts[a1].y - verts[a2].y) < 1.e-6 && verts[b1].y < verts[b2].y)
                                per_case = 8;
                            if (verts[a2].y < verts[a1].y && abs(verts[b1].y - verts[b2].y) < 1.e-6)
                                per_case = 9;

                            bool fully = ((x2 + a) >= verts[boxes[iboxdop].get_points()[1]].x) ? true : false;

                            point v1;
                            v1.init((fully) ? verts[boxes[iboxdop].get_points()[1]].x : x2, verts[boxes[iboxdop].get_points()[1]].y);
                            point v2;
                            v2.init(verts[a1].x, std::max(verts[a1].y, verts[a2].y));
                            point v3;
                            v3.init((fully) ? verts[boxes[iboxdop].get_points()[1]].x : x2, std::max(verts[a1].y, verts[a2].y));
                            point v4;
                            v4.init(verts[boxes[iboxdop].get_points()[1]].x, std::max(verts[a1].y, verts[a2].y));
                            point v5;
                            v5.init(verts[a1].x, std::min(verts[b1].y, verts[b2].y));
                            point v6;
                            v6.init((fully) ? verts[boxes[iboxdop].get_points()[1]].x : x2, std::min(verts[b1].y, verts[b2].y));
                            point v7;
                            v7.init(verts[boxes[iboxdop].get_points()[1]].x, std::min(verts[b1].y, verts[b2].y));
                            point v8;
                            v8.init((fully) ? verts[boxes[iboxdop].get_points()[1]].x : x2, verts[b2].y);

                            // Нижние
                            if (per_case == 1 || per_case == 4 || per_case == 9)
                            {
                                box curbox1;
                                int v[4];
                                v[0] = boxes[iboxdop].get_points()[0];
                                v[1] = addgetvert(v1);
                                v[2] = addgetvert(v3);
                                v[3] = addgetvert(v2);
                                curbox1.init(v, boxes[iboxdop].eps, boxes[iboxdop].get_a());
                                boxes.push_back(curbox1);
                                if (!fully)
                                {
                                    v[0] = addgetvert(v1);
                                    v[1] = boxes[iboxdop].get_points()[1];
                                    v[2] = addgetvert(v4);
                                    v[3] = addgetvert(v3);
                                    curbox1.init(v, boxes[iboxdop].eps, boxes[iboxdop].get_a());
                                    boxes.push_back(curbox1);
                                }
                            }

                            // верхние 
                            if (per_case == 2 || per_case == 4 || per_case == 8)
                            {
                                box curbox1;
                                int v[4];
                                v[0] = addgetvert(v5);
                                v[1] = addgetvert(v6);
                                v[2] = addgetvert(v8);
                                v[3] = boxes[iboxdop].get_points()[3];
                                curbox1.init(v, boxes[iboxdop].eps, boxes[iboxdop].get_a());
                                boxes.push_back(curbox1);
                                if (!fully)
                                {
                                    v[0] = addgetvert(v6);
                                    v[1] = addgetvert(v7);
                                    v[2] = boxes[iboxdop].get_points()[2];
                                    v[3] = addgetvert(v8);
                                    curbox1.init(v, boxes[iboxdop].eps, boxes[iboxdop].get_a());
                                    boxes.push_back(curbox1);
                                }
                            }

                            // самые важные
                            box curbox2;
                            int v[4];
                            v[0] = addgetvert(v2);
                            v[1] = addgetvert(v3);
                            v[2] = addgetvert(v6);
                            v[3] = addgetvert(v5);
                            curbox2.init(v, boxes[iboxdop].eps, a);
                            curbox2.set_que(curque + 1);
                            boxes.push_back(curbox2);
                            if (!fully)
                            {
                                v[0] = addgetvert(v3);
                                v[1] = addgetvert(v4);
                                v[2] = addgetvert(v7);
                                v[3] = addgetvert(v6);
                                curbox2.init(v, boxes[iboxdop].eps, boxes[iboxdop].get_a());
                                boxes.push_back(curbox2);
                            }

                            for_delete.push_back(iboxdop);
                        }
                    }

                    // 3. Обработка верхней границы
                    a1 = boxes[ibox].get_points()[3];
                    b1 = boxes[ibox].get_points()[2];
                    y1 = verts[a1].y;
                    if (fabs(y1 - ymax) > 1.e-6)
                    {
                        double y2 = y1 + a * 15;
                        for (int iboxdop = 0; iboxdop < boxessize; ++iboxdop)
                        {
                            if (iboxdop == ibox) continue;
                            if (abs(y1 - verts[boxes[iboxdop].get_points()[0]].y) > 1.e-6)
                                continue;
                            if (verts[boxes[iboxdop].get_points()[2]].x - 1.e-6 < verts[a1].x)
                                continue;
                            if (verts[boxes[iboxdop].get_points()[3]].x > verts[b1].x - 1.e-6)
                                continue;
                            if (a + 1.e-6 > boxes[iboxdop].get_a())
                                continue;
                            // Сюда дошло значит пересекаются
                            int a2 = boxes[iboxdop].get_points()[0];
                            int b2 = boxes[iboxdop].get_points()[1];
                            int per_case = 0;
                            // 1. вариант пересечения [a2 < a1  b2 < b1]
                            if (verts[a2].x < verts[a1].x && verts[b2].x < verts[b1].x)
                                per_case = 1;
                            if (verts[a1].x < verts[a2].x && verts[b1].x < verts[b2].x)
                                per_case = 2;
                            if (fabs(verts[a1].x - verts[a2].x) < 1.e-6 && abs(verts[b1].x - verts[b2].x) < 1.e-6)
                                per_case = 3;
                            if (verts[a2].x < verts[a1].x && verts[b1].x < verts[b2].x)
                                per_case = 4;
                            if (fabs(verts[a1].x - verts[a2].x) < 1.e-6 && verts[b2].x < verts[b1].x)
                                per_case = 5;
                            if (verts[a1].x < verts[a2].x && abs(verts[b1].x - verts[b2].x) < 1.e-6)
                                per_case = 6;
                            if (verts[a1].x < verts[a2].x && verts[b2].x < verts[b1].x)
                                per_case = 7;
                            if (fabs(verts[a1].x - verts[a2].x) < 1.e-6 && verts[b1].x < verts[b2].x)
                                per_case = 8;
                            if (verts[a2].x < verts[a1].x && abs(verts[b1].x - verts[b2].x) < 1.e-6)
                                per_case = 9;


                            bool fully = ((y2 + a) >= verts[boxes[iboxdop].get_points()[2]].y) ? true : false;

                            point v1;
                            v1.init(std::max(verts[a1].x, verts[a2].x), verts[a1].y);
                            point v2;
                            v2.init(std::min(verts[b1].x, verts[b2].x), verts[a1].y);
                            point v3;
                            v3.init(verts[a2].x, (fully) ? verts[boxes[iboxdop].get_points()[2]].y : y2);
                            point v4;
                            v4.init(std::max(verts[a1].x, verts[a2].x), (fully) ? verts[boxes[iboxdop].get_points()[2]].y : y2);
                            point v5;
                            v5.init(std::min(verts[b1].x, verts[b2].x), (fully) ? verts[boxes[iboxdop].get_points()[2]].y : y2);
                            point v6;
                            v6.init(verts[b2].x, (fully) ? verts[boxes[iboxdop].get_points()[2]].y : y2);
                            point v7;
                            v7.init(std::max(verts[a1].x, verts[a2].x), verts[boxes[iboxdop].get_points()[2]].y);
                            point v8;
                            v8.init(std::min(verts[b1].x, verts[b2].x), verts[boxes[iboxdop].get_points()[2]].y);

                            // Левые боковые
                            if (per_case == 1 || per_case == 4 || per_case == 9)
                            {
                                box curbox1;
                                int v[4];
                                v[0] = boxes[iboxdop].get_points()[0];
                                v[1] = addgetvert(v1);
                                v[2] = addgetvert(v4);
                                v[3] = addgetvert(v3);
                                curbox1.init(v, boxes[iboxdop].eps, boxes[iboxdop].get_a());
                                boxes.push_back(curbox1);
                                if (!fully)
                                {
                                    v[0] = addgetvert(v3);
                                    v[1] = addgetvert(v4);
                                    v[2] = addgetvert(v7);
                                    v[3] = boxes[iboxdop].get_points()[3];
                                    curbox1.init(v, boxes[iboxdop].eps, boxes[iboxdop].get_a());
                                    boxes.push_back(curbox1);
                                }
                            }

                            // Правые боковые
                            if (per_case == 2 || per_case == 4 || per_case == 8)
                            {
                                box curbox1;
                                int v[4];
                                v[0] = addgetvert(v2);
                                v[1] = boxes[iboxdop].get_points()[1];
                                v[2] = addgetvert(v6);
                                v[3] = addgetvert(v5);
                                curbox1.init(v, boxes[iboxdop].eps, boxes[iboxdop].get_a());
                                boxes.push_back(curbox1);
                                if (!fully)
                                {
                                    v[0] = addgetvert(v5);
                                    v[1] = addgetvert(v6);
                                    v[2] = boxes[iboxdop].get_points()[2];
                                    v[3] = addgetvert(v8);
                                    curbox1.init(v, boxes[iboxdop].eps, boxes[iboxdop].get_a());
                                    boxes.push_back(curbox1);
                                }
                            }

                            // самые важные
                            box curbox2;
                            int v[4];
                            v[0] = addgetvert(v1);
                            v[1] = addgetvert(v2);
                            v[2] = addgetvert(v5);
                            v[3] = addgetvert(v4);
                            curbox2.init(v, boxes[iboxdop].eps, a);
                            curbox2.set_que(curque + 1);
                            boxes.push_back(curbox2);
                            if (!fully)
                            {
                                v[0] = addgetvert(v4);
                                v[1] = addgetvert(v5);
                                v[2] = addgetvert(v8);
                                v[3] = addgetvert(v7);
                                curbox2.init(v, boxes[iboxdop].eps, boxes[iboxdop].get_a());
                                boxes.push_back(curbox2);
                            }

                            for_delete.push_back(iboxdop);
                        }
                    }

                    // 4. Обработка левой стенки
                    a1 = boxes[ibox].get_points()[0];
                    b1 = boxes[ibox].get_points()[3];
                    x1 = verts[a1].x;
                    if (fabs(x1 - xmin) > 1.e-6)
                    {
                        double x2 = x1 - 15 * a;
                        for (int iboxdop = 0; iboxdop < boxessize; ++iboxdop)
                        {
                            if (iboxdop == ibox) continue;
                            if (abs(x1 - verts[boxes[iboxdop].get_points()[1]].x) > 1.e-6)
                                continue;
                            if (verts[boxes[iboxdop].get_points()[3]].y - 1.e-6 < verts[a1].y)
                                continue;
                            if (verts[boxes[iboxdop].get_points()[0]].y > verts[b1].y - 1.e-6)
                                continue;
                            if (a + 1.e-6 > boxes[iboxdop].get_a())
                                continue;
                            // Сюда дошло значит пересекаются
                            int a2 = boxes[iboxdop].get_points()[1];
                            int b2 = boxes[iboxdop].get_points()[2];
                            int per_case = 0;
                            // 1. вариант пересечения [a2 < a1  b2 < b1]
                            if (verts[a2].y < verts[a1].y && verts[b2].y < verts[b1].y)
                                per_case = 1;
                            if (verts[a1].y < verts[a2].y && verts[b1].y < verts[b2].y)
                                per_case = 2;
                            if (fabs(verts[a1].y - verts[a2].y) < 1.e-6 && abs(verts[b1].y - verts[b2].y) < 1.e-6)
                                per_case = 3;
                            if (verts[a2].y < verts[a1].y && verts[b1].y < verts[b2].y)
                                per_case = 4;
                            if (fabs(verts[a1].y - verts[a2].y) < 1.e-6 && verts[b2].y < verts[b1].y)
                                per_case = 5;
                            if (verts[a1].y < verts[a2].y && abs(verts[b1].y - verts[b2].y) < 1.e-6)
                                per_case = 6;
                            if (verts[a1].y < verts[a2].y && verts[b2].y < verts[b1].y)
                                per_case = 7;
                            if (fabs(verts[a1].y - verts[a2].y) < 1.e-6 && verts[b1].y < verts[b2].y)
                                per_case = 8;
                            if (verts[a2].y < verts[a1].y && abs(verts[b1].y - verts[b2].y) < 1.e-6)
                                per_case = 9;

                            bool fully = ((x2 - a) <= verts[boxes[iboxdop].get_points()[0]].x) ? true : false;

                            point v1;
                            v1.init((fully) ? verts[boxes[iboxdop].get_points()[0]].x : x2, verts[a2].y);
                            point v2;
                            v2.init(verts[boxes[iboxdop].get_points()[0]].x, std::max(verts[a1].y, verts[a2].y));
                            point v3;
                            v3.init((fully) ? verts[boxes[iboxdop].get_points()[0]].x : x2, std::max(verts[a1].y, verts[a2].y));
                            point v4;
                            v4.init(verts[a1].x, std::max(verts[a1].y, verts[a2].y));
                            point v5;
                            v5.init(verts[boxes[iboxdop].get_points()[0]].x, std::min(verts[b1].y, verts[b2].y));
                            point v6;
                            v6.init((fully) ? verts[boxes[iboxdop].get_points()[0]].x : x2, std::min(verts[b1].y, verts[b2].y));
                            point v7;
                            v7.init(verts[a1].x, std::min(verts[b1].y, verts[b2].y));
                            point v8;
                            v8.init((fully) ? verts[boxes[iboxdop].get_points()[0]].x : x2, verts[b2].y);

                            // Нижние
                            if (per_case == 1 || per_case == 4 || per_case == 9)
                            {
                                box curbox1;
                                int v[4];
                                v[0] = addgetvert(v1);
                                v[1] = boxes[iboxdop].get_points()[1];
                                v[2] = addgetvert(v4);
                                v[3] = addgetvert(v3);
                                curbox1.init(v, boxes[iboxdop].eps, boxes[iboxdop].get_a());
                                boxes.push_back(curbox1);
                                if (!fully)
                                {
                                    v[0] = boxes[iboxdop].get_points()[0];
                                    v[1] = addgetvert(v1);
                                    v[2] = addgetvert(v3);
                                    v[3] = addgetvert(v2);
                                    curbox1.init(v, boxes[iboxdop].eps, boxes[iboxdop].get_a());
                                    boxes.push_back(curbox1);
                                }
                            }

                            // верхние 
                            if (per_case == 2 || per_case == 4 || per_case == 8)
                            {
                                box curbox1;
                                int v[4];
                                v[0] = addgetvert(v6);
                                v[1] = addgetvert(v7);
                                v[2] = boxes[iboxdop].get_points()[2];
                                v[3] = addgetvert(v8);
                                curbox1.init(v, boxes[iboxdop].eps, boxes[iboxdop].get_a());
                                boxes.push_back(curbox1);
                                if (!fully)
                                {
                                    v[0] = addgetvert(v5);
                                    v[1] = addgetvert(v6);
                                    v[2] = addgetvert(v8);
                                    v[3] = boxes[iboxdop].get_points()[3];
                                    curbox1.init(v, boxes[iboxdop].eps, boxes[iboxdop].get_a());
                                    boxes.push_back(curbox1);
                                }
                            }

                            // самые важные
                            box curbox2;
                            int v[4];
                            v[0] = addgetvert(v3);
                            v[1] = addgetvert(v4);
                            v[2] = addgetvert(v7);
                            v[3] = addgetvert(v6);
                            curbox2.init(v, boxes[iboxdop].eps, a);
                            curbox2.set_que(curque + 1);
                            boxes.push_back(curbox2);
                            if (!fully)
                            {
                                v[0] = addgetvert(v2);
                                v[1] = addgetvert(v3);
                                v[2] = addgetvert(v6);
                                v[3] = addgetvert(v5);
                                curbox2.init(v, boxes[iboxdop].eps, boxes[iboxdop].get_a());
                                boxes.push_back(curbox2);
                            }

                            for_delete.push_back(iboxdop);
                        }
                    }

                    if (for_delete.size() != 0)
                    {
                        boxes[ibox].set_que(-1);
                        if (boxes[ibox].eps != boxes[ibox].eps)
                            boxes[ibox].set_a(0.0);
                        break;
                    }

                }
            }

            if (for_delete.size() == 0 && flag)
                break;
            else if (for_delete.size() == 0)
            {
                curque++;
                flag = true;
            }
            else if (flag)
            {
                flag = false;
            }

            std::sort(for_delete.rbegin(), for_delete.rend());
            for (int i = 0; i < for_delete.size(); ++i)
            {
                auto it = boxes.begin() + for_delete[i];
                boxes.erase(it);
            }

            for_delete.clear();
        }

    return 0;
}

box task::get_box(int id)
{
    return boxes[id];
}

int task::get_box_size()
{
    return boxes.size();
}

point task::get_point(int id)
{
    return verts[id];
}

line task::get_vert_line(int id)
{
    return vert_lines[id];
}

line task::get_horz_line(int id)
{
    return horz_lines[id];
}

int task::get_vert_lines_size()
{
    return vert_lines.size();
}

int task::get_horz_lines_size()
{
    return horz_lines.size();
}

double task::get_max_x()
{
    double maxx = verts[0].x;
    for (int ibox = 0; ibox < boxes.size(); ++ibox)
    {
        if (maxx < verts[boxes[ibox].p[1]].x)
        {
            maxx = verts[boxes[ibox].p[1]].x;
        }
    }
    return maxx;
}

double task::get_min_x()
{
    double minx = verts[0].x;
    for (int ibox = 0; ibox < boxes.size(); ++ibox)
    {
        if (minx > verts[boxes[ibox].p[0]].x)
        {
            minx = verts[boxes[ibox].p[0]].x;
        }
    }
    return minx;
}

double task::get_max_y()
{
    double maxy = verts[0].y;
    for (int ibox = 0; ibox < boxes.size(); ++ibox)
    {
        if (maxy < verts[boxes[ibox].p[2]].y)
        {
            maxy = verts[boxes[ibox].p[2]].y;
        }
    }
    return maxy;
}

double task::get_min_y()
{
    double miny = verts[0].y;
    for (int ibox = 0; ibox < boxes.size(); ++ibox)
    {
        if (miny > verts[boxes[ibox].p[0]].y)
        {
            miny = verts[boxes[ibox].p[0]].y;
        }
    }
    return miny;
}

int task::get_nets_size()
{
    int res = 0;
    for (int ibox = 0; ibox < boxes.size(); ++ibox)
        if (boxes[ibox].is_conductor())
            res++;
    return res;
}

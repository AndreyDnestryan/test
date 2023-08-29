#include <iostream>
#include <vector>

using uint = unsigned int;

using namespace std;

template <typename T> T sqr(T a) {
    return a*a;
}

/*template <typename T> T min(T a, T b) {
    return a < b ? a : b;
}*/ 

struct Circle { //Создаем класс с функциями
    double x, y, r;
    
    bool is_point_inside(double px, double py);
    
    double min_x() { return x - r; };
    double max_x() { return x + r; };
    double min_y() { return y - r; };
    double max_y() { return y + r; };
};

bool Circle::is_point_inside(double px, double py) {
    if (sqr(x-px) + sqr(y-py) <= sqr(r)) return true;
    return false;
}

int main() {
    //const double RECT_SIDE = 1e-4; // размер ячейки
    //const double RECT_AREA = sqr(RECT_SIDE); // прощадь ячейки
    
    uint N;
    vector<Circle> circles;
    double total_area = 0.;
    
    cin >> N; // считываем входные данные
    circles.resize(N);
    double sum = 0.;
    for (auto &c : circles) {
        cin >> c.x >> c.y >> c.r;
		sum += c.r;
    }
	
    double RECT_SIDE = 1e-5/(10.*sum); // 2*pi*sqrt(2)*sum < epsilon = 1e-5
    double RECT_AREA = sqr(RECT_SIDE);
    
    double min_x = circles[0].min_x(), min_y = circles[0].min_y();
    double max_x = circles[0].max_x(), max_y = circles[0].max_y();
    for (uint i = 1; i < circles.size(); i++) { // выбираем прямоугольник, содержащий все круги
        auto &c = circles[i];            
        min_x = min(min_x, c.min_x());
        max_x = max(max_x, c.max_x());
        min_y = min(min_y, c.min_y());
        max_y = max(max_y, c.max_y());
        //cout << min_x << ", " << max_x << "; " << min_y << ", " << max_y << endl;
    }
    
    for (double x = min_x - RECT_SIDE/2; x <= max_x + RECT_SIDE/2; x += RECT_SIDE) {
        for (double y = min_y - RECT_SIDE/2; y <= max_y + RECT_SIDE/2; y += RECT_SIDE) {
            for (auto c : circles) {
                if (c.is_point_inside(x, y)) { // проверяем ячейки на принадлежность кругам
                    total_area += RECT_AREA;
                    break;
                }
            }
        }
    }
    
    cout << "Total area: " << total_area << endl;
    return 0;
}

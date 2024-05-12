#include <iostream>

int main(int argc, char** argv) {
    unsigned int n_x, n_y;
    short int I_max;
    double x_L, y_L, x_R, y_R;

    if (argc == 8) {
        n_x = static_cast<unsigned int>(std::strtol(argv[1], nullptr, 10));
        n_y = static_cast<unsigned int>(std::strtol(argv[2], nullptr, 10));
        x_L = std::strtod(argv[3], nullptr);
        y_L = std::strtod(argv[4], nullptr);
        x_R = std::strtod(argv[5], nullptr);
        y_R = std::strtod(argv[6], nullptr);
        I_max = static_cast<unsigned int>(std::strtol(argv[7], nullptr, 10));
    }
    else {
        std::cerr << "Incorrect number of parameters. Expected 7 parameters." << std::endl;
    }

    return 0;
}


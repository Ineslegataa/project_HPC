#include <iostream>
#include <cstdlib>
#include <complex>
#include <vector>
#include <cstring>

extern "C" void write_pgm_image(void* image, int maxval, int xsize, int ysize, const char* image_name);

class Mandelbrot_set {
private:
    // aggiungere valori di default
    unsigned int n_x, n_y;
    short int I_max=0;
    double x_L, y_L, x_R, y_R;
    double delta_x,delta_y;
    std::vector<std::complex<double>> grid;

    void initialize_grid() {
    delta_x = (x_R - x_L) / (n_x);
    delta_y = (y_R - y_L) / (n_y);
    grid.resize(n_x * n_y);
    for (size_t i = 0; i < n_x; ++i) {
        for (size_t j = 0; j < n_y; ++j) {
            double real_part = x_L + i * delta_x;
            double imag_part = y_L + j * delta_y;
            grid[i * n_y + j] = std::complex<double>(real_part, imag_part); // Access using row-major order
        }
    }
    }


public:
    Mandelbrot_set(unsigned int nx, unsigned int ny, double xl, double yl, double xr, double yr, short int imax)
        : n_x(nx), n_y(ny), x_L(xl), y_L(yl), x_R(xr), y_R(yr), I_max(imax) {initialize_grid();}

    void display_parameters() const {
        std::cout << "Resolution: " << n_x << "x" << n_y << std::endl;
        std::cout << "Bounds: [" << x_L << ", " << y_L << "] to [" << x_R << ", " << y_R << "]" << std::endl;
        std::cout << "Maximum iterations: " << I_max << std::endl;
    }

    int pixel_color(size_t I, size_t J) const {

        std::complex<double> c = grid[I * n_y + J];
        std::complex<double> z = 0;
        bool check_mandelbrot = true;

        for (int n = 0; n < I_max; ++n) {
            if (std::norm(z) > 4) return n;
            z = z * z + c;
        }
        return 0;

    }

    std::vector<unsigned char> generate_image() {
        std::vector<unsigned char> image(n_x * n_y);
        for (size_t i = 0; i < n_x; ++i) {
            for (size_t j = 0; j < n_y; ++j) {
                image[i * n_y + j] = static_cast<unsigned char>(pixel_color(i, j));
            }
        }
        return image;
    }

};

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

    Mandelbrot_set mandelbrot_set(n_x,n_y,x_L,y_L,x_R,y_R,I_max);
    mandelbrot_set.display_parameters();
    auto image = mandelbrot_set.generate_image();
    write_pgm_image(image.data(), 255, n_x, n_y, "mandelbrot.pgm");

    return 0;
}


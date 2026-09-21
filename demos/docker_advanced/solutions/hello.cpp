#include <iostream>

int main(int argc, char **argv) {
    const char *message = argc > 1 ? argv[1] : "Hello from a Setonix-ready image";
    std::cout << message << '\n';
    return 0;
}

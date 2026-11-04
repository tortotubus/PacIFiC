#include "editor.hpp"
#include <libxml/parser.h>
#include <iostream>

/**
 * @brief Run the interactive editor or perform a read-only document
 * validation.
 */
int main(int argc, char** argv) {
    if (argc == 2 && std::string(argv[1]) == "--help") {
        std::cout << "Usage: grains-xml [--check] document.xml\n"
                     "Grains3D terminal XML editor. Missing paths start new Grains3D documents. See Grains/Tools/XML/README.md for keys and scope.\n";
        return 0;
    }
    bool check = argc == 3 && std::string(argv[1]) == "--check";
    if (argc != 2 && !check) { std::cerr << "Usage: grains-xml [--check] document.xml\n"; return 2; }
    xmlInitParser();
    int result = 0;
    try {
        editor::Document document(argv[check ? 2 : 1], !check);
        if (check) { result = document.validate() ? 0 : 1; std::cout << document.message << '\n'; }
        else { editor::Editor application(document); application.run(); }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n'; result = 1;
    }
    xmlCleanupParser();
    return result;
}

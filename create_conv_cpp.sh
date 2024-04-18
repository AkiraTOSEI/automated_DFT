#!/bin/bash

# INCARファイルのパスを指定
CPP_FILE="convergence_checker.cpp"
L2_NORM_FILE="l2_norm.cpp"
CONV_SUMMARY="convergence_summary.cpp"

# INCARファイルの内容を定義
cat > $CPP_FILE << EOF
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>

// Function to parse the times_val which could be a fraction or a float
float parseFractionOrFloat(const std::string& input) {
    size_t slashPos = input.find('/');
    if (slashPos != std::string::npos) {
        // Input is a fraction, split it into numerator and denominator
        float numerator = std::stof(input.substr(0, slashPos));
        float denominator = std::stof(input.substr(slashPos + 1));
        return numerator / denominator;  // Return the result of the division
    }
    // Otherwise, assume it is a float
    return std::stof(input);
}

int main(int argc, char* argv[]) {
    // Verify the correct number of command line arguments are provided
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <pressure.dat> <encut.dat> <times_val> <criteria_val> <subtract_mode>" << std::endl;
        return 1;
    }

    // Read command line arguments
    const char* pressure_file = argv[1];
    float times_val = parseFractionOrFloat(argv[2]);
    float criteria_val = std::stof(argv[3]);
    std::string subtract_mode = argv[4];

    // Open the pressure.dat file
    std::ifstream pressure_in(pressure_file);
    if (!pressure_in) {
        std::cerr << "Error opening pressure.dat." << std::endl;
        return 1;
    }

    // Read all pressures from file
    std::vector<float> pressures;
    float value;
    while (pressure_in >> value) {
        pressures.push_back(value);
    }
    if (pressures.empty()) {
        std::cerr << "No data in pressure file." << std::endl;
        return 1;
    }

    // Determine the reference value to subtract based on subtract_mode
    float subtract_value = 0;
    if (subtract_mode == "last") {
        subtract_value = pressures.back();  // Use the last pressure value
    } else if (subtract_mode != "zero") {
        std::cerr << "Invalid subtract_mode. Use 'last' or 'zero'." << std::endl;
        return 1;
    }

    // Calculate results and output both the scaled values and the comparison results
    for (auto p : pressures) {
        float scaled_value = p * times_val;  // Calculate p times times_val
        float diff = std::abs(p - subtract_value) * times_val;  // Calculate absolute scaled difference
        bool result = diff < criteria_val;  // Compare scaled difference to criteria
        std::cout << scaled_value << " " << diff << " " << (result ? "true" : "false") << std::endl;
    }

    return 0;
}


EOF



cat > $L2_NORM_FILE << EOF
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>

int main(int argc, char* argv[]) {
    // Check if the file name is provided
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <force.dat>" << std::endl;
        return 1;
    }

    // Open the file
    std::ifstream file(argv[1]);
    if (!file) {
        std::cerr << "Error opening file." << std::endl;
        return 1;
    }

    std::string line;
    double x, y, z;

    // Read each line of the file
    while (getline(file, line)) {
        std::istringstream iss(line);

        // Read the three values from the line
        if (!(iss >> x >> y >> z)) {
            std::cerr << "Failed to read three values from line." << std::endl;
            continue;
        }

        // Calculate the L2 norm of the vector (x, y, z)
        double l2_norm = sqrt(x * x + y * y + z * z);

        // Output the L2 norm
        std::cout << l2_norm << std::endl;
    }

    file.close();
    return 0;
}

EOF



cat > $CONV_SUMMARY << EOF
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <vector>
#include <string>

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <encut.dat> <convergence.dat>" << std::endl;
        return 1;
    }

    const char* encut_file = argv[1];
    const char* convergence_file = argv[2];

    // Open the convergence.dat file
    std::ifstream convergence_in(convergence_file);
    if (!convergence_in) {
        std::cerr << "Error opening convergence.dat." << std::endl;
        return 1;
    }

    std::string line;
    int line_number = 0;
    int min_true_line = std::numeric_limits<int>::max();
    while (getline(convergence_in, line)) {
        std::istringstream iss(line);
        bool value, all_true = true;
        while (iss >> std::boolalpha >> value) {  // Reading with boolalpha to interpret "true" and "false" as boolean
            if (!value) {
                all_true = false;
                break;
            }
        }
        if (all_true && line_number < min_true_line) {
            min_true_line = line_number;  // Store the first all-true line number
        }
        line_number++;
    }
    convergence_in.close();

    if (min_true_line == std::numeric_limits<int>::max()) {
        std::cerr << "No all-true row found in convergence.dat." << std::endl;
        return 1;
    }

    // Open the encut.dat file
    std::ifstream encut_in(encut_file);
    if (!encut_in) {
        std::cerr << "Error opening encut.dat." << std::endl;
        return 1;
    }

    // Skip to the line in encut.dat that corresponds to min_true_line
    line_number = 0;
    std::string encut_value;
    while (getline(encut_in, encut_value)) {
        if (line_number == min_true_line) {
            std::cout << encut_value << std::endl;
            break;
        }
        line_number++;
    }

    encut_in.close();

    return 0;
}

EOF

echo "convergence.cpp, l2_norm.cpp, convergence_summary.cpp files have been successfully created."


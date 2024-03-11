#!/bin/bash

# INCARファイルのパスを指定
CPP_FILE="convergence_checker.cpp"

# Check if a file name is provided as the first argument
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <Threshold>"
  exit 1
fi

thres="$1"

# INCARファイルの内容を定義
cat > $CPP_FILE << EOF
#include <fstream>
#include <iostream>
#include <vector>

bool is_converged(const std::vector<double>& energies) {
  if (energies.size() <= 2) {
    return false;
  }

  bool converged_once = false;
  int consecutive_converged = 0;

  for (size_t i = 1; i < energies.size(); ++i) {
    double diff = (energies[i - 1] - energies[i]) / energies[i];

    if (diff >= -1*$thres) {
      consecutive_converged++;

      if (consecutive_converged == 2) {
        return true;
      }
    } else {
      consecutive_converged = 0;
    }
  }

  return false;
}

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <kp_e.data>" << std::endl;
    return 1;
  }

  std::ifstream infile(argv[1]);
  if (!infile.is_open()) {
    std::cerr << "Error opening file " << argv[1] << std::endl;
    return 1;
  }

  std::vector<double> energies;
  double energy;
  while (infile >> energy) {
    energies.push_back(energy);
  }

  bool converged = is_converged(energies);

  std::cout << (converged ? "True" : "False") << std::endl;

  return 0;
}

EOF

echo "convergence.cpp file has been successfully created."


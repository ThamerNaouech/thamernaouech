---
date: '2024-08-08T11:50:54.000Z'
title:  Optimization of Composite Laminate Structures for Enhanced Strength and Stability
tagline: 
preview: >-
  
image: >-
  https://github.com/ThamerNaouech/DL/blob/main/Structure-of-CFRP-laminates-15.png?raw=true
---

# Introduciton
Composite materials, especially Carbon Fiber Reinforced Polymers (CFRP), are extensively used in aerospace engineering due to their exceptional strength-to-weight ratio and flexibility in design. However, designing optimal stacking sequences for these materials can be challenging due to the need to meet stringent safety, stability, and material usage requirements. This project focuses on developing an algorithm to optimize the stacking sequence of composite laminates, ensuring they are structurally sound, efficient, and compliant with aerospace standards. By leveraging computational techniques, the project demonstrates how engineering problems can be solved systematically and effectively.


# Project Overview
The goal of this project is to create a C++ program that optimizes the stacking sequence of composite laminates by:

**1- Generating Possible Configurations:** The algorithm systematically produces all potential combinations of ply orientations {0°, 45°, −45°, 90°} for a given number of layers.

**2- Filtering for Balance and Symmetry:** Ensures that the stacking sequence is both balanced and symmetric to maintain stability and compliance with design principles.

**3- Performing Structural Calculations:** For each configuration, the algorithm calculates critical parameters such as:
* Transformed stiffness matrices (𝑄_bar).
* Critical buckling loads and shear responses.
* Reserve Factor (RF), ensuring a safety margin of at least 1.

**4- Iterative Optimization:** The program iteratively evaluates configurations, increasing the number of layers until the optimal design is achieved, minimizing material usage while meeting all requirements.

**5- Delivering Results:** The output is a stacking sequence that is lightweight, safe, and stable under given loads, along with key structural properties and reserve factors for verification.

This project represents a practical application of composite theory and computational engineering, providing a framework for efficient material design in high-performance industries like aerospace.

![Transfer Learning](https://github.com/ThamerNaouech/DL/blob/main/pann7.jpg?raw=true)



# Implementation


## 1. Understanding the Problem
Composite laminates, made up of layers of carbon fiber reinforced polymer (CFRP), are widely used in aerospace due to their strength, lightweight properties, and durability. The goal of this project is to:
* Design a stacking sequence that is balanced and symmetric.
* Ensure safety and strength under applied loads.
* Minimize material usage while meeting aerospace standards.




## 2. Define Key Concepts
Before coding, it’s important to understand the terms:
- **Ply orientations**: Layers of material can be oriented at specific angles: 
  {0°, 45°, −45°, 90°}.
- **Balanced and symmetric**: The stacking sequence should have equal positive and negative orientations and be symmetric around the midplane.
- **Reserve Factor (RF)**: A measure of the safety margin under applied loads (must be ≥ 1).
- **Buckling and shear loads**: Critical parameters to evaluate structural stability.






## 3. Input Parameters
The algorithm requires the following inputs:
- **Material properties**:  E_1, E_2, G_12, nu_12  (elastic modulus, shear modulus, and Poisson's ratio).
- **Applied loads**:  N_x, N_y, N_xy  (in-plane loads).
- **Laminate dimensions**:  a, b, t  (length, width, and ply thickness).




## 4. Plan the Algorithm
The algorithm consists of these steps:
1. **Generate possible ply orientations**: Create all combinations of {0°, 45°, −45°, 90°} up to a given number of layers.
2. **Filter for balanced and symmetric sequences**: Retain only those combinations that are symmetric around the midplane and balanced in orientation.
3. **Calculate transformed stiffness matrix Q_bar**: For each ply orientation, compute the transformed stiffness matrix using classical laminate theory.
4. **Evaluate reserve factor (RF)**: Determine the structural response under given loads and calculate the RF for each stacking sequence.
5. **Optimize for material usage**: Iterate through combinations to find the stacking sequence with the minimum number of layers that meets all criteria (RF, strength, and safety).








## 5. Implementation in C++



### Step 1: Setup Inputs
Create a structure for input parameters like material properties, loads, and laminate dimensions.
```jsx
      struct MaterialProperties {
        double E1, E2, G12, nu12;
      };
      struct Loads {
        double Nx, Ny, Nxy;
      };
      struct Dimensions {
        double a, b, t;
      };

```
 


 
### Step 2: Generate Ply Combinations
Write a recursive function to generate all possible stacking sequences.
```jsx
      void generateCombinations(std::vector<int>& combination, int depth, std::vector<std::vector<int>>&allCombinations) {
        if (depth == 0) {
          allCombinations.push_back(combination);
          return;
        }
        for (int angle : {0, 45, -45, 90}) {
          combination.push_back(angle);
          generateCombinations(combination, depth - 1, allCombinations);
          combination.pop_back();
        }
      }

```




### Step 3: Filter for Symmetry and Balance
Write a function to check if a sequence is symmetric and balanced.
```jsx
      bool isBalancedSymmetric(const std::vector<int>& sequence) {
        int n = sequence.size();
        for (int i = 0; i < n / 2; ++i) {
          if (sequence[i] != sequence[n - 1 - i]) return false;
        }
        // Additional checks for balance...
        return true;
      }

```



### Step 4: Compute Q_bar
Implement matrix calculations for transformed stiffness based on ply angles.
```jsx
      Matrix computeQBar(double angle, const MaterialProperties& matProps) {
        // Transform stiffness matrix for a given angle
        // (Matrix math omitted for brevity)
      }

```



### Step 5: Evaluate Reserve Factor
Evaluate the structural response under given loads and calculate RF.
```jsx
      double calculateReserveFactor(const std::vector<int>& stackingSequence, const MaterialProperties& matProps, const Loads& loads) {
        // Compute buckling, shear, and RF
        // (Details depend on composite theory)
      }

```

### Step 6: Iterate for Optimization
Run the algorithm iteratively, increasing the number of layers until the design meets all criteria.






## 6. Results
* **Output:** A stacking sequence that is balanced, symmetric, and meets safety standards with minimal material usage.
* **Visualization:** Optionally, you can plot the reserve factor against the number of layers or generate a diagram of the stacking sequence.




# Complete code
```jsx
    #include <iostream>
    #include <vector>
    #include <cmath>
    #include <iomanip>
    #include <limits>
    #include <unordered_map>

    // Constants
    const double PI = 3.14159265358979323846;
    double a, b, E1, E2, G12, t, v12; // Variables for dimensions and material properties
    const std::vector<int> possible_angles = { 0, 45, -45, 90 };  // Allowed ply orientations in degrees
    double Nx, Ny, Nxy; // External loads in x, y, and shear directions
    const double MIN_PERCENTAGE = 0.1;  // Minimum 10% for each used ply orientation

    // Q matrix elements
    double v21;
    double denom;
    double Q11, Q12, Q22, Q66;

    // Function to generate all combinations of ply orientations for a given number of plies
    void generate_combinations(int n, std::vector<int>& current_combination, int index, std::vector<std::vector<int>>& all_combinations) {
      if (index == n) {
        // Once a full combination is generated, add it to the list
        all_combinations.push_back(current_combination);
        return;
      }
      // Try each possible angle for the current ply and recurse
      for (int angle : possible_angles) {
        current_combination[index] = angle;
        generate_combinations(n, current_combination, index + 1, all_combinations);
      }
    }

    // Function to check if a combination meets the 10% rule for the used ply orientations
    bool meets_percentage_criteria(const std::vector<int>& combination) {
      std::unordered_map<int, int> ply_count; // count occurrences of each ply angle
      int total_plies = combination.size();

      // Count occurrences of each used ply orientation
      for (int angle : combination) {
        ply_count[angle]++;
      }

      // Check that each used orientation is at least 10% of the total plies
      for (std::unordered_map<int, int>::iterator it = ply_count.begin(); it != ply_count.end(); ++it) {
        if (it->second < total_plies * MIN_PERCENTAGE) {
          return false;
        }
      }

      return true;
    }

    // Function to find balanced and symmetric combinations
    std::vector<std::vector<int>> find_balanced_symmetric_combinations(int n) {
      int half_n = n / 2; // Half number of plies for symmetry
      bool is_odd = n % 2 != 0; // Check if odd number of plies
      std::vector<std::vector<int>> all_combinations;
      std::vector<int> current_combination(half_n);
      // Generate all combinations of half the plies number
      generate_combinations(half_n, current_combination, 0, all_combinations);

      std::vector<std::vector<int>> balanced_combinations; // Store the valid balanced and symmetric combinations
      for (const auto& combination : all_combinations) {
        int sum = 0;
        // Calculate the sum of angles for balance check
        for (int angle : combination) {
          sum += angle;
        }
        if (is_odd) {
          for (int middle_angle : possible_angles) {
            if (sum * 2 + middle_angle == 0) {
              std::vector<int> full_combination = combination;
              full_combination.push_back(middle_angle);
              full_combination.insert(full_combination.end(), combination.rbegin(), combination.rend());

              if (n < 10 || meets_percentage_criteria(full_combination)) {
                balanced_combinations.push_back(full_combination);
              }
            }
          }
        }
        else {
          if (sum == 0) {
            std::vector<int> full_combination = combination;
            full_combination.insert(full_combination.end(), combination.rbegin(), combination.rend());

            if (n < 10 || meets_percentage_criteria(full_combination)) {
              balanced_combinations.push_back(full_combination);
            }
          }
        }
      }
      return balanced_combinations;
    }

    // Function to calculate Q_bar matrix elements for a ply at a given orientation
    std::vector<std::vector<double>> calculate_Q_bar(double theta_rad) {
      std::vector<std::vector<double>> Q_bar(3, std::vector<double>(3, 0.0));

      double cos_theta = cos(theta_rad);
      double sin_theta = sin(theta_rad);
      double cos_2theta = cos_theta * cos_theta;
      double sin_2theta = sin_theta * sin_theta;
      double cos_4theta = cos_2theta * cos_2theta;
      double sin_4theta = sin_2theta * sin_2theta;

      // Transform the stiffness elements into the ply's local coordinate system
      double Q11_bar = Q11 * cos_4theta + 2 * (Q12 + 2 * Q66) * sin_2theta * cos_2theta + Q22 * sin_4theta;
      double Q12_bar = (Q11 + Q22 - 4 * Q66) * sin_2theta * cos_2theta + Q12 * (cos_4theta + sin_4theta);
      double Q22_bar = Q11 * sin_4theta + 2 * (Q12 + 2 * Q66) * sin_2theta * cos_2theta + Q22 * cos_4theta;
      double Q16_bar = (Q11 - Q12 - 2 * Q66) * sin_theta * cos_theta * cos_2theta + (Q12 - Q22 + 2 * Q66) * sin_theta * sin_theta * sin_theta * cos_theta;
      double Q26_bar = (Q11 - Q12 - 2 * Q66) * sin_theta * sin_2theta * cos_theta + (Q12 - Q22 + 2 * Q66) * cos_theta * cos_theta * cos_theta * sin_theta;
      double Q66_bar = (Q11 + Q22 - 2 * Q12 - 2 * Q66) * sin_2theta * cos_2theta + Q66 * (cos_4theta + sin_4theta);
      
      // Fill the Q_bar matrix
      Q_bar[0][0] = Q11_bar;
      Q_bar[0][1] = Q12_bar;
      Q_bar[0][2] = Q16_bar;
      Q_bar[1][0] = Q12_bar;
      Q_bar[1][1] = Q22_bar;
      Q_bar[1][2] = Q26_bar;
      Q_bar[2][0] = Q16_bar;
      Q_bar[2][1] = Q26_bar;
      Q_bar[2][2] = Q66_bar;

      return Q_bar;
    }

    // Function to find m, n
    void find_best_m_n(const std::vector<std::vector<double>>& D_matrix, double alpha, double beta, int& m, int& n) {
      double min_sig_crit = std::numeric_limits<double>::max(); // Initialize to maximum value
      int best_m = 1;
      int best_n = 1;

      // Iterate over possible m, n values (1 to 5)
      for (int current_m = 1; current_m <= 5; ++current_m) {
        for (int current_n = 1; current_n <= 5; ++current_n) {
          double m_alpha = current_m / alpha;
          double m_beta_n = beta * current_n * current_n;

          if (m_alpha * m_alpha + m_beta_n != 0) {  // Ensure no division by zero
            double sig_crit = (PI * PI / (b * b * t)) *
              (1 / (m_alpha * m_alpha + m_beta_n)) *
              (D_matrix[0][0] * pow(m_alpha, 4) +
                2 * (D_matrix[0][1] + D_matrix[2][2]) * pow(m_alpha * current_n, 2) +
                D_matrix[1][1] * pow(current_n, 4));

            // Check if this combination gives a lower sigma_crit value
            if (std::abs(sig_crit) < std::abs(min_sig_crit)) {
              min_sig_crit = sig_crit;
              best_m = current_m;
              best_n = current_n;
            }
          }
        }
      }

      m = best_m;
      n = best_n;
    }

    int main() {

      // Prompt user for inputs for material and dimensional properties
      std::cout << "Enter the value for a: ";
      std::cin >> a;

      std::cout << "Enter the value for b: ";
      std::cin >> b;

      std::cout << "Enter the value for t: ";
      std::cin >> t;

      std::cout << "Enter the value for E1: ";
      std::cin >> E1;

      std::cout << "Enter the value for E2: ";
      std::cin >> E2;

      std::cout << "Enter the value for G12: ";
      std::cin >> G12;

      std::cout << "Enter the value for v12: ";
      std::cin >> v12;

      // Calculate Q matrix elements based on inputs
      v21 = (v12 * E2) / E1;
      denom = 1 - v12 * v21;
      Q11 = E1 / denom;
      Q12 = (v12 * E2) / denom;
      Q22 = E2 / denom;
      Q66 = G12;

      // Prompt user for load values
      std::cout << "Enter the value for Nx: ";
      std::cin >> Nx;

      std::cout << "Enter the value for Ny: ";
      std::cin >> Ny;

      std::cout << "Enter the value for Nxy: ";
      std::cin >> Nxy;

      int num_stacks = 1; // Start with 1 ply and increase
      double RF_max = 0;
      std::vector<int> best_combination; // Store the best ply combination

      // Main loop to find the optimal combination of plies
      while (true) {

        double a_curr = a;
        double b_curr = b;
        // Find balanced and symmetric ply combinations for the current stack count
        auto balanced_symmetric_combinations = find_balanced_symmetric_combinations(num_stacks);
        
        // Variables for intermediate calculations
        double sig_x, sig_y, alpha, beta, tau_xy;

        // Determine if Nx or Ny is the dominant load, and adjust calculations accordingly
        if (std::abs(Nx) >= std::abs(Ny)) {
          sig_x = Nx / (t * num_stacks);
          sig_y = Ny / (t * num_stacks);
          alpha = a_curr / b_curr;
          beta = sig_y / sig_x;
          tau_xy = Nxy / (t * num_stacks);
        }
        else {
          sig_x = Ny / (t * num_stacks);
          sig_y = Nx / (t * num_stacks);
          std::swap(a_curr, b_curr); // Correctly swap a and b
          alpha = a_curr / b_curr;
          beta = sig_y / sig_x;
          tau_xy = Nxy / (t * num_stacks);
        }

        double current_max_RF = 0; // Track the max RF for current stack
        std::vector<int> current_best_combination; // Track the best combination for current stack

        for (const auto& combination : balanced_symmetric_combinations) {
          int n = combination.size();
          std::vector<double> z(n + 1);
          double total_thickness = n * t;
          z[0] = -total_thickness / 2;

          for (int i = 1; i <= n; ++i) {
            z[i] = z[i - 1] + t;
          }

          // Initialize A and D matrices for stiffness calculations
          std::vector<std::vector<double>> A_matrix(3, std::vector<double>(3, 0.0));
          std::vector<std::vector<double>> D_matrix(3, std::vector<double>(3, 0.0));

          for (int k = 0; k < n; ++k) {
            double theta_rad = combination[k] * PI / 180;  // Convert angle to radians
            auto Q_bar = calculate_Q_bar(theta_rad);  // Get the transformed stiffness matrix

            for (int i = 0; i < 3; ++i) {
              for (int j = 0; j < 3; ++j) {
                A_matrix[i][j] += Q_bar[i][j] * t;
                D_matrix[i][j] += Q_bar[i][j] * (pow(z[k + 1], 3) - pow(z[k], 3)) / 3.0;
              }
            }
          }

          double lambda = std::sqrt(D_matrix[0][0] * D_matrix[1][1]) / (D_matrix[0][1] + 2 * D_matrix[2][2]);
          double tau_crit;
          // Determine critical shear stress based on lambda
          if (lambda < 1) {
            tau_crit = (4 / (t * b_curr * b_curr)) *
              (std::sqrt(D_matrix[1][1] * (D_matrix[0][1] + 2 * D_matrix[2][2]))) *
              (11.7 + 0.532 * lambda + 0.938 * lambda * lambda);
          }
          else {
            tau_crit = (4 / (t * b_curr * b_curr)) *
              std::sqrt(std::sqrt(D_matrix[0][0] * std::pow(D_matrix[1][1], 3))) *
              (8.12 + 5.05 / lambda);
          }
          double RF_shear = std::abs(tau_crit / (1.5 * tau_xy));

          double RF_buckling, sig_crit, RF;
          if (beta == 0) {
            sig_crit = (12 * D_matrix[2][2] / (b_curr * b_curr * t)) +
              (std::sqrt(D_matrix[0][0] / D_matrix[1][1]) / (alpha * alpha * t));
            RF_buckling = std::abs(sig_crit / (1.5 * sig_x));
            RF = std::pow((1 / RF_buckling) + std::pow((1 / RF_shear), 2), -1);
          }
          else {
            int m = 0, n = 0;
            find_best_m_n(D_matrix, alpha, beta, m, n);

            sig_crit = (PI * PI / (b * b * t)) *
              (1 / ((m / alpha) * (m / alpha) + (n * n * beta))) *
              (D_matrix[0][0] * pow((m / alpha), 4) +
                2 * (D_matrix[0][1] + D_matrix[2][2]) * pow((m / alpha) * n, 2) +
                D_matrix[1][1] * pow(n, 4));
            RF_buckling = std::abs(sig_crit / (1.5 * sig_x));
            RF = std::pow((1 / RF_buckling) + std::pow((1 / RF_shear), 2), -1);
          }

          // Update the best combination based on the current maximum RF
          if (RF > current_max_RF) {
            current_max_RF = RF;
            current_best_combination = combination;
          }
        }
        // If the current maximum RF exceeds the previous maximum, update it
        if (current_max_RF > RF_max) {
          RF_max = current_max_RF;
          best_combination = current_best_combination;
        }
        // If RF_max is greater than or equal to 1, the best combination has been found
        if (RF_max >= 1) {
          break;
        }
        num_stacks++; // Increase the number of stacks and check again
      }
      
      // Output the final results
      std::cout << "RF_max: " << RF_max << std::endl;
      std::cout << "Best Combination: ";
      for (int angle : best_combination) {
        std::cout << angle << " ";
      }
      std::cout << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;

      return 0;
    }
```


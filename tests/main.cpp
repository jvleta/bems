#include <iostream>
#include <fstream>
#include <string>

#include "gtest/gtest.h"
#include "bem.h"

TEST(CELAP, Example1) {
  run_celap_ex1();
  std::ifstream current_output("../tests/celap/celap_example1.out");
  std::ifstream expected_output("../tests/celap/celap_example1.out_std");
  bool found_difference = false;
  std::string str1, str2;
  while (std::getline(current_output, str1) && std::getline(expected_output, str2)) {
    std::cout << str1 + " | " + str2 + "\n";
    if (str1 != str2) {
     found_difference = true;
   }
  }
  EXPECT_TRUE(!found_difference);
}

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}



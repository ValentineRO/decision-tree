#include "test.h"
#include <ctime>
using namespace std;

int main(){
  //std::srand(std::time(0));

  learning_test("exemple_papier",4);
  learning_test("iris",4);
  learning_test("wine",4);
  learning_test("blood_donation",4);
  learning_test("breast_cancer",4);

  return 0;
}

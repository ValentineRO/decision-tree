#include "test.h"
#include <vector>
#include <algorithm>
using namespace std;

int main(){
  
  try {
    string dataset_name = "";
    //test(dataset_name);
    //learning_test(dataset_name,4);
  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }

  return 0;
}

GUROBIDIR =/home/hurev/gurobi911/linux64/include
GUROBILIBDIR =/home/hurev/gurobi911/linux64/lib
CXXFLAGS = -g -Wall -Wextra
LD=g++

LDFLAGS=-g -m64 -fPIC -D_GLIBC_DEBUG -DNOSHORTS -fopenmp 
# -fpermissive
EXEC=decision_tree

all:$(EXEC)

$(EXEC): main.o test.o train_tree.o decontracting_algorithm.o models.o solution.o model_utilities.o clusters.o objective.o following_constraints.o tree_structuring_constraints.o branching_constraints.o counting_errors.o tree.o data_formatting.o
	$(LD) main.cpp $(LDFLAGS)  -I$(GUROBIDIR) -o  $(EXEC) test.o train_tree.o decontracting_algorithm.o models.o solution.o model_utilities.o clusters.o objective.o following_constraints.o tree_structuring_constraints.o branching_constraints.o counting_errors.o tree.o data_formatting.o -ldl -lpthread -L$(LIB64DIR) -llapack -lblas -lgfortran -lm -L$(GUROBILIBDIR) -lgurobi_g++5.2 -lgurobi91

main.o : main.cpp test.h
	g++ -c -g main.cpp -I$(GUROBIDIR)

test.o : test.cpp test.h train_tree.h $(GUROBIDIR)/gurobi_c++.h
	g++ -c -g test.cpp -I$(GUROBIDIR)

train_tree.o : train_tree.cpp train_tree.h decontracting_algorithm.h $(GUROBIDIR)/gurobi_c++.h
	g++ -c -g train_tree.cpp -I$(GUROBIDIR)

decontracting_algorithm.o : decontracting_algorithm.cpp decontracting_algorithm.h models.h $(GUROBIDIR)/gurobi_c++.h
	g++ -c -g decontracting_algorithm.cpp  -I$(GUROBIDIR)

models.o : models.cpp models.h model_utilities.h objective.h following_constraints.h branching_constraints.h tree_structuring_constraints.h counting_errors.h clusters.h
	g++ -c -g models.cpp -I$(GUROBIDIR)

clusters.o : clusters.cpp clusters.h model_utilities.h $(GUROBIDIR)/gurobi_c++.h
	g++ -c -g clusters.cpp -I$(GUROBIDIR)

objective.o : objective.cpp objective.h $(GUROBIDIR)/gurobi_c++.h
	g++ -c -g objective.cpp -I$(GUROBIDIR)

following_constraints.o : following_constraints.cpp following_constraints.h  $(GUROBIDIR)/gurobi_c++.h
	g++ -c -g following_constraints.cpp -I$(GUROBIDIR)

tree_structuring_constraints.o : tree_structuring_constraints.cpp tree_structuring_constraints.h $(GUROBIDIR)/gurobi_c++.h
	g++ -c -g tree_structuring_constraints.cpp -I$(GUROBIDIR)

branching_constraints.o : branching_constraints.cpp branching_constraints.h $(GUROBIDIR)/gurobi_c++.h
	g++ -c -g branching_constraints.cpp -I$(GUROBIDIR)

counting_errors.o : counting_errors.cpp counting_errors.h $(GUROBIDIR)/gurobi_c++.h
	g++ -c -g counting_errors.cpp -I$(GUROBIDIR)

model_utilities.o : model_utilities.cpp model_utilities.h solution.h $(GUROBIDIR)/gurobi_c++.h
	g++ -c -g model_utilities.cpp -I$(GUROBIDIR)

solution.o: solution.cpp solution.h tree.h $(GUROBIDIR)/gurobi_c++.h
	g++ -c -g solution.cpp -I$(GUROBIDIR)

tree.o: tree.cpp tree.h data_formatting.h $(GUROBIDIR)/gurobi_c++.h
	g++ -c -g tree.cpp -I$(GUROBIDIR)

data_formatting.o: data_formatting.cpp data_formatting.h
	g++ -c -g data_formatting.cpp 

clean :
	rm *.o

run: decision_tree
	./decision_tree

all: Compile Link Execute

Compile:
		g++ -c -std=c++11 main_project3.cpp -larmadillo

Link:
		g++ -std=c++11 main_project3.o -o main_project3.exe  -larmadillo

Execute:
		./main_project3.exe

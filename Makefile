all: Compile Link Execute


Compile:
		g++ -c -std=c++11 ex10project3.cpp -larmadillo

Link:
		g++ -std=c++11 ex10project3.o -o ex10project3.exe  -larmadillo

Execute:
		./ex10project3.exe

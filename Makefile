10: Compile_10 Link_10 Execute_10
9 : Compile_9 Link_9 Execute_9


Compile_10:
		g++ -c -std=c++11 ex10project3.cpp -larmadillo

Link_10:
		g++ -std=c++11 ex10project3.o -o ex10project3.exe  -larmadillo

Execute_10:
		./ex10project3.exe

Compile_9:
		g++ -c -std=c++11 main_project3.cpp -larmadillo

Link_9:
		g++ -std=c++11 main_project3.o -o main_project3.exe  -larmadillo

Execute_9:
		./main_project3.exe

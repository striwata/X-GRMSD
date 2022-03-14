CC = g++ -O3 -g
CFLAGS = -std=c++11

test: test.o main_function.o hung.o local_ICP.o data_form.o basic_function.o cube.o
	$(CC) -o test.exe test.o main_function.o hung.o local_ICP.o data_form.o basic_function.o cube.o 

data_form.o: data_form.cpp data_form.h
	$(CC) -c data_form.cpp -o data_form.o

basic_function.o: basic_function.cpp basic_function.h
	$(CC) -c basic_function.cpp -o basic_function.o

hung.o: Hungarian.cpp Hungarian.h
	$(CC) -c Hungarian.cpp -o hung.o

cube.o: cube.cpp cube.h
	$(CC) $(CFLAGS) -c cube.cpp -o cube.o

local_ICP.o: local_ICP.cpp local_ICP.h data_form.h basic_function.h
	$(CC) -c local_ICP.cpp -o local_ICP.o

main_function.o: main_function.cpp main_function.h local_ICP.h data_form.h basic_function.h cube.h
	$(CC) $(CFLAGS) -c main_function.cpp -o main_function.o

test.o: test.cpp main_function.h
	$(CC) -c test.cpp -o test.o 

all:
	g++ -o test main.cpp -std=c++17

run: all
	./test

mesh:
	g++ -o test test.cpp -std=c++17

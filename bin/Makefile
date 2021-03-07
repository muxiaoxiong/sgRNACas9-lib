#Off-Spotter make file
all: Results Load_Memory Table_Creation Detach_Memory

Results: Results.cpp
	g++ Results.cpp -m64 -O3 -o Results

Load_Memory: Load_Memory.cpp
	g++ Load_Memory.cpp -m64 -O3 -o Load_Memory

Table_Creation: Table_Creation.cpp
	g++ Table_Creation.cpp -m64 -O3 -o Table_Creation

Detach_Memory: Detach_Memory.cpp
	g++ Detach_Memory.cpp -m64 -O3 -o Detach_Memory

clean:
	rm Detach_Memory Table_Creation Load_Memory Results

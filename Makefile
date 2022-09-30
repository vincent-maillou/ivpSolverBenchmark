all: RK4_LargeSystem.out 

RK4_LargeSystem.out: RK4_LargeSystem.cc
	g++ -std=c++17 -I ../eigen  RK4_LargeSystem.cc -o RK4_LargeSystem.out 

clean:
	rm *.out

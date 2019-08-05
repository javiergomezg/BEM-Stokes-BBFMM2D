 cc = g++-6
OBJSB = H2_2D_Tree.o H2_2D_Node.o kernel_Base.o kernel_Types.o

# Con esta regla se compila el programa Stokes
all: $(OBJSB)

H2_2D_Tree.o: H2_2D_Tree.cpp header/H2_2D_Tree.hpp
	$(cc) -c H2_2D_Tree.cpp

H2_2D_Node.o: H2_2D_Node.cpp header/H2_2D_Node.hpp
	$(cc) -c H2_2D_Node.cpp

kernel_Base.o: kernel_Base.cpp header/kernel_Base.hpp
	$(cc) -c kernel_Base.cpp

kernel_Types.o: kernel_Types.cpp header/kernel_Types.hpp
	$(cc) -c kernel_Types.cpp

clean:
	rm -f $(OBJSB) 

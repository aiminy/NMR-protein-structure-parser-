all:main_read_pdb main_write_title
objects1 =atom.o main_read_pdb.o aa.o protein.o chain.o Vector3D.o aa_summary.o model.o 
objects2 =main_write_title.o
main_read_pdb:$(objects1)
	g++ -o main_read_pdb $(objects1)
main_write_title:$(objects2)
	g++ -o main_write_title $(objects2)
main_read_pdb.o:main_read_pdb.cpp linked_list.h
	g++ -c main_read_pdb.cpp
main_write_title.o:main_write_title.cpp
	g++ -c main_write_title.cpp
atom.o:atom.cpp atom.h
	g++ -c atom.cpp
Vector3D.o:Vector3D.h
	g++ -c Vector3D.cpp
aa.o:aa.cpp aa.h linked_list.h 
	g++ -c aa.cpp
chain.o:chain.cpp chain.h aa.h linked_list.h
	g++ -c chain.cpp
model.o:model.cpp model.h aa.h linked_list.h
	g++ -c model.cpp
aa_summary.o:aa_summary.cpp aa_summary.h
	g++ -c aa_summary.cpp
#protein.o:protein.cpp protein.h aa_summary.h
#	g++ -c protein.cpp
protein.o:protein.cpp protein.h aa_summary.h
	R CMD SHLIB protein.cpp
clean:
	rm -f main_read_pdb *.o

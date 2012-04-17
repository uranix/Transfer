#include <mesh.h>
#include <iostream>
#include <windows.h>

int main() {
	Mesh m("D:\\Projects\\NetGen\\radiation.vol");
	std::cout << "Mesh " << (m.check()?"OK":"Bad!") << std::endl;
	m.saveVtk("D:\\Projects\\NetGen\\radiation.vtk");
	system("pause");
	return 0;
}
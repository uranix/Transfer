#include "Config.h"

#include <stdio.h>
#include <string.h>

#define die(x, args...) do { fprintf(stderr, "%s:%d: " x, cfgfile, line , ##args); throw 0; } while (0)

Config::Config(const char *cfgfile) {
	FILE *f = fopen(cfgfile, "r");
	char buf[1024];
	if (!f) {
		fprintf(stderr, "failed to open file `%s'\n", cfgfile);
		throw 0;
	}

	nreg = -1;
	int line = 0;
	int rid = -1;
	dev = 0;
	dump = 0;
	strcpy(outfn, "solution.vtk");

	while (fgets(buf, 1024, f)) {
		line ++;
		if (buf[0] == '#')
			continue;
		if (buf[0] == '\n')
			continue;
		if (buf[0] == '[') {
			if (sscanf(buf, "[region %d]", &rid) != 1)
				die("scanf failed to parse region section\n");
			if (nreg < 0)
				die("section appeared before numregs was set\n");
			if (rid > nreg || rid < 0)
				die("region id (%d) is out of bounds (0-%d)\n", rid, nreg);
			continue;
		}
		if (!strncmp(buf, "kappa", 5)) {
			if (rid < 0) 
				die("kappa parameter outside section");
			if (sscanf(buf, "kappa = %lf", kappa + rid) != 1) 
				die("scanf failed to parse kappa parameter\n");
			continue;
		}
		if (!strncmp(buf, "Ip", 2)) {
			if (rid < 0) 
				die("Ip parameter outside section");
			if (sscanf(buf, "Ip = %lf", Ip + rid) != 1) 
				die("scanf failed to parse Ip parameter\n");
			continue;
		}
		if (!strncmp(buf, "meshfile = ", 11)) {
			if (rid >= 0) 
				die("meshfile parameter in section");
			buf[strlen(buf)-1] = 0;
			strncpy(meshfn, buf + 11, 1023);
			meshfn[1023] = 0;
			continue;
		}
		if (!strncmp(buf, "outfile = ", 10)) {
			if (rid >= 0) 
				die("outfile parameter in section");
			buf[strlen(buf)-1] = 0;
			strncpy(outfn, buf + 10, 1023);
			outfn[1023] = 0;
			continue;
		}
		if (!strncmp(buf, "numregs", 7)) {
			if (rid >= 0) 
				die("numregs parameter in section");
			if (sscanf(buf, "numregs = %d", &nreg) != 1) 
				die("scanf failed to parse numregs parameter\n");
			Ip = new double[nreg + 1];
			kappa = new double[nreg + 1];
			continue;
		}
		if (!strncmp(buf, "angorder", 8)) {
			if (rid >= 0) 
				die("angorder parameter in section");
			if (sscanf(buf, "angorder = %d", &maxk) != 1)
				die("scanf failed to parse angorder parameter\n");
			continue;
		}
		if (!strncmp(buf, "device", 6)) {
			if (rid >= 0) 
				die("device parameter in section");
			if (sscanf(buf, "device = %d", &dev) != 1) 
				die("scanf failed to parse device parameter\n");
			continue;
		}
		if (!strncmp(buf, "dump", 4)) {
			if (rid >= 0) 
				die("dump parameter in section");
			if (sscanf(buf, "dump = %d", &dump) != 1) 
				die("scanf failed to parse dump parameter\n");
			continue;
		}
	}

	fclose(f);
	printf("Config parameters:\n");
	printf("\tnregions = %d\n", nreg);
	printf("\tangorder = %d\n", maxk);
	printf("\tdevice = %d\n", dev);
	printf("\tdump = %s\n", dump?"true":"false");
	printf("\tkappa = { ");
	for (int i = 1; i < nreg; i++)
		printf("%2.6e, ", kappa[i]);
	printf("%2.6e}\n", kappa[nreg]);
	printf("\tIp = { ");
	for (int i = 1; i < nreg; i++)
		printf("%2.6e, ", Ip[i]);
	printf("%2.6e}\n", Ip[nreg]);
	printf("\tMesh file : %s\n", meshfn);
	printf("\tSolution file : %s\n", outfn);

}

#undef die

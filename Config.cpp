#include "Config.h"

#include <stdio.h>
#include <string.h>

#define die(x, args...) fprintf(stderr, "%s:%d: " x, cfgfile, line , ##args)

Config::Config(const char *cfgfile) {
	FILE *f = fopen(cfgfile, "r");
	char buf[1024];
	if (!f)
		throw 123;

	nreg = -1;
	int line = -1;
	int rid = -1;

	while (fgets(buf, 1024, f)) {
		line ++;
		if (buf[0] == '#')
			continue;
		if (buf[0] == '\n')
			continue;
		if (buf[0] == '[') {
			if (sscanf(buf, "[region %d]", &rid) != 1) {
				die("scanf failed to parse section\n");
				throw 123;
			}
			if (nreg < 0) {
				die("section appeared before numregs was set\n");
				throw 124;
			}
			if (rid > nreg || rid < 0) {
				die("region id (%d) is out of bounds (0-%d)\n", rid, nreg);
				throw 125;
			}
			continue;
		}
		if (!strncmp(buf, "kappa", 5)) {
			if (rid < 0) {
				die("kappa parameter outside section");
				throw 126;
			}
			if (sscanf(buf, "kappa = %lf", kappa + rid) != 1) {
				die("scanf failed to parse parameter\n");
				throw 127;
			}
			continue;
		}
		if (!strncmp(buf, "Ip", 2)) {
			if (rid < 0) {
				die("Ip parameter outside section");
				throw 128;
			}
			if (sscanf(buf, "Ip = %lf", Ip + rid) != 1) {
				die("scanf failed to parse parameter\n");
				throw 129;
			}
			continue;
		}
		if (!strncmp(buf, "meshfile = ", 11)) {
			if (rid >= 0) {
				die("meshfile parameter in section");
				throw 130;
			}
			strncpy(meshfn, buf + 11, 1023);
			meshfn[1023] = 0;
			continue;
		}
		if (!strncmp(buf, "numregs", 7)) {
			if (rid >= 0) {
				die("numregs parameter in section");
				throw 131;
			}
			if (sscanf(buf, "numregs = %d", &nreg) != 1) {
				die("scanf failed to parse parameter\n");
				throw 132;
			}
			continue;
		}
		if (!strncmp(buf, "angorder", 8)) {
			if (rid >= 0) {
				die("angorder parameter in section");
				throw 133;
			}
			if (sscanf(buf, "angorder = %d", &maxk) != 1) {
				die("scanf failed to parse parameter\n");
				throw 134;
			}
			continue;
		}
		if (!strncmp(buf, "device", 6)) {
			if (rid >= 0) {
				die("device parameter in section");
				throw 135;
			}
			if (sscanf(buf, "device = %d", &maxk) != 1) {
				die("scanf failed to parse parameter\n");
				throw 136;
			}
			continue;
		}
	}

	fclose(f);
}

#undef die

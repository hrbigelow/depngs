#ifndef _DEP_H
#define _DEP_H

int main_comp(int argc, char ** argv);

int main_mode(int argc, char ** argv);

int main_dist(int argc, char ** argv);

int main_simp(int argc, char ** argv);

int main_anomaly(int argc, char ** argv);

int main_simc(int argc, char ** argv);

int main_bqs(int argc, char ** argv);

int main_bqslocus(int argc, char ** argv);

int main_bqs2jpd(int argc, char ** argv);

extern "C" {
    int main_pug(int argc, char ** argv);
    int main_diststats(int argc, char **argv);
}

#endif // _DEP_H

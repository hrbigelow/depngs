#include "kbtree.h"
#include <stdio.h>

int cmp(const int *a, const int *b)
{
    if (*a < *b) return -1;
    if (*b < *a) return 1;
    return 0;
}

KBTREE_INIT(itree, int *, cmp)

int main(int argc, char **argv)
{
    int q = strtod(argv[1], NULL);

    int p[] = { 10, 20, 30, 30, 30, 40, 50, 60 };
    unsigned n_p = sizeof(p) / sizeof(p[0]);

    kbtree_t(itree) *t;
    t = kb_init(itree, 128);

    unsigned i;
    for (i = 0; i != n_p; ++i)
        kb_put(itree, t, p + i);

    fprintf(stdout, "kb_size = %u\n", kb_size(t));
    int **l, **u;
    kb_interval(itree, t, &q, &l, &u);
    fprintf(stdout, "l = %i (%i), u = %i (%i)\n", 
            (*l && **l ? **l : -1), 
            (*l ? *l - p : -1), 
            (*u && **u ? **u : -1), 
            (*u ? *u - p : -1));

    int **eq;
    eq = kb_get(itree, t, &q);
    fprintf(stdout, "eq = %i\n", (eq && *eq) ? **eq : -1);

    kb_destroy(itree, t);
    return 0;
}

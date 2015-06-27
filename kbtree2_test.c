#include "kbtree.h"
#include <stdio.h>

int cmp(int a, int b)
{
    if (a < b) return -1;
    if (b < a) return 1;
    return 0;
}

KBTREE_INIT(itree, int, cmp)

int main(int argc, char **argv)
{
    int q = strtod(argv[1], NULL);

    int p[] = { 10, 20, 30, 30, 30, 40, 50, 60 };
    unsigned n_p = sizeof(p) / sizeof(p[0]);

    kbtree_t(itree) *t;
    t = kb_init(itree, 128);

    unsigned i;
    for (i = 0; i != n_p; ++i)
        kb_put(itree, t, p[i]);

    int *l, *u;
    kb_interval(itree, t, q, &l, &u);
    fprintf(stdout, "l = %i, u = %i\n", *l, *u);
    kb_destroy(itree, t);
    return 0;
}

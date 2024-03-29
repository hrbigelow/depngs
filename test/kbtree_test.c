#include "kbtree.h"
#include <stdio.h>


typedef struct {
    uint64_t u,v;
} pair_t;

/* intervals are sorted in the usual way.  those which have a full
   containment relationship are deemed equal.*/
int cmp_ival(pair_t a, pair_t b)
{
    if (a.v <= b.u) return -1;
    if (b.v <= a.u) return 1;
    if (a.u < b.u && a.v < b.v) return -1;
    if (b.u < a.u && b.v < a.v) return 1;
    return 0;
}

KBTREE_INIT(ival, pair_t, cmp_ival)

int main(int argc, char **argv)
{
    /* pair_t q; */
    /* q.u = strtod(argv[1], NULL); */
    /* q.v = strtod(argv[2], NULL); */

    pair_t p[] = { { 10, 20 }, { 30, 40 }, { 50, 60 }, { 70, 71 } };
    unsigned n_p = sizeof(p) / sizeof(p[0]);

    kbtree_t(ival) *t;
    t = kb_init(ival, 128);

    unsigned i;
    for (i = n_p; i != 0; --i)
        kb_put(ival, t, p[i - 1]);

    kbitr_t itr;
    kb_itr_first_ival(t, &itr);
    pair_t pn = kb_itr_key(pair_t, &itr);
    printf("{ %Zu, %Zu }\n", pn.u, pn.v);

    int ret = 1;
    while (ret) {
        pn = kb_itr_key(pair_t, &itr);
        printf("{ %Zu, %Zu }\n", pn.u, pn.v);
        ret = kb_itr_next_ival(t, &itr);
    }

    /* pair_t *l, *u; */
    /* kb_interval(ival, t, q, &l, &u); */
    /* l ? printf("l = { %Zu, %Zu }, ", l->u, l->v) : printf("l = <NULL>, "); */
    /* u ? printf("u = { %Zu, %Zu }", u->u, u->v) : printf("u = <NULL>"); */
    /* printf("\n"); */
    kb_destroy(ival, t);
    return 0;
}

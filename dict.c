/* functions for searching a dictionary of strings mapped to
   integers */

#include "cache.h"
#include <stddef.h>
#include <string.h>
#include <stdlib.h>

struct dict {
    const char *key;
    unsigned long index;
};

static struct dict *index_buf;
static char *key_buf = NULL; /* holds zero-terminated keys that each dict->key
                                points into */

static size_t index_nr = 0;
static size_t index_alloc = 0;
static size_t num_key_chars = 0;
static size_t num_chars_alloc = 0;

/* internal sorting function */
int dict_less_key(const void *pa, const void *pb)
{
    const struct dict
        *a = (struct dict *)pa,
        *b = (struct dict *)pb;

    return strcmp(a->key, b->key);
}




/* add item to the dictionary. does not check for duplicate keys or
   duplicate index values. */
void dict_add_item(const char* key, unsigned long index)
{
    size_t keysize = strlen(key) + 1;
    ALLOC_GROW(index_buf, index_nr + 1, index_alloc);
    ALLOC_GROW(key_buf, num_key_chars + keysize, num_chars_alloc);

    struct dict item;
    strncpy(key_buf + num_key_chars, key, keysize);
    item.key = key_buf + num_key_chars;
    num_key_chars += keysize;

    item.index = index;
    index_buf[index_nr++] = item;
}

/* call this after all items added, to get the index ready for
   searches */
void dict_build()
{
    qsort(index_buf, index_nr, sizeof(index_buf[0]), dict_less_key);
}


/* release static resources.  after calling this, it is okay to add
   items / build again*/
void dict_free()
{
    free(key_buf);
    key_buf = NULL;
    free(index_buf);
    index_buf = NULL;
    index_nr = 0;
    index_alloc = 0;
    num_key_chars = 0;
    num_chars_alloc = 0;
}


/* return index value of key, or -1 if not found */
long dict_search(const char *key)
{
    struct dict q, *res;
    q.key = key; /* dict does not own the string */
    res = (struct dict *)
        bsearch(&q, index_buf, index_nr, sizeof(struct dict), dict_less_key);
    return res ? (long)res->index : -1;
}

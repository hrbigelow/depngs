/* add item to the dictionary. does not check for duplicate keys or
   duplicate index values. */
void dict_add_item(const char* key, unsigned long index);

/* call this after all items added, to get the index ready for
   searches */
void dict_build();

/* release static resources.  after calling this, it is okay to add
   items / build again*/
void dict_free();

/* return index value of key, or -1 if not found */
long dict_search(const char *key);

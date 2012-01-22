#ifndef SET_H
#define SET_H
#include <glib.h>

#define set_new() g_hash_table_new_full(g_direct_hash, g_direct_equal, NULL, NULL)
#define set_insert(set, element) g_hash_table_insert(set, element, element)
#define set_include(set, element) g_hash_table_lookup_extended(set, element, NULL, NULL)
#define set_remove(set, element) g_hash_table_remove(set, element)
#define set_size(set) g_hash_table_size(set)
#define set_destroy(set) g_hash_table_destroy(set)

#define set_iter_init(iter, set) g_hash_table_iter_init(&iter, set);
#define set_iter_next(iter, key, value) g_hash_table_iter_next(&iter, &key, &value)

typedef GHashTable set;
typedef GHashTableIter set_iter;

#endif

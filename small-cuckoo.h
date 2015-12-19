/** -*- mode: C; c-file-style: "k&r" -*-
 * Implements memory-mappable Cuckoo hash table for less than 64k keys.
 * @file small-cuckoo.h
 */
#pragma once

#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct small_cuckoo {
     size_t table_size;
     uint16_t *table;
     uint16_t n_entries, entries_len;
     struct {
          uint64_t key;
          uint64_t value;
     } *entries;
} small_cuckoo;

typedef struct small_cuckoo_iter {
     small_cuckoo *sc;
     uint16_t i;
} small_cuckoo_iter;

extern small_cuckoo small_cuckoo_new(size_t initial_size);
extern void small_cuckoo_insert(small_cuckoo *sc, uint64_t key, uint64_t value);
extern bool small_cuckoo_find(small_cuckoo *sc, uint64_t key, uint64_t *value);
extern void small_cuckoo_free(small_cuckoo *sc);
extern void small_cuckoo_serialize(int fd, small_cuckoo *sc);
extern void small_cuckoo_deserialize(int fd, small_cuckoo *sc);

extern void small_cuckoo_iterate(small_cuckoo *sc, small_cuckoo_iter *iter);
extern bool small_cuckoo_iter_has_next(small_cuckoo_iter *iter);
extern void small_cuckoo_iter_next(small_cuckoo_iter *iter, uint64_t *key, uint64_t *value);




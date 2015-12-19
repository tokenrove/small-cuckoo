/** -*- mode: C; c-file-style: "k&r" -*-
 * An implementation of Cuckoo hashing for small tables.
 *
 * @see Pagh, Rasmus; Rodler, Flemming Friche (2001). "Cuckoo
 * Hashing". Algorithms — ESA 2001. Lecture Notes in Computer
 * Science. 2161. pp. 121–133. doi:10.1007/3-540-44676-1_10. ISBN
 * 978-3-540-42493-2.
 */

#include "small-cuckoo.h"
#include "ensure.h"
#include "bithacks.h"

/* Larson's hash function.
 * Described in Per-Ake Larson, Dynamic Hash Tables, CACM 31(4), April 1988, pp. 446--457.
 * Acceptable according to <http://www.strchr.com/hash_functions>. */
uint16_t larsons_hash(uint64_t key)
{
     uint32_t h = 0xdeadbeef;
     uint8_t *s = (uint8_t *)&key;
     enum { M = 101 };
     h = h * M + *s++;
     h = h * M + *s++;
     h = h * M + *s++;
     h = h * M + *s++;
     h = h * M + *s++;
     h = h * M + *s++;
     h = h * M + *s++;
     return h ^ (h>>16);
}

static uint16_t hash_1(size_t n, uint64_t key)
{
     return (larsons_hash(key) & ((n>>1)-1))<<1;
}

/* Use CRC32 if we have it in hardware, Bob Jenkins's stuff otherwise.
 * Acceptable according to <http://www.strchr.com/hash_functions>. */
#ifdef __SSE4_2__

#include <x86intrin.h>

static uint16_t hash_2(size_t n, uint64_t key)
{
     uint32_t h;
#ifdef __x86_64__
     h = _mm_crc32_u64(-1, key);
#else
     h = _mm_crc32_u32(-1, ((uint32_t*)&key)[0]);
     h = _mm_crc32_u32(h, ((uint32_t*)&key)[1]);
#endif
     h ^= (h>>16);
     return 1 + ((h & ((n>>1)-1))<<1);
}

#else

/* Hash function due to Bob Jenkins (original code in the public
 * domain).  I have removed most of Bob's comments from this code, but
 * lots of information about it can be found at
 * <http://burtleburtle.net/bob/hash>. */
#define hashsize(n) ((uint32_t)1<<(n))
#define hashmask(n) (hashsize(n)-1)
#define rot(x,k) (((x)<<(k)) | ((x)>>(32-(k))))

/* mix 3 32-bit values reversibly. */
#define mix(a,b,c)                              \
     {                                          \
          a -= c;  a ^= rot(c, 4);  c += b;     \
          b -= a;  b ^= rot(a, 6);  a += c;     \
          c -= b;  c ^= rot(b, 8);  b += a;     \
          a -= c;  a ^= rot(c,16);  c += b;     \
          b -= a;  b ^= rot(a,19);  a += c;     \
          c -= b;  c ^= rot(b, 4);  b += a;     \
     }

/* final mixing of 3 32-bit values (a,b,c) into c */
#define final(a,b,c)                            \
     {                                          \
          c ^= b; c -= rot(b,14);               \
          a ^= c; a -= rot(c,11);               \
          b ^= a; b -= rot(a,25);               \
          c ^= b; c -= rot(b,16);               \
          a ^= c; a -= rot(c,4);                \
          b ^= a; b -= rot(a,14);               \
          c ^= b; c -= rot(b,24);               \
     }

static uint32_t hashword(
const uint32_t *k,                   /* the key, an array of uint32_t values */
size_t          length,               /* the length of the key, in uint32_ts */
uint32_t        initval)         /* the previous hash, or an arbitrary value */
{
     uint32_t a,b,c;

     /* Set up the internal state */
     a = b = c = 0xdeadbeef + (((uint32_t)length)<<2) + initval;

     for (; length > 3; length -= 3, k += 3) {
          a += k[0]; b += k[1]; c += k[2];
          mix(a,b,c);
     }

     switch(length) {
     case 3: c+=k[2];
     case 2: b+=k[1];
     case 1: a+=k[0];
          final(a,b,c);
     }

     return c;
}

static uint16_t hash_2(size_t n, uint64_t key)
{
     uint32_t h;
     h = hashword((const uint32_t *)&key, 2, 0x55555555);
     h ^= (h>>16);
     return 1 + ((h & ((n>>1)-1))<<1);
}

#endif


small_cuckoo small_cuckoo_new(size_t initial_size)
{
     small_cuckoo sc = {0};
     sc.table_size = 1<<(ceil_pow2(initial_size)+1);
     ENSURE(sc.table = calloc(sc.table_size, sizeof sc.table[0]));
     sc.n_entries = 1;          /* Entry 0 is special. */
     sc.entries_len = 1+initial_size;
     ENSURE(sc.entries = malloc(sc.entries_len * sizeof sc.entries[0]));
     return sc;
}

static void insert(small_cuckoo *sc, uint16_t i);

static void double_size(small_cuckoo *sc)
{
     uint16_t *prev_table = sc->table;
     sc->table_size <<= 1;
     ENSURE(sc->table = calloc(sc->table_size, sizeof sc->table[0]));
     for (unsigned i = 0; i < sc->table_size>>1; ++i) {
          uint16_t k = prev_table[i];
          if (k) insert(sc, k);
     }
     free(prev_table);
}

enum { MAX_LOOPS = 20 };

static void insert(small_cuckoo *sc, uint16_t i)
{
     uint16_t h;
     for (size_t n = MAX_LOOPS; n > 0; --n) {
#define X(fn)                                                  \
          h = fn(sc->table_size, sc->entries[i].key);          \
          i ^= sc->table[h];                                   \
          sc->table[h] ^= i;                                   \
          i ^= sc->table[h];                                   \
          if (i == 0) return;

          X(hash_1);
          X(hash_2);
#undef X
     }

     double_size(sc);
     insert(sc, i);
}

void small_cuckoo_insert(small_cuckoo *sc, uint64_t key, uint64_t value)
{
     uint16_t i = sc->n_entries;
     ENSURE(i > 0);
     ++sc->n_entries;
     if (sc->n_entries >= sc->entries_len) {
          sc->entries_len <<= 1;
          ENSURE(sc->entries = realloc(sc->entries, sc->entries_len * sizeof sc->entries[0]));
     }
     sc->entries[i].key = key;
     sc->entries[i].value = value;
     insert(sc, i);
}

bool small_cuckoo_find(small_cuckoo *sc, uint64_t key, uint64_t *value)
{
     uint16_t i;
#define X(h)                                            \
     i = sc->table[h];                                  \
     if (i && sc->entries[i].key == key) {              \
          if (value) *value= sc->entries[i].value;      \
          return true;                                  \
     }
     X(hash_1(sc->table_size, key));
     X(hash_2(sc->table_size, key));
     return false;
#undef X
}

void small_cuckoo_free(small_cuckoo *sc)
{
     if (sc->table) free(sc->table);
     if (sc->entries) free(sc->entries);
     *sc = (small_cuckoo){0};
}

/* We only write out the entries, not the table; it gets reconstructed
 * when we read the metadata.
 */
void small_cuckoo_serialize(int fd, small_cuckoo *sc)
{
#define WRITE_UNDER(t,x,n) do { uint32_t u = t(x); ENSURE(n == write(fd, &u, n)); } while(0)
     WRITE_UNDER(htole16, sc->n_entries, 2);
     for (uint16_t i = 0; i < sc->n_entries; ++i) {
          WRITE_UNDER(htole64, sc->entries[i].key, 8);
          WRITE_UNDER(htole64, sc->entries[i].value, 8);
     }
#undef WRITE_UNDER
}

void small_cuckoo_deserialize(int fd, small_cuckoo *sc)
{
     *sc = (small_cuckoo){0};
#define READ(v,n) ENSURE(n == read(fd, v, n))
#define READ_AND(then,v,n) do { uint32_t u = 0; READ(&u,n); v = then(u); } while(0)
     READ_AND(le16toh, sc->n_entries, 2);
     sc->table_size = 1<<(ceil_pow2(sc->n_entries)+1);
     ENSURE(sc->table = malloc(sc->table_size * sizeof sc->table[0]));
     ENSURE(sc->entries = malloc(sc->n_entries * sizeof sc->entries[0]));
     for (uint16_t i = 0; i < sc->n_entries; ++i) {
          READ_AND(le64toh, sc->entries[i].key, 8);
          READ_AND(le64toh, sc->entries[i].value, 8);
     }
#undef READ_AND
#undef READ

     for (uint16_t i = 0; i < sc->n_entries; ++i)
          insert(sc, i);
}

void small_cuckoo_iterate(small_cuckoo *sc, small_cuckoo_iter *iter)
{
     *iter = (small_cuckoo_iter){ .sc = sc, .i = 0 };
}

bool small_cuckoo_iter_has_next(small_cuckoo_iter *iter)
{
     for (; iter->i < iter->sc->table_size; ++iter->i) {
          if (iter->sc->table[iter->i]) return true;
     }
     return false;
}

extern void small_cuckoo_iter_next(small_cuckoo_iter *iter, uint64_t *key, uint64_t *value)
{
     for (; iter->i < iter->sc->table_size; ++iter->i) {
          uint16_t j = iter->sc->table[iter->i];
          if (j) {
               if (key) *key = iter->sc->entries[j].key;
               if (value) *value = iter->sc->entries[j].value;
               ++iter->i;
               return;
          }
     }
     ENSURE(false);
}


#ifdef UNIT_TEST

#include <tap.h>
#include <time.h>

/* Fowler-Noll-Vo hash, per http://isthe.com/chongo/tech/comp/fnv/ */
static uint64_t fnv_hash(uint8_t *data, size_t n)
{
     uint64_t h = 14695981039346656037ULL;
     for (; n > 0; ++data, --n) {
          h ^= *data;
          h += (h<<1) + (h<<4) + (h<<5) + (h<<7) + (h<<8) + (h<<40); /* alternately: h *= 1099511628211ULL; */
     }
     return h;
}

/* The idea for this equation comes from section 7.6 of Aho, Sethi,
 * and Ullmann; Compilers: Principles, Techniques, and Tools (2002).
 */
static double evaluate_hash_quality(uint64_t *b, size_t n)
{
     double sum = 0.0;
     for (size_t j = 0; j < n; ++j) {
          sum += (b[j] * (b[j]+1.0)/2.0);
     }
     return sum/(1.5*(double)n-.5);
}

void test_basic_ops_randomized()
{
     int t = time(NULL);
     note("%s: seed %d", __func__, t);
     srand(t);

     /* Note: should be a power of two for hash quality tests. */
     enum { TEST_BASIC_N_ELEMENTS = 1024 };
     uint64_t keys[TEST_BASIC_N_ELEMENTS] = {0};
     uint64_t values[TEST_BASIC_N_ELEMENTS] = {0};
     uint64_t hash_quality_test[2][TEST_BASIC_N_ELEMENTS] = {{0},{0}};

     small_cuckoo sc;
     sc = small_cuckoo_new(0);
     for (int i = 0; i < TEST_BASIC_N_ELEMENTS; i++) {
          keys[i] = rand();
          keys[i] = fnv_hash((uint8_t *)&keys[i], 8);
          values[i] = rand();
          small_cuckoo_insert(&sc, keys[i], values[i]);
          int n = TEST_BASIC_N_ELEMENTS;
          ++hash_quality_test[0][hash_1(n<<1, keys[i])>>1];
          ++hash_quality_test[1][hash_2(n<<1, keys[i])>>1];
     }

     int success = 1;
     for (int i = 0; i < TEST_BASIC_N_ELEMENTS; i++) {
          uint64_t v;
          success &= small_cuckoo_find(&sc, keys[i], &v);
          success &= v == values[i];
     }
     ok(success, "all keys found and values match");

     small_cuckoo_iter iter;
     small_cuckoo_iterate(&sc, &iter);
     success = 1;
     while (small_cuckoo_iter_has_next(&iter) && success) {
          uint64_t k, v;
          small_cuckoo_iter_next(&iter, &k, &v);
          int i;
          for (i = 0; i < TEST_BASIC_N_ELEMENTS; i++) {
               if (k == keys[i]) {
                    success &= v == values[i];
                    keys[i] = 0;
                    break;
               }
          }
          if (i == TEST_BASIC_N_ELEMENTS) success = 0;
     }
     ok(success, "iterator finds all entries inserted");

     small_cuckoo_free(&sc);

     for (int i = 0; i < 2; ++i) {
          double quality = evaluate_hash_quality(hash_quality_test[i], TEST_BASIC_N_ELEMENTS);
          note("estimated quality of hash %d is %f\n", i+1, quality);
          /* Quality below 0.5 would be great but should be
           * impossible, so we put that test in to catch testing
           * errors. */
          ok(quality > 0.5 && quality < 1.05, "hash quality acceptable");
     }
}

void test_basic_ops_incremental()
{
     note(__func__);

     /* Note: should be a power of two for hash quality tests. */
     enum { TEST_BASIC_N_ELEMENTS = 1024 };
     uint64_t hash_quality_test[2][TEST_BASIC_N_ELEMENTS] = {{0},{0}};

     small_cuckoo sc;
     sc = small_cuckoo_new(0);
     for (int i = 0; i < TEST_BASIC_N_ELEMENTS; i++) {
          small_cuckoo_insert(&sc, i, i);
          int n = TEST_BASIC_N_ELEMENTS;
          ++hash_quality_test[0][hash_1(n<<1, i)>>1];
          ++hash_quality_test[1][hash_2(n<<1, i)>>1];
     }

     int success = 1;
     for (uint32_t i = 0; i < TEST_BASIC_N_ELEMENTS; i++) {
          uint64_t v;
          success &= small_cuckoo_find(&sc, i, &v);
          success &= v == i;
     }
     ok(success, "all keys found and values match");

     small_cuckoo_iter iter;
     small_cuckoo_iterate(&sc, &iter);
     success = 1;
     while (small_cuckoo_iter_has_next(&iter) && success) {
          uint64_t k, v;
          small_cuckoo_iter_next(&iter, &k, &v);
     }
     ok(success, "iterator finds all entries inserted");

     small_cuckoo_free(&sc);

     for (int i = 0; i < 2; ++i) {
          double quality = evaluate_hash_quality(hash_quality_test[i], TEST_BASIC_N_ELEMENTS);
          note("estimated quality of hash %d is %f\n", i+1, quality);
          /* Quality below 0.5 would be great but should be
           * impossible, so we put that test in to catch testing
           * errors. */
          ok(quality > 0.5 && quality < 1.05, "hash quality acceptable");
     }
}

int main()
{
     struct {
          void (*fn)();
          int count;
     } tests[] = {
          {test_basic_ops_randomized, 4},
          {test_basic_ops_incremental, 4}
     };

     int i, count = 0, n = (sizeof tests)/(sizeof tests[0]);
     for (i = 0; i < n; i++)
          count += tests[i].count;
     plan(count, "small-cuckoo");
     for (i = 0; i < n; i++)
          tests[i].fn();
     done_testing();

}

#endif



// ============================================================================
// Type Definitions
// ============================================================================

/*!
 * @struct Pscalar
 */
typedef struct {
  int i;
} Pscalar;

/*!
 * @struct Pvector
 */
typedef struct {
#if dimension == 1
  Pscalar x;
#elif dimension == 2
  Pscalar x;
  Pscalar y;
#else // dimension == 3
  Pscalar x;
  Pscalar y;
  Pscalar z;
#endif
} Pvector;

/*!
 * @struct PAttributes
 */
typedef struct {
  char *name;
  Pvector v;
  bool nodump;
  bool dirty;
} PAttributes;

// ============================================================================
// Global Declarations
// ============================================================================

PAttributes *_pattribute = NULL;
Pscalar *pall = NULL;
size_t _pattribute_len = 0;
size_t npvar = 0;

Pvector ppos;
Pvector pvel;
Pvector pforce;

// ============================================================================
// Function Declarations
// ============================================================================

void init_psolver();
void particle_fields_update_all();
static inline double *_pval(Pscalar s, Particle *n);

Pscalar _init_pscalar(const char *name);

int plist_len(Pscalar *list);
Pscalar *plist_append(Pscalar *list, Pscalar sc);
Pscalar *plist_prepend(Pscalar *list, Pscalar sc);
Pscalar *plist_add(Pscalar *list, Pscalar sc);
int plist_lookup(Pscalar *l, Pscalar s1);
Pscalar *plist_copy(Pscalar *l);
Pscalar *plist_concat(Pscalar *l1, Pscalar *l2);
void plist_print(Pscalar *l, FILE *fp);

Pvector _init_pvector(const char *name);

int pvectors_len(Pvector *list);
Pvector *pvectors_append(Pvector *list, Pvector v);
Pvector *pvectors_add(Pvector *list, Pvector vv);
Pvector *pvectors_copy(Pvector *l);

// ============================================================================
// Macros
// ============================================================================

/*!
 * @def new_pscalar
 */
// clang-format off
@define new_pscalar(name) name = _init_pscalar (#name);
// clang-format on

/*!
 * @def new_pvector
 */
// clang-format off
@define new_pvector(name) name = _init_pvector (#name);
// clang-format on

/*!
 * @def pval
 * @note Assume node exists in scope and is a pointer of Particle
 */
// clang-format off
#define pval(s) (*_pval((s), particle))
#define pval_other(s) (*_pval((s), particle_other))
// clang-format on

/*!
 * @def pname
 * @note Assume node exists in scope and is a pointer of Particle
 */
// clang-format off
#define pname(s) (char*)(_pattribute[s.i].name)
// clang-format on

/*!
 * @def pdirty
 * @note Assume node exists in scope and is a pointer of Particle
 */
// clang-format off
#define pdirty(s) (_pattribute[(s).i].dirty)
// clang-format on

/*!
 * @def pnodump
 * @note Assume node exists in scope and is a pointer of Particle
 */
// clang-format off
#define pnodump(s) (_pattribute[(s).i].nodump)
// clang-format on

/*!
 * @def foreach_pscalar
 */
macro foreach_pscalar(Pscalar *list = pall) {
  {
    Pscalar *_i = (Pscalar *)(list);
    if (_i)
      for (Pscalar s = *_i; (&s)->i >= 0; s = *++_i) {
        // clang-format off
        {...}
        // clang-format on
      }
  }
}

/*!
 * @def foreach_pvector
 */
macro foreach_pvector(Pvector *list) {
  {
    Pvector *_i = (Pvector *)(list);
    if (_i)
      for (Pvector v = *_i; (&v)->x.i >= 0; v = *++_i) {
        // clang-format off
        {...}
        // clang-format on
      }
  }
}

// ============================================================================
// Function definitions
// ============================================================================

/*!
 * @brief
 */
void init_psolver() {
  new_pvector(pvel);
  new_pvector(pforce);
  new_pvector(ppos);

  particle_fields_update_all();
}

void particle_fields_update_all() {
  size_t n = npvar;
  if (pall)
    free(pall);
  pall = (Pscalar *)malloc(sizeof(Pscalar) * (n + 1));
  assert(pall);
  for (size_t i = 0; i < n; i++)
    pall[i].i = i;
  pall[n].i = -1;
}

// ============================================================================
// Pscalar
// ============================================================================

/*!
 * @brief
 *
 * @relates Pscalar
 */
Pscalar _init_pscalar(const char *name) {
  Pscalar s = {.i = npvar++};

  if ((size_t)(s.i + 1) > _pattribute_len) {
    size_t old_len = _pattribute_len;
    _pattribute_len = (size_t)(s.i + 1);
    if (_pattribute == NULL)
      _pattribute = (PAttributes *)calloc(_pattribute_len, sizeof(PAttributes));
    else
      _pattribute = (PAttributes *)realloc(
          _pattribute, _pattribute_len * sizeof(PAttributes));
    assert(_pattribute);
    for (size_t i = old_len; i < _pattribute_len; i++)
      _pattribute[i] = (PAttributes){0};
  }

  if (_pattribute[s.i].name)
    free(_pattribute[s.i].name);
  _pattribute[s.i].name = strdup(name ? name : "");

  foreach_dimension() { _pattribute[s.i].v.x.i = -1; }

  particle_fields_update_all();

  return s;
}

/*!
 * @brief
 * @relates Pscalar
 */
static inline double *_pval(Pscalar s, Particle *p) {
  // if (set_dirty)
  //   _pattribute[s.i].dirty = true;

  return &((double *)((char *)p + sizeof(Particle)))[s.i];
}

/*!
 * @brief
 *
 * @relates Pscalar
 */
int plist_len(Pscalar *list) {
  if (!list)
    return 0;
  int ns = 0;
  foreach_pscalar(list) ns++;
  return ns;
}

/*!
 * @brief
 *
 * @relates Pscalar
 */
Pscalar *plist_append(Pscalar *list, Pscalar sc) {
  int len = plist_len(list);
  qrealloc(list, len + 2, Pscalar);
  list[len] = sc;
  list[len + 1].i = -1;
  return list;
}

/*!
 * @brief
 *
 * @relates Pscalar
 */
Pscalar *plist_prepend(Pscalar *list, Pscalar sc) {
  int len = plist_len(list);
  qrealloc(list, len + 2, Pscalar);
  for (int i = len; i >= 1; i--)
    list[i] = list[i - 1];
  list[0] = sc;
  list[len + 1].i = -1;
  return list;
}

/*!
 * @brief
 * @relates Pscalar
 */
Pscalar *plist_add(Pscalar *list, Pscalar sc) {
  foreach_pscalar(list) {
    if (s.i == sc.i)
      return list;
  }
  return plist_append(list, sc);
}

/*!
 * @brief
 * @relates Pscalar
 */
int plist_lookup(Pscalar *l, Pscalar s1) {
  if (l != NULL)
    foreach_pscalar(l) if (s1.i == s.i) return true;
  return false;
}

/*!
 * @brief
 * @relates Pscalar
 */
Pscalar *plist_copy(Pscalar *l) {
  Pscalar *list = NULL;
  if (l != NULL)
    foreach_pscalar(l) list = plist_append(list, s);
  return list;
}

/*!
 * @brief
 * @relates Pscalar
 */
Pscalar *plist_concat(Pscalar *l1, Pscalar *l2) {
  Pscalar *l3 = plist_copy(l1);
  foreach_pscalar(l2) l3 = plist_append(l3, s);
  return l3;
}

/*!
 * @brief
 * @relates Pscalar
 */
void plist_print(Pscalar *l, FILE *fp) {
  int i = 0;
  foreach_pscalar(l) fprintf(fp, "%s%s", i++ == 0 ? "{" : ",", pname(s));
  fputs(i > 0 ? "}\n" : "{}\n", fp);
}

// ============================================================================
// Pvector
// ============================================================================

/*!
 * @brief
 * @relates Pvector
 */
Pvector _init_pvector(const char *name) {
  struct {
    char *x, *y, *z;
  } ext = {".x", ".y", ".z"};

  Pvector v = {0};

  foreach_dimension() {
    if (name) {
      char cname[strlen(name) + 3];
      strcat(strcpy(cname, name), ext.x);
      v.x = _init_pscalar(cname);
    } else {
      v.x = _init_pscalar(NULL);
    }
  }

  foreach_dimension() { _pattribute[v.x.i].v = v; }

  return v;
}

/*!
 * @brief
 * @relates Pvector
 */
int pvectors_len(Pvector *list) {
  if (!list)
    return 0;
  int nv = 0;
  foreach_pvector(list) nv++;
  return nv;
}

/*!
 * @brief
 * @relates Pvector
 */
Pvector *pvectors_append(Pvector *list, Pvector v) {
  int len = pvectors_len(list);
  qrealloc(list, len + 2, Pvector);
  list[len] = v;
  list[len + 1] = (Pvector){{-1}};
  return list;
}

/*!
 * @brief
 * @relates Pvector
 */
Pvector *pvectors_add(Pvector *list, Pvector vv) {
  foreach_pvector(list) {
    bool id = true;
    foreach_dimension() {
      if (v.x.i != vv.x.i)
        id = false;
      if (id)
        return list;
    }
  }
  return pvectors_append(list, vv);
}

/*!
 * @brief
 * @relates Pvector
 */
Pvector *pvectors_copy(Pvector *l) {
  Pvector *list = NULL;
  if (l != NULL)
    foreach_pvector(l) list = pvectors_append(list, v);

  return list;
}

#pragma once

#include "capsule-mesh-prototype-ft.h"


// ============================================================================
// Type Definitions
// ============================================================================

typedef struct Capsules {
  CapsuleMesh **caps;
  int nbcaps;
  int capacity;
} Capsules;

// ============================================================================
// Globals
// ============================================================================

Capsules allCaps = {0};

// ============================================================================
// Macros
// ============================================================================

#ifndef NCAPS
#define NCAPS 1
#endif
 
#define CAPS(i) (*allCaps.caps[i])
#define CAPS_COUNT() (allCaps.nbcaps)

// ============================================================================
// Function Declarations
// ============================================================================

void free_one_caps(CapsuleMesh *mesh);
void free_all_caps(Capsules *caps);
void initialize_capsules();
void initialize_active_capsule(CapsuleMesh *mesh, int cap_id, int cap_type);

void capsule_manager_reserve(Capsules *caps, int capacity);
void capsule_manager_set_count(Capsules *caps, int count);
CapsuleMesh *capsule_manager_add(Capsules *caps);

// ============================================================================
// Function Definitions
// ============================================================================

/*!
 * @brief
 */
void capsule_manager_reserve(Capsules *caps, int capacity) {
  if (capacity <= caps->capacity)
    return;

  int old_capacity = caps->capacity;
  caps->caps =
      (CapsuleMesh **)realloc(caps->caps, (size_t)capacity * sizeof(CapsuleMesh *));
  assert(caps->caps);
  caps->capacity = capacity;

  for (int i = old_capacity; i < caps->capacity; i++) {
    caps->caps[i] = (CapsuleMesh *)malloc(sizeof(CapsuleMesh));
    assert(caps->caps[i]);
    initialize_empty_capsule(caps->caps[i]);
  }
}

/*!
 * @brief
 */
void capsule_manager_set_count(Capsules *caps, int count) {
  assert(count >= 0);
  capsule_manager_reserve(caps, count > 0 ? count : 1);

  for (int i = count; i < caps->nbcaps; i++) {
    if (caps->caps[i]->isactive)
      free_one_caps(caps->caps[i]);
    initialize_empty_capsule(caps->caps[i]);
  }
  for (int i = caps->nbcaps; i < count; i++)
    initialize_empty_capsule(caps->caps[i]);
  caps->nbcaps = count;
}

/*!
 * @brief
 */
CapsuleMesh *capsule_manager_add(Capsules *caps) {
  if (caps->nbcaps >= caps->capacity) {
    int next_capacity = caps->capacity > 0 ? 2 * caps->capacity : 1;
    capsule_manager_reserve(caps, next_capacity);
  }

  CapsuleMesh *mesh = caps->caps[caps->nbcaps++];
  initialize_empty_capsule(mesh);
  return mesh;
}

/*!
 * @brief
 */
void free_one_caps(CapsuleMesh *mesh) {
  if (!mesh)
    return;
  capsule_mesh_clear_storage(mesh);
  initialize_empty_capsule(mesh);
}

/*!
 * @brief
 */
void free_all_caps(Capsules *caps) {
  for (int i = 0; i < caps->capacity; i++) {
    if (i < caps->nbcaps && caps->caps[i]->isactive)
      free_one_caps(caps->caps[i]);
    free(caps->caps[i]);
  }
  free(caps->caps);
  caps->caps = NULL;
  caps->nbcaps = 0;
  caps->capacity = 0;
}

/*!
 * @brief
 */
void initialize_capsules() {
  if (!allCaps.caps)
    capsule_manager_set_count(&allCaps, NCAPS);
  for (int i = 0; i < allCaps.nbcaps; i++) {
    if (CAPS(i).isactive)
      initialize_empty_capsule(&CAPS(i));
  }
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face() a.x[] = 0.;
  }
}

/*!
 * @brief
 */
void initialize_active_capsule(CapsuleMesh *mesh, int cap_id, int cap_type) {
  initialize_empty_capsule(mesh);
  mesh->cap_type = cap_type;
  mesh->cap_id = cap_id;
  mesh->isactive = true;
}

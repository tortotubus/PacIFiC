#pragma once

#include "particle/Particle.h"

// ============================================================================
// Type definitions
// ============================================================================

/*!
 * @struct ParticleList
 * @brief Non-owning dynamic array of pointers to @ref Particle.
 *
 * @details The container stores borrowed pointers; it does not allocate or
 * destroy the referenced Particle objects. It simply manages the pointer array.
 */
typedef struct {
  Particle **ptrs; /*!< borrowed pointers to nodes */
  size_t size;     /*!< active count */
  size_t capacity; /*!< allocated length of ptrs[] */
} ParticleList;

// ============================================================================
// Function declarations
// ============================================================================

int particle_list_init(ParticleList *list, size_t initial_capacity);
void particle_list_free(ParticleList *list);
int particle_list_reserve(ParticleList *list, size_t new_capacity);
int particle_list_push(ParticleList *list, Particle *particle);
Particle *particle_list_remove_swap(ParticleList *list, size_t i);
static inline Particle *particle_list_get(ParticleList *list, size_t i);
void particle_list_clear(ParticleList *list);
static size_t particle_list_next_cap_(size_t capacity);

// ============================================================================
// Function definitions
// ============================================================================

/*!
 * @brief
 */
static size_t particle_list_next_cap_(size_t capacity) {
  return capacity ? capacity * 2u : 8u;
}

/*!
 * @brief Initialize an ParticleList with an optional initial capacity.
 *
 * @param list Pointer to the ParticleList to initialize. Must not be NULL.
 * @param initial_capacity Initial number of pointers to allocate. If 0, no
 * memory is allocated.
 *
 * @return 0 on success, -1 on failure (invalid input or memory allocation
 * failure).
 *
 * @details Sets the size to 0 and capacity to 0. If initial_capacity > 0,
 * allocates memory for that many Particle pointers.
 *
 * @relates ParticleList
 */
int particle_list_init(ParticleList *list, size_t initial_capacity) {
  if (!list)
    return -1;
  list->ptrs = NULL;
  list->size = 0;
  list->capacity = 0;

  if (initial_capacity == 0)
    return 0;
  list->ptrs = (Particle **)malloc(initial_capacity * sizeof(*list->ptrs));
  if (!list->ptrs)
    return -1;
  list->capacity = initial_capacity;
  return 0;
}

/*!
 * @brief Free all allocated memory and reset the ParticleList.
 *
 * @param list Pointer to the ParticleList to free. If NULL, this function does
 * nothing.
 *
 * @details Deallocates the internal pointer array, resets size and capacity to
 * 0, and sets ptrs to NULL. Does not free the Particle objects themselves.
 *
 * @relates ParticleList
 */
void particle_list_free(ParticleList *list) {
  if (!list)
    return;
  free(list->ptrs);
  list->ptrs = NULL;
  list->size = 0;
  list->capacity = 0;
}

/*!
 * @brief Reserve capacity for at least new_capacity elements in the
 * ParticleList.
 *
 * @param list Pointer to the ParticleList. Must not be NULL.
 * @param new_capacity New capacity to reserve.
 *
 * @return 0 on success or if new_capacity <= current capacity, -1 on failure
 * (memory allocation error).
 *
 * @details If new_capacity is less than or equal to the current capacity, no
 * action is taken. Otherwise, reallocates the internal pointer array to
 * accommodate new_capacity elements.
 *
 * @relates ParticleList
 */
int particle_list_reserve(ParticleList *list, size_t new_capacity) {
  if (!list)
    return -1;
  if (new_capacity <= list->capacity)
    return 0;
  Particle **particle =
      (Particle **)realloc(list->ptrs, new_capacity * sizeof(*particle));
  if (!particle)
    return -1;
  list->ptrs = particle;
  list->capacity = new_capacity;
  return 0;
}

/*!
 * @brief Add an Particle pointer to the end of the ParticleList.
 *
 * @param list Pointer to the ParticleList. Must not be NULL.
 * @param particle Pointer to the Particle to add. Must not be NULL.
 *
 * @return 0 on success, -1 on failure (invalid input or memory allocation
 * failure).
 *
 * @details If the list is at capacity, automatically reserves more capacity
 * (doubling the current capacity, or setting to 8 if empty) before adding the
 * particle.
 *
 * @relates ParticleList
 */
int particle_list_push(ParticleList *list, Particle *particle) {
  if (!list || !particle)
    return -1;
  if (list->size == list->capacity) {
    size_t new_capacity = particle_list_next_cap_(list->capacity);
    if (particle_list_reserve(list, new_capacity) != 0)
      return -1;
  }
  list->ptrs[list->size++] = particle;
  return 0;
}

/*!
 * @brief Remove an element at index i using swap-with-last removal.
 *
 * @param list Pointer to the ParticleList. Must not be NULL.
 * @param i Index of the element to remove. Must be < size.
 *
 * @return Pointer to the removed Particle, or NULL if invalid input.
 *
 * @details Swaps the element at index i with the last element and decrements
 * size. This is efficient but does not preserve list order. The caller is
 * responsible for deciding how or if to destroy the returned Particle object.
 *
 * @relates ParticleList
 */
Particle *particle_list_remove_swap(ParticleList *list, size_t i) {
  if (!list || i >= list->size)
    return NULL;
  Particle *dead = list->ptrs[i];
  list->ptrs[i] = list->ptrs[list->size - 1];
  list->size--;
  return dead; /* caller decides how/if to destroy it */
}

/*!
 * @brief Retrieve the Particle pointer at index i.
 *
 * @param list Pointer to the ParticleList. If NULL, returns NULL.
 * @param i Index of the element to retrieve. Must be < size.
 *
 * @return Pointer to the Particle at index i, or NULL if invalid input or out
 * of bounds.
 *
 * @relates ParticleList
 */
Particle *particle_list_get(ParticleList *list, size_t i) {
  return (list && i < list->size) ? list->ptrs[i] : NULL;
}

/*!
 * @brief Clear the list without freeing its storage.
 *
 * @param list Pointer to the ParticleList. If NULL, this function does nothing.
 *
 * @details Sets size to 0 and leaves capacity/ptrs unchanged. Does not free
 * any Particle objects (non-owning container).
 *
 * @relates ParticleList
 */
void particle_list_clear(ParticleList *list) {
  if (!list)
    return;
  list->size = 0;
}
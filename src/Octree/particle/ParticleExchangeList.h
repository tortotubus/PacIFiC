#pragma once

#include "particle/ParticleList.h"

// ============================================================================
// Type definitions
// ============================================================================

/**
 * @struct ParticleExchangeList
 * @brief Non-owning list of particles exchanged with one peer MPI rank.
 */
typedef struct ParticleExchangeList {
  int pid; /**< peer rank, -1 until first insertion */
  unsigned char *buffs;
  size_t nbuff;
  ParticleList particles;
} ParticleExchangeList;

// ============================================================================
// Function declarations
// ============================================================================

int particle_exchange_list_init(ParticleExchangeList *list,
                                size_t initial_capacity);
void particle_exchange_list_free(ParticleExchangeList *list);
void particle_exchange_list_clear(ParticleExchangeList *list);
int particle_exchange_list_set_pid(ParticleExchangeList *list, int pid);
int particle_exchange_list_push(ParticleExchangeList *list, Particle *particle);
int particle_exchange_list_push_unique(ParticleExchangeList *list,
                                       Particle *particle);
void particle_exchange_list_init_buffer(ParticleExchangeList *list,
                                        size_t nbuff);
void particle_exchange_list_free_buffer(ParticleExchangeList *list);

// ============================================================================
// Function definitions
// ============================================================================

/*!
 * @brief
 */
int particle_exchange_list_init(ParticleExchangeList *list,
                                size_t initial_capacity) {
  if (!list)
    return -1;
  list->pid = -1;
  list->buffs = NULL;
  list->nbuff = 0;
  return particle_list_init(&list->particles, initial_capacity);
}

/*!
 * @brief
 */
void particle_exchange_list_free(ParticleExchangeList *list) {
  if (!list)
    return;
  particle_exchange_list_free_buffer(list);
  particle_list_free(&list->particles);
  list->pid = -1;
}

/*!
 * @brief
 */
void particle_exchange_list_clear(ParticleExchangeList *list) {
  if (!list)
    return;
  particle_list_clear(&list->particles);
}

/*!
 * @brief
 */
int particle_exchange_list_set_pid(ParticleExchangeList *list, int pid) {
  if (!list)
    return -1;
  if (list->pid < 0) {
    list->pid = pid;
    return 0;
  }
  return (list->pid == pid) ? 0 : -1;
}

/*!
 * @brief
 */
int particle_exchange_list_push(ParticleExchangeList *list,
                                Particle *particle) {
  if (!list || !particle)
    return -1;

  return particle_list_push(&list->particles, particle);
}

/*!
 * @brief
 */
int particle_exchange_list_push_unique(ParticleExchangeList *list,
                                       Particle *particle) {
  if (!list || !particle)
    return -1;

  for (size_t i = list->particles.size; i > 0; i--) {
    if (list->particles.ptrs[i - 1] == particle)
      return 0;
  }

  return particle_list_push(&list->particles, particle);
}

/*!
 * @brief
 */
void particle_exchange_list_init_buffer(ParticleExchangeList *list,
                                        size_t nbuffs) {
  if (!list)
    return;
  particle_exchange_list_free_buffer(list);
  list->nbuff = nbuffs;
  if (nbuffs > 0) {
    list->buffs = (unsigned char *)calloc(nbuffs, sizeof(*list->buffs));
    assert(list->buffs);
  }
}

/*!
 * @brief
 */
void particle_exchange_list_free_buffer(ParticleExchangeList *list) {
  if (!list)
    return;
  free(list->buffs);
  list->buffs = NULL;
  list->nbuff = 0;
}

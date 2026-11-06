#ifndef GRAINS_TIMER_BACKEND_HH
#define GRAINS_TIMER_BACKEND_HH

namespace grains::detail
{
    /** Synchronize GPU work for wall-clock timing; a no-op in CPU-only builds. */
    void synchronizeTimerDevice() noexcept;
}

#endif

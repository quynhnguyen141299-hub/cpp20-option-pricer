#pragma once
/// @file arena.hpp
/// Monotonic arena allocator for path-generation scratch memory.
///
/// Why: The MC hot loop allocates/frees vectors of doubles millions of times.
/// malloc/free overhead + cache pollution from scattered heap allocations
/// dominates at high path counts.  A per-thread arena pre-allocates one
/// contiguous block and bumps a pointer — O(1) alloc, zero fragmentation,
/// cache-line-friendly sequential access.
///
/// RAII: the arena owns its memory and frees it in the destructor.
/// Move-only: avoids double-free on the backing buffer.

#include <cstddef>
#include <cstdlib>
#include <memory>
#include <cassert>
#include <new>
#include <span>

namespace pricer {

class Arena {
public:
    /// Construct with capacity in bytes.  Default 4 MB per thread.
    explicit Arena(std::size_t capacity = 4ULL * 1024 * 1024)
        : capacity_(capacity)
        , buf_(static_cast<std::byte*>(std::aligned_alloc(64, capacity)))  // cache-line aligned
        , offset_(0)
    {
        if (!buf_) throw std::bad_alloc{};
    }

    // RAII: free on destruction
    ~Arena() { std::free(buf_); }

    // Move-only
    Arena(Arena&& o) noexcept
        : capacity_(o.capacity_), buf_(o.buf_), offset_(o.offset_)
    {
        o.buf_ = nullptr;
        o.capacity_ = 0;
        o.offset_ = 0;
    }

    Arena& operator=(Arena&& o) noexcept {
        if (this != &o) {
            std::free(buf_);
            capacity_ = o.capacity_;
            buf_      = o.buf_;
            offset_   = o.offset_;
            o.buf_      = nullptr;
            o.capacity_ = 0;
            o.offset_   = 0;
        }
        return *this;
    }

    Arena(const Arena&) = delete;
    Arena& operator=(const Arena&) = delete;

    /// Allocate n doubles from the arena.  Returns a span for bounds safety.
    [[nodiscard]] std::span<double> alloc_doubles(std::size_t n) noexcept {
        constexpr std::size_t align = alignof(double);
        std::size_t aligned_offset = (offset_ + align - 1) & ~(align - 1);
        std::size_t bytes = n * sizeof(double);

        if (aligned_offset + bytes > capacity_) {
            return {};  // caller must handle: arena exhausted
        }

        auto* ptr = reinterpret_cast<double*>(buf_ + aligned_offset);
        offset_ = aligned_offset + bytes;
        return {ptr, n};
    }

    /// Reset to reuse all memory (no dealloc, just rewind the pointer).
    void reset() noexcept { offset_ = 0; }

    [[nodiscard]] std::size_t used()     const noexcept { return offset_; }
    [[nodiscard]] std::size_t capacity() const noexcept { return capacity_; }

private:
    std::size_t capacity_;
    std::byte*  buf_;
    std::size_t offset_;
};

} // namespace pricer

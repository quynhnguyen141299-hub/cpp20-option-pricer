#pragma once
/// @file thread_pool.hpp
/// Minimal thread pool using std::jthread and stop_token.
///
/// Why std::jthread over raw threads:
///   - Cooperative cancellation via stop_token (no need for atomic flags)
///   - Auto-join in destructor (RAII — no detached zombie threads)
///   - Simpler shutdown path than manual CV + flag patterns
///
/// Design: fixed-size pool, lock-free-ish task queue (std::mutex + deque).
/// Each worker checks stop_token before dequeuing.  Submitting returns
/// a std::future for the result.

#include <thread>
#include <mutex>
#include <condition_variable>
#include <deque>
#include <functional>
#include <future>
#include <vector>
#include <type_traits>

namespace pricer {

class ThreadPool {
public:
    explicit ThreadPool(unsigned n_threads = 0) {
        if (n_threads == 0)
            n_threads = std::max(1u, std::thread::hardware_concurrency());

        workers_.reserve(n_threads);
        for (unsigned i = 0; i < n_threads; ++i) {
            workers_.emplace_back([this](std::stop_token stoken) {
                while (!stoken.stop_requested()) {
                    std::function<void()> task;
                    {
                        std::unique_lock lock(mtx_);
                        cv_.wait(lock, stoken, [this] { return !tasks_.empty(); });
                        if (stoken.stop_requested() && tasks_.empty()) return;
                        if (tasks_.empty()) continue;
                        task = std::move(tasks_.front());
                        tasks_.pop_front();
                    }
                    task();
                }
            });
        }
    }

    // RAII: jthreads request stop + join automatically in destructor
    ~ThreadPool() {
        {
            std::scoped_lock lock(mtx_);
            // Drain remaining tasks
        }
        for (auto& w : workers_) w.request_stop();
        cv_.notify_all();
        // jthread destructor joins automatically
    }

    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    /// Submit a callable, get a future back.
    template <typename F, typename... Args>
        requires std::invocable<F, Args...>
    [[nodiscard]] auto submit(F&& f, Args&&... args)
        -> std::future<std::invoke_result_t<F, Args...>>
    {
        using R = std::invoke_result_t<F, Args...>;
        auto task = std::make_shared<std::packaged_task<R()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );
        auto fut = task->get_future();
        {
            std::scoped_lock lock(mtx_);
            tasks_.emplace_back([task = std::move(task)] { (*task)(); });
        }
        cv_.notify_one();
        return fut;
    }

    [[nodiscard]] unsigned size() const noexcept {
        return static_cast<unsigned>(workers_.size());
    }

private:
    std::vector<std::jthread>          workers_;
    std::deque<std::function<void()>>  tasks_;
    std::mutex                         mtx_;
    std::condition_variable_any        cv_;
};

} // namespace pricer

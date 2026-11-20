// Copyright (C) 2012 Jakob Progsch, Václav Zeman
// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This file implements the methods of thead_pool.hpp

#include "dca/parallel/stdthread/thread_pool/thread_pool.hpp"

#include "mpi.h"

#include "dca/parallel/stdthread/thread_pool/affinity.hpp"

namespace dca {
namespace parallel {

ThreadPool::ThreadPool(size_t n_threads) : stop_(false), active_id_(0) {
  core_count_ = get_core_count();
  master_affinity_ = get_affinity();

  enlarge(n_threads);
}

ThreadPool::~ThreadPool() {
  stop_ = true;
  for (auto& condition : condition_)
    condition->notify_one();
  for (std::thread& worker : workers_)
    worker.join();
}

void ThreadPool::enlarge(size_t n_threads) {
  if (workers_.size() >= n_threads)
    return;

  if (n_threads > n_threads_)
    n_threads_ = n_threads;

  // Lock all queues to avoid modifications while enlarging the pool.
  std::vector<std::unique_lock<std::mutex>> locks;
  for (auto& mutex_ptr : queue_mutex_)
    locks.emplace_back(std::unique_lock<std::mutex>(*mutex_ptr));

  // Create enough resources for the new workers.
  queue_mutex_.reserve(n_threads);
  condition_.reserve(n_threads);
  tasks_.reserve(n_threads);

  // Release the locks on the existing queues.
  locks.clear();

  for (size_t id = workers_.size(); id < n_threads; ++id) {
    queue_mutex_.emplace_back(std::make_unique<std::mutex>());
    condition_.emplace_back(std::make_unique<std::condition_variable>());
    tasks_.emplace_back(std::make_unique<std::queue<std::packaged_task<void()>>>());

    // Start a loop on each new thread.
    workers_.emplace_back(&ThreadPool::workerLoop, this, id);
  }
}

void ThreadPool::workerLoop(int id) {
  // Give thread 'id' its own block of the cores of this process:
  // [id * core_count_ / n_threads_, (id + 1) * core_count_ / n_threads_).
  //
  // Only the cores of this process are used: shifting them by 'core_count_ * id' places the thread
  // on the cores of another rank as soon as a rank owns a few cores only. With a single core there
  // is nothing to split, so the affinity is left as it is.
  if (core_count_ > 1 && n_threads_ > 1) {
    const int n_threads = static_cast<int>(n_threads_);
    const int first = id * core_count_ / n_threads;
    const int last = (id + 1) * core_count_ / n_threads;

    std::vector<int> affinities;
    for (int i = first; i < last && i < static_cast<int>(master_affinity_.size()); ++i)
      affinities.push_back(master_affinity_[i]);

    if (!affinities.empty())
      set_affinity(affinities);
  }

  while (true) {
    std::packaged_task<void()> task;

    {  // Acquire new task.
      std::unique_lock<std::mutex> lock(*queue_mutex_[id]);
      condition_[id]->wait(lock, [this, id] { return stop_ || !tasks_[id]->empty(); });
      // If all the work is done and no more will be posted, return.
      if (stop_ && tasks_[id]->empty())
        return;
      task = std::move(tasks_[id]->front());
      tasks_[id]->pop();
    }
    // Execute task.
    task();
  }
}

}  // namespace parallel
}  // namespace dca

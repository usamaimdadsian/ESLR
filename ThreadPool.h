#ifndef THREADPOOL_H
#define THREADPOOL_H

#include <vector>
#include <queue>
#include <atomic>
#include <future>
#include <condition_variable>
#include <thread>
#include <functional>
#include <stdexcept>

#define  THREADPOOL_MAX_NUM 16
//#define  THREADPOOL_AUTO_GROW

class ThreadPool
{
	using Task = std::function<void()>;
	std::vector<std::thread> _pool;
	std::queue<Task> _tasks;
	std::mutex _lock;
	std::condition_variable _task_cv;
	std::atomic<bool> _stop{ false };

public:
	inline ThreadPool(unsigned short size = 4) { addThread(size); }
	inline ~ThreadPool()
	{
		_stop.store(true);
		_task_cv.notify_all();
		for (std::thread& thread : _pool) {
			if (thread.joinable())
				thread.join();
		}
	}

public:

	template<class F, class... Args>
	auto commit(F&& f, Args&&... args) ->std::future<decltype(f(args...))>
	{
		if (_stop.load())
			throw std::runtime_error("commit on ThreadPool is stopped.");

		using RetType = decltype(f(args...)); // typename std::result_of<F(Args...)>::type, 
		auto task = std::make_shared<std::packaged_task<RetType()>>(
			std::bind(std::forward<F>(f), std::forward<Args>(args)...)
			);
		std::future<RetType> future = task->get_future();
		{
			std::lock_guard<std::mutex> lock{ _lock };
			_tasks.emplace([task]() {
				(*task)();
				});
		}

		_task_cv.notify_one();

		return future;
	}

private:
	void addThread(unsigned short size)
	{
		for (; _pool.size() < THREADPOOL_MAX_NUM && size > 0; --size)
		{
			_pool.emplace_back([this]
				{
					while (!_stop.load())
					{
						Task task;

						{
							std::unique_lock<std::mutex> lock{ this->_lock };
							this->_task_cv.wait(lock, [this] {
								return this->_stop.load() || !this->_tasks.empty();
								});
							if (this->_stop.load() && this->_tasks.empty())
								return;
							task = move(this->_tasks.front());
							_tasks.pop();
						}
						task();
					}
				});
		}
	}
};

#endif // !THREADPOOL_H


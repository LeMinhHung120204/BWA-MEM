import threading
from dataclasses import dataclass


@dataclass
class ktp_worker_t:
    pl: "ktp_t"
    index: int
    step: int = 0
    data: any = None

class ktp_t:
    def __init__(self, shared, func, n_workers, n_steps):
        self.shared = shared
        self.func = func
        self.index = 0
        self.n_workers = n_workers
        self.n_steps = n_steps
        self.workers = [ktp_worker_t(self, i, 0, None) for i in range(n_workers)]
        self.mutex = threading.Lock()
        self.cv = threading.Condition(self.mutex)

def ktp_worker(data : ktp_worker_t):
    w = data  # ktp_worker_t instance
    p = w.pl  # ktp_t instance
    while w.step < p.n_steps:
        with p.mutex:
            while any(p.workers[i].step <= w.step and p.workers[i].index < w.index for i in range(p.n_workers) if w != p.workers[i]):
                p.cv.wait()
        
        # working on w->step
        w.data = p.func(p.shared, w.step, w.data if w.step else None)
        
        # update step and notify other workers
        with p.mutex:
            w.step = (w.step + 1) % p.n_steps if w.step == p.n_steps - 1 or w.data else p.n_steps
            if w.step == 0:
                w.index = p.index
                p.index += 1
            p.cv.notify_all()

def kt_pipeline(n_threads: int, func, shared_data, n_steps: int):
    if n_threads < 1:
        n_threads = 1
    
    aux = ktp_t(shared_data, func, n_threads, n_steps)
    
    threads = [threading.Thread(target=ktp_worker, args=(worker,)) for worker in aux.workers]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    
    with aux.mutex:
        aux.mutex = None
        aux.cv = None

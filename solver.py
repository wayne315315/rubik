import itertools
from multiprocessing import Process, Queue, cpu_count

import numpy as np
from tqdm import tqdm

from r3 import index
from r3 import rs

# DFS
def dfs(index_q, level=8):
    stack = [(index_q, ())]
    ans = None

    pbar = tqdm()
    pbar.update(1)
    while stack:
        indices, seq = stack.pop()
        for r in rs:
            seq_ = seq + (r,)
            indices_ = r(indices)
            if np.all(indices_ == index):
                ans = seq_
                break
            if len(seq_) <= level:
                pbar.update(1)
                stack.append((indices_, seq_))
    pbar.close()
    return ans


def generate_seq(rs, level_max, verbose=False):
    for level in range(level_max+1):
        if verbose:
            print(f"Level: {level}")
        for seq in itertools.product(rs, repeat=level):
            yield seq

def worker(i, q2, index_q, level_max, verbose):
    n = cpu_count()
    seq_gen = enumerate(generate_seq(rs, level_max)) if i else tqdm(enumerate(generate_seq(rs, level_max, verbose=verbose)))
    for j, seq in seq_gen:
        if j % n == i:
            indices = index_q.copy()
            for r in seq:
                indices = r(indices)
            if np.all(indices == index):
                q2.put(seq)
                break

def util(q2, q3):
    while True:
        seq = q2.get()
        if seq is None:
            break
        q3.put(seq)


def brute_force_multi(index_q, level_max=60, verbose=False):
    q2 = Queue()
    q3 = Queue()
    ps = [Process(target=worker, args=(i, q2, index_q, level_max, verbose) ) for i in range(cpu_count())]
    p_util = Process(target=util, args=(q2, q3))
    p_util.start()
    for p in ps:
        p.start()

    while all(p.is_alive() for p in ps):
        pass

    for p in ps:
        p.terminate()
    for p in ps:
        p.join()

    ans = []
    while True:
        try:
            seq = q3.get_nowait()
        except:
            q2.put(None)
            break
        else:
            ans.append(seq)
    p_util.join()
    print("")

    return ans

if __name__ == "__main__":
    import random
    k = random.randrange(3,7)
    seq = random.choices(rs, k=k)
    print("Question", k, seq)
    index_q = index
    for r in seq:
        index_q = r(index_q)
    ###
    answers =  brute_force_multi(index_q)
    ###
    
    print("")
    print("Answer")
    for ans in answers:
        indices = index_q.copy()
        for r in ans:
            indices = r(indices)
        print(ans, np.all(indices == index))

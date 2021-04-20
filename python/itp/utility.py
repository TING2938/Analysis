#!/usr/bin/env python3
import concurrent.futures as futures
import subprocess as sub
import os


def replaceText(before, after, fnm):
    os.system(f"sed -i 's/{before}/{after}/g' {fnm}")


def writeToFile(fnm, content):
    file = open(fnm, "w")
    file.write(content)
    file.close()


def run(cmd: str, check:bool=False):
    h = sub.run(cmd, shell=True, stdout=sub.PIPE, stderr=sub.PIPE, check=check)
    return h.stdout.decode('utf-8')


def processExec(_func, _iter):
    pool = futures.ProcessPoolExecutor(max_workers=len(_iter))
    all_tasks = []
    for i in _iter:
        all_tasks.append(pool.submit(_func, i))
    futures.wait(all_tasks)


def threadExec(_func, _iter):
    pool = futures.ThreadPoolExecutor(max_workers=len(_iter))
    all_tasks = []
    for i in _iter:
        all_tasks.append(pool.submit(_func, i))
    futures.wait(all_tasks)


def normalExec(_func, _iter):
    for i in _iter:
        _func(i)


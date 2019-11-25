from multiprocessing import Process
import os


# from https://stackoverflow.com/questions/49123439/python-how-to-run-process-in-detached-mode
def detachify(func):
    """Decorate a function so that its calls are async in a detached process.

    Usage
    -----

    .. code::
            import time

            @detachify
            def f(message):
                time.sleep(5)
                print(message)

            f('Async and detached!!!')

    """

    # create a process fork and run the function
    def forkify(*args, **kwargs):
        if os.fork() != 0:
            return
        func(*args, **kwargs)

    # wrapper to run the forkified function
    def wrapper(*args, **kwargs):
        proc = Process(target=lambda: forkify(*args, **kwargs))
        proc.start()
        proc.join()
        return

    return wrapper

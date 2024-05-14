from multiprocessing import Pool
from functools import partial
from tqdm import tqdm


def ptqdm(function, iterable, processes, zipped=False, chunksize=1, desc=None, disable=False, **kwargs):
    """
    Run a function in parallel with a tqdm progress bar and an arbitrary number of iterables and arguments.
    Multiple iterables can be packed into a tuple and passed to the 'iterable argument'. The iterables must be the first arguments in the function that is run in parallel.
    Results are always ordered and the performance is the same as of Pool.map.
    :param function: The function that should be parallelized.
    :param iterable: The iterable passed to the function.
    :param processes: The number of processes used for the parallelization.
    :param zipped: If multiple iterables are packed into a tuple. The iterables will be unpacked and passed as separate arguments to the function.
    :param chunksize: The iterable is based on the chunk size chopped into chunks and submitted to the process pool as separate tasks.
    :param desc: The description displayed by tqdm in the progress bar.
    :param disable: Disables the tqdm progress bar.
    :param kwargs: Any additional arguments that should be passed to the function.
    """
    if kwargs:
        function_wrapper = partial(wrapper, function=function, zipped=zipped, **kwargs)
    else:
        function_wrapper = partial(wrapper, function=function, zipped=zipped)

    if zipped:
        length = len(iterable[0])
        iterable = zip(*iterable)
    else:
        length = len(iterable)

    results = [None] * length
    with Pool(processes=processes) as p:
        with tqdm(desc=desc, total=length, disable=disable) as pbar:
            for i, result in p.imap_unordered(function_wrapper, enumerate(iterable), chunksize=chunksize):
                results[i] = result
                pbar.update()
    return results


def wrapper(enum_iterable, function, zipped, **kwargs):
    i = enum_iterable[0]
    if zipped:
        result = function(*enum_iterable[1], **kwargs)
    else:
        result = function(enum_iterable[1], **kwargs)
    return i, result
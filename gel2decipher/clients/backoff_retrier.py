import requests
import urllib2
import logging
import time
import random


def wrapper(func, retries):
    """
    This wrapper implements a truncated binary exponential backoff algorithm between retries.
    (https://en.wikipedia.org/wiki/Exponential_backoff#Binary_exponential_backoff)
    It only captures exceptions raised by the packages requests and urllib2:
    * requests.exceptions.RequestException
    * requests.exceptions.ConnectionError
    * urllib2.URLError
    Other exceptions will override any retries.

    :param func:       the wrapped function
    :param retries:    the maximum number of retries. -1 are infinite retries
    :return:           the return of the wrapped function if any
    """

    def retry(*args, **kwargs):

        results = None
        retries_count = 0
        backoff_iteration = 1
        truncate_iteration = 8
        success = False
        while not success:
            try:
                results = func(*args, **kwargs)
                success = True
            except (requests.exceptions.ConnectionError, requests.exceptions.Timeout,
                    requests.exceptions.TooManyRedirects, urllib2.URLError), ex:
                # do we want here any other requests.exception??
                logging.error(str(ex))
                # retries a fixed number of times
                if retries != -1 and retries_count >= retries:
                    raise ex
                retries_count += 1
                # waits for an increasing random time
                random_sleep = random.randrange(0, (2 ** backoff_iteration) - 1)
                logging.debug("Retrying connection after %s seconds" % str(random_sleep))
                time.sleep(random_sleep)
                # when it reaches the maximum value that it may wait it stops increasing time
                if backoff_iteration < truncate_iteration:
                    backoff_iteration += 1
        return results

    return retry

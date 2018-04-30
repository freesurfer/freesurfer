import datetime
import logging

logger = logging.getLogger(__name__)

class ProcessTimer:
    def __init__(self, message=None):
        if message is None:
            message = "process start"
        logger.info(message)
        self.start_time = datetime.datetime.now()

    @property
    def elapsed_time(self):
        return datetime.datetime.now() - self.start_time

    def mark_time(self, message):
        logger.info('%s:%s', message, str(self.elapsed_time))

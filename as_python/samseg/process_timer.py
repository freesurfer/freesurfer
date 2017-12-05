import datetime
import logging

logger = logging.getLogger(__name__)

class ProcessTimer:
    def __init__(self):
        self.start_time = datetime.datetime.now()

    def mark_time(self, message):
        elapsed_time = datetime.datetime.now() - self.start_time
        logger.info('%s:%s', message, str(elapsed_time))

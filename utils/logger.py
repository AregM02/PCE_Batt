import logging

def create_logger(name: str, level: str = 'DEBUG'):
    """Creates a logger with the given name and level."""

    log_level = getattr(logging, level)

    logger = logging.getLogger(name)
    logger.setLevel(log_level)
    console_handler = logging.StreamHandler()
    logger.addHandler(console_handler)

    return logger
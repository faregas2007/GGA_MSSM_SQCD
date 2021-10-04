import logging.config
import sys
from pathlib import Path

#import pretty_errors
#from rich.logging import RichHandler


# Directories
base_dir = Path(__file__).parent.parent.absolute()
print(base_dir)
config_dir = Path(base_dir, 'config')
input_dir = Path(base_dir, 'input') 
amp_dir =  Path(base_dir, 'amp')
working_dir = Path(base_dir, 'working')
core_dir = Path(base_dir, 'core')
logs_dir = Path(base_dir, 'logs')

# create dirs
amp_dir.mkdir(parents=True, exist_ok=True)
config_dir.mkdir(parents=True, exist_ok=True)
core_dir.mkdir(parents=True, exist_ok=True)
working_dir.mkdir(parents=True, exist_ok=True)
input_dir.mkdir(parents=True, exist_ok=True)
logs_dir.mkdir(parents=True, exist_ok=True)

"""
# logger
logging_config={
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {
        "minimal": {"format", "%{message}s"},
        "detailed": {
            "format":"%(levelname)s %(asctime)s [%(filename)s:%(funcName)s:%(lineno)d]\n%(message)s\n"
        },
    },
    "handlers":{
        "console":{
            "class": "logging.StreamHandler",
            "stream": sys.stdout,
            "formatter": "minimal",
            "level": logging.DEBUG,
        },
        "info": {
            "class": "logging.handlers.RotatingFileHandler",
            "filename":Path(logs_dir, "info.log"),
            "maxBytes":10485760,
            "backupCount":10,
            "formatter": "detailed",
            "level":logging.INFO,
        },
        "error": {
            "class": "logging.handlers.RotatingFileHandler",
            "filename":Path(logs_dir, "error.log"),
            "maxBytes":10385760,
            "backupCount":10,
            "formatter":"detailed",
            "level":logging.ERROR,
        },
    },
    "loggers": {
        "root": {
            "handlers": ["console", "info", "error"],
            "level": logging.INFO,
            "propagate":True,
        },
    },
}
logging.config.dictConfig(logging_config)
logger = logging.getLogger("root")
logger.handlers[0] = RichHandler(markup=True)
"""
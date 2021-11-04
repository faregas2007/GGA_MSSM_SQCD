import sys, os
from pathlib import Path

import logging
import logging.config

base_dir = Path(__file__).parent.parent.absolute()
config_dir = Path(base_dir, 'config')
src_dir = Path(base_dir, 'src')
model_dir = Path(base_dir, 'model')
data_dir = Path(base_dir, 'data')
working_dir = Path(base_dir, 'working')
logs_dir = Path(base_dir, 'logs')
test_dir = Path(base_dir, 'tests')
xs_dir = Path(src_dir, 'xs')
rens_dir = Path(src_dir, 'renschemes')
amp_dir = Path(base_dir, 'amp')
# create dir
config_dir.mkdir(parents=True, exist_ok=True)
src_dir.mkdir(parents=True, exist_ok=True)
model_dir.mkdir(parents=True, exist_ok=True)
data_dir.mkdir(parents=True, exist_ok=True)
working_dir.mkdir(parents=True, exist_ok=True)
logs_dir.mkdir(parents=True, exist_ok=True)
test_dir.mkdir(parents=True, exist_ok=True)
xs_dir.mkdir(parents=True, exist_ok=True)
rens_dir.mkdir(parents=True, exist_ok=True)
amp_dir.mkdir(parents=True, exist_ok=True)

class CustomFilter(logging.Filter):

    COLOR = {
        "DEBUG":"GREEN",
        "INFO":"BLUE",
        "WARNING":"YELLOW",
        "ERROR":"RED",
        "CRITICAL":"RED",
    }

    def filter(self, record):
        record.color = CustomFilter.COLOR[record.levelname]
        return True
 
logging.basicConfig(
    filename = Path(logs_dir, 'runs.log'),
    filemode = 'w',
    level=logging.DEBUG,
    format="%(asctime)s - [%(levelname)s] - [%(color)s] - %(name)s - (%(filename)s).%(funcName)s(%(lineno)d) - %(message)s",
)

logger = logging.getLogger(__name__)
logger.addFilter(CustomFilter())

logger.debug('debug message, color is green')
logger.info('info message, color is blue')
logger.warning('warning message, color is yellow')
logger.error('error message, color is red')
logger.critical('error message, color is red')
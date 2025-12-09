"""
The QC module includes quality control functions from Pecos, see https://pecos.readthedocs.io for more details.
"""

from pecos.monitoring import (
    check_timestamp,
    check_missing,
    check_corrupt,
    check_range,
    check_delta,
    check_outlier,
)

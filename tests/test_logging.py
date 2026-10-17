import logging
import math

from pyXSteam.XSteam import XSteam


def test_package_logging_is_namespaced_and_uses_null_handler(caplog):
    package_logger = logging.getLogger("pyXSteam")
    assert any(
        isinstance(handler, logging.NullHandler)
        for handler in package_logger.handlers
    )

    with caplog.at_level(logging.DEBUG, logger="pyXSteam"):
        steam_table = XSteam()
        result = steam_table.tsat_p(-1.0)

    assert math.isnan(result)
    assert steam_table.logger.name == "pyXSteam.XSteam"
    assert any(
        record.name == "pyXSteam.XSteam"
        and record.levelno == logging.DEBUG
        and "out of range" in record.getMessage()
        for record in caplog.records
    )
    assert not any(
        record.levelno >= logging.WARNING and "out of range" in record.getMessage()
        for record in caplog.records
    )

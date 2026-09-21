"""The library must announce exports through `logging` without touching the app's logging setup."""
import logging
import subprocess
import sys

from hypercomplex import build_table, export_csv

TABLE_ARGS = dict(mode="integer", csv_mode="long")


def test_export_returns_the_path_and_writes_the_file(tmp_path):
    path = str(tmp_path / "t.csv")
    assert export_csv(path, build_table("standard", 2), **TABLE_ARGS) == path
    assert (tmp_path / "t.csv").read_text().startswith("i,j,sign,index")


def test_message_is_logged_on_the_package_logger_at_info(tmp_path, caplog):
    path = str(tmp_path / "t.csv")
    with caplog.at_level(logging.INFO, logger="hypercomplex"):
        export_csv(path, build_table("standard", 2), **TABLE_ARGS)
    records = [r for r in caplog.records if r.name.startswith("hypercomplex")]
    assert len(records) == 1
    assert records[0].levelno == logging.INFO
    assert path in records[0].getMessage()


def test_export_does_not_touch_root_logger_configuration(tmp_path):
    root = logging.getLogger()
    handlers, level = list(root.handlers), root.level
    export_csv(str(tmp_path / "t.csv"), build_table("standard", 2), **TABLE_ARGS)
    assert list(root.handlers) == handlers
    assert root.level == level


def test_package_logger_has_a_null_handler():
    handlers = logging.getLogger("hypercomplex").handlers
    assert any(isinstance(h, logging.NullHandler) for h in handlers)


def test_clean_interpreter_stays_silent_and_unconfigured(tmp_path):
    """Fresh process, no app logging setup: nothing printed, root logger untouched."""
    code = (
        "import logging, sys\n"
        "from hypercomplex import build_table, export_csv\n"
        "root = logging.getLogger()\n"
        "before = (len(root.handlers), root.level)\n"
        f"export_csv({str(tmp_path / 'c.csv')!r}, build_table('standard', 2), mode='integer', csv_mode='long')\n"
        "after = (len(root.handlers), root.level)\n"
        "print('UNCHANGED' if before == after else f'CHANGED {before}->{after}')\n"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert result.stdout.strip() == "UNCHANGED"
    assert result.stderr == ""


def test_application_can_opt_in_to_the_message(tmp_path):
    """The documented way for an application to see the message."""
    code = (
        "import logging\n"
        "logging.basicConfig(level=logging.INFO, format='%(name)s: %(message)s')\n"
        "from hypercomplex import build_table, export_csv\n"
        f"export_csv({str(tmp_path / 'c.csv')!r}, build_table('standard', 2), mode='integer', csv_mode='long')\n"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert "Table exported to" in result.stderr
    assert "hypercomplex.printer.cd_table_printer" in result.stderr

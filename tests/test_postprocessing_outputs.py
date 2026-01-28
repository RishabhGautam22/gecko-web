import os
import tempfile

from gecko_web.postprocessing import run_postprocessing


def _write_min_dictionary(path: str):
    # Minimal dictionary line: code formula class MW flag nC nH nN nO nS
    content = """
! Minimal dictionary
TST C1H2O1 1 46.0 0 1 2 0 1 0
TST2 C2H4O1 1 72.0 0 2 4 0 1 0
""".strip()
    with open(path, "w") as f:
        f.write(content + "\n")


def test_postprocessing_generates_outputs_with_synthetic_data():
    with tempfile.TemporaryDirectory() as tmpdir:
        dict_path = os.path.join(tmpdir, "dictionary.out")
        _write_min_dictionary(dict_path)

        result = run_postprocessing(tmpdir, "testvoc", allow_synthetic=True)

        assert result["status"] == "success"
        assert os.path.exists(os.path.join(tmpdir, "aerosol_data.csv"))
        assert os.path.exists(os.path.join(tmpdir, "postprocessing_summary.json"))
        assert os.path.exists(os.path.join(tmpdir, "plot_audit.json"))
        assert os.path.exists(os.path.join(tmpdir, "data_quality_appendix.json"))
        assert os.path.exists(os.path.join(tmpdir, "data_quality_appendix.pdf"))

        # At least one plot should exist (main plot always generated)
        assert os.path.exists(os.path.join(tmpdir, "plot.png"))


def test_plot_audit_contains_entries():
    with tempfile.TemporaryDirectory() as tmpdir:
        dict_path = os.path.join(tmpdir, "dictionary.out")
        _write_min_dictionary(dict_path)

        result = run_postprocessing(tmpdir, "testvoc", allow_synthetic=True)
        plot_audit = result.get("summary", {}).get("plot_audit", {})

        assert "vbs_plot" in plot_audit
        assert "yield_plot" in plot_audit
        assert "main_plot" in plot_audit
        assert plot_audit["main_plot"]["generated"] is True

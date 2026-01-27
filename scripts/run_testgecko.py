#!/usr/bin/env python3
"""Run GECKO-A/BoxModel validations and build TestGeko dataset."""
import csv
import hashlib
import json
import shutil
import time
import urllib.request
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, TypedDict, Union


class BoxmodelMetrics(TypedDict):
    present: bool
    max_time: Optional[float]
    unique_species: int
    species_last: Dict[str, float]
    row_count: int


class GeckoMetrics(TypedDict, total=False):
    numsp: Optional[int]
    numre: Optional[int]
    num_n: Optional[int]
    numhv: Optional[int]
    numo2: Optional[int]
    nummeo2: Optional[int]
    numain: Optional[int]
    numaou: Optional[int]
    numwin: Optional[int]
    numwou: Optional[int]


RunHashes = Dict[str, str]


class RunMetrics(TypedDict):
    gecko: GeckoMetrics
    boxmodel: BoxmodelMetrics
    hashes: RunHashes


MetricValue = Union[str, float, int, None]


class ComparisonEntry(TypedDict):
    reference: MetricValue
    current: MetricValue
    match: bool


class BoxmodelComparison(TypedDict):
    max_time: ComparisonEntry
    unique_species: ComparisonEntry
    species_last: Dict[str, ComparisonEntry]


class ComparisonResult(TypedDict):
    gecko: Dict[str, ComparisonEntry]
    boxmodel: BoxmodelComparison
    hashes: Dict[str, ComparisonEntry]
    all_hashes_match: bool

VOC_INPUTS = {
    "isoprene": "CH3Cd(=CdH2)CdH=CdH2",
    "alpha-pinene": "C12HCH2CH(C1(CH3)CH3)CH2CdH=Cd2CH3",
    "beta-pinene": "C12HCH2CH(C1(CH3)CH3)CH2CH2Cd2=CdH2",
    "myrcene": "CH3Cd(CH3)=CdHCH2CH2Cd(=CdH2)CdH=CdH2",
    "limonene": "C1H2CH2Cd(CH3)=CdHCH2C1HCd(CH3)=CdH2",
}

TARGET_SPECIES = ["OH", "O3", "NO", "NO2", "HO2"]
BASE_URL = "http://localhost:8000"
RUN_TYPES = ["reference", "current"]
ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "output"
TEST_GEKO = ROOT / "TestGeko"


def http_json(url: str, data: Optional[bytes] = None) -> Dict[str, Any]:
    req = urllib.request.Request(url, data=data)
    if data is not None:
        req.add_header("Content-Type", "application/json")
    with urllib.request.urlopen(req, timeout=60) as resp:
        return json.loads(resp.read().decode("utf-8"))


def submit_job(voc: str) -> Tuple[str, Dict[str, Any]]:
    payload = json.dumps({"voc_name": voc, "job_type": "generator"}).encode("utf-8")
    reply = http_json(f"{BASE_URL}/api/jobs", payload)
    job_id = reply["job_id"]
    while True:
        status = http_json(f"{BASE_URL}/api/jobs/{job_id}")
        if status["status"] in {"completed", "failed"}:
            return job_id, status
        time.sleep(5)


def clean_species(raw: str) -> str:
    text = raw.strip()
    if text.startswith("b'") and text.endswith("'"):
        text = text[2:-1]
    return text.strip().upper()


def parse_boxmodel_metrics(run_dir: Path) -> BoxmodelMetrics:
    csv_path = run_dir / "aerosol_data.csv"
    metrics: BoxmodelMetrics = {
        "present": csv_path.exists(),
        "max_time": None,
        "unique_species": 0,
        "species_last": {},
        "row_count": 0,
    }
    if not csv_path.exists():
        return metrics

    last_values: Dict[str, Tuple[float, float]] = {}
    species_seen: set[str] = set()
    max_time = 0.0
    with csv_path.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            metrics["row_count"] += 1
            try:
                time_val = float(row.get("time", "0"))
                conc = float(row.get("concentration", "0"))
            except ValueError:
                continue
            species = clean_species(row.get("Species", ""))
            species_seen.add(species)
            prev_time = last_values.get(species, (float("-inf"), 0.0))[0]
            if time_val >= prev_time:
                last_values[species] = (time_val, conc)
            if time_val > max_time:
                max_time = time_val

    metrics["max_time"] = max_time
    metrics["unique_species"] = len(species_seen)
    metrics["species_last"] = {
        sp: last_values[sp][1]
        for sp in TARGET_SPECIES
        if sp in last_values
    }
    return metrics


def parse_gecko_metrics(run_dir: Path) -> GeckoMetrics:
    metrics: GeckoMetrics = {}
    out_files = sorted(run_dir.glob("*.out"))
    for out_file in out_files:
        text = out_file.read_text(errors="ignore")
        if "numsp" not in text:
            continue
        for key in [
            "numsp",
            "numre",
            "num_n",
            "numhv",
            "numo2",
            "nummeo2",
            "numain",
            "numaou",
            "numwin",
            "numwou",
        ]:
            metrics[key] = extract_first_int(text, key)
        break
    return metrics


def extract_first_int(text: str, key: str) -> Optional[int]:
    import re

    pattern = re.compile(rf"{key}\s*=\s*(\d+)")
    match = pattern.search(text)
    if match:
        return int(match.group(1))
    return None


def hash_selected_files(run_dir: Path) -> RunHashes:
    important_files = {
        "reactions.dum",
        "dictionary.out",
        "gasspe.dum",
        "partspe.dum",
        "wallspe.dum",
        "listprimary.dat",
        "size.dum",
        "vdiffusion.dat",
        "henry_gromhe.dat",
        "Tg.dat",
        "pvap.nan.dat",
        "aerosol_data.csv",
    }
    hashes: RunHashes = {}
    for path in run_dir.iterdir():
        if path.is_dir():
            continue
        if path.name in important_files:
            hashes[path.name] = sha256_file(path)
    return hashes


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def collect_metrics(run_dir: Path) -> RunMetrics:
    return {
        "gecko": parse_gecko_metrics(run_dir),
        "boxmodel": parse_boxmodel_metrics(run_dir),
        "hashes": hash_selected_files(run_dir),
    }


def ensure_dirs():
    if TEST_GEKO.exists():
        shutil.rmtree(TEST_GEKO)
    for sub in ["reference", "current", "results", "inputs"]:
        (TEST_GEKO / sub).mkdir(parents=True, exist_ok=True)


def save_input_descriptor(voc: str, run_type: str) -> None:
    content: Dict[str, Union[str, float]] = {
        "voc": voc,
        "cheminput": VOC_INPUTS[voc],
        "timestamp": time.time(),
        "run_type": run_type,
    }
    target = TEST_GEKO / "inputs" / f"{voc}_{run_type}.json"
    target.write_text(json.dumps(content, indent=2))


def compare_metrics(ref: RunMetrics, cur: RunMetrics) -> ComparisonResult:
    gecko_comp: Dict[str, ComparisonEntry] = {
        key: {
            "reference": ref["gecko"].get(key),
            "current": cur["gecko"].get(key),
            "match": ref["gecko"].get(key) == cur["gecko"].get(key),
        }
        for key in sorted(set(ref["gecko"]) | set(cur["gecko"]))
    }

    box_comp: BoxmodelComparison = {
        "max_time": numeric_compare(ref["boxmodel"].get("max_time"), cur["boxmodel"].get("max_time")),
        "unique_species": numeric_compare(
            ref["boxmodel"].get("unique_species"),
            cur["boxmodel"].get("unique_species"),
            allow_float=False,
        ),
        "species_last": {},
    }
    ref_box = ref["boxmodel"]
    cur_box = cur["boxmodel"]
    all_species = sorted(set(ref_box.get("species_last", {})) | set(cur_box.get("species_last", {})))
    for sp in all_species:
        box_comp["species_last"][sp] = numeric_compare(
            ref_box.get("species_last", {}).get(sp),
            cur_box.get("species_last", {}).get(sp),
        )

    hash_comp: Dict[str, ComparisonEntry] = {}
    all_files = sorted(set(ref["hashes"]) | set(cur["hashes"]))
    for fname in all_files:
        hash_comp[fname] = {
            "reference": ref["hashes"].get(fname),
            "current": cur["hashes"].get(fname),
            "match": ref["hashes"].get(fname) == cur["hashes"].get(fname),
        }
    return {
        "gecko": gecko_comp,
        "boxmodel": box_comp,
        "hashes": hash_comp,
        "all_hashes_match": all(item["match"] for item in hash_comp.values()),
    }


def numeric_compare(
    ref_val: MetricValue,
    cur_val: MetricValue,
    allow_float: bool = True,
) -> ComparisonEntry:
    if ref_val is None or cur_val is None:
        return {"reference": ref_val, "current": cur_val, "match": ref_val == cur_val}
    if not isinstance(ref_val, (int, float)) or not isinstance(cur_val, (int, float)):
        return {"reference": ref_val, "current": cur_val, "match": ref_val == cur_val}
    if not allow_float:
        return {"reference": ref_val, "current": cur_val, "match": ref_val == cur_val}
    if ref_val == 0:
        match = abs(cur_val) <= 1e-6
    else:
        match = abs((cur_val - ref_val) / ref_val) <= 1e-4
    return {"reference": ref_val, "current": cur_val, "match": match}


def write_report(comparisons: Dict[str, ComparisonResult]) -> None:
    lines = [
        "# TestGeko Validation Report",
        "",
        "| VOC | GECKO Metrics Match | BoxModel Metrics Match | File Hashes Match |",
        "| --- | --- | --- | --- |",
    ]
    for voc in VOC_INPUTS:
        gecko_ok = all(item["match"] for item in comparisons[voc]["gecko"].values())
        box_ok = (
            comparisons[voc]["boxmodel"]["max_time"]["match"]
            and comparisons[voc]["boxmodel"]["unique_species"]["match"]
            and all(val["match"] for val in comparisons[voc]["boxmodel"]["species_last"].values())
        )
        hash_ok = comparisons[voc]["all_hashes_match"]
        lines.append(
            f"| {voc} | {'✅' if gecko_ok else '❌'} | {'✅' if box_ok else '❌'} | {'✅' if hash_ok else '❌'} |"
        )

    (TEST_GEKO / "report.md").write_text("\n".join(lines))


def main():
    ensure_dirs()
    overall_metrics: Dict[str, Dict[str, RunMetrics]] = {voc: {} for voc in VOC_INPUTS}

    for run_type in RUN_TYPES:
        for voc in VOC_INPUTS:
            print(f"Running {run_type} job for {voc}...")
            job_id, status = submit_job(voc)
            if status["status"] != "completed":
                raise RuntimeError(f"Job {job_id} for {voc} failed: {status.get('error')}")

            src = DATA_DIR / job_id
            dest = TEST_GEKO / run_type / voc
            if dest.exists():
                shutil.rmtree(dest)
            shutil.copytree(src, dest)
            save_input_descriptor(voc, run_type)
            stats = collect_metrics(dest)
            overall_metrics[voc][run_type] = stats

    comparisons: Dict[str, ComparisonResult] = {}
    for voc in VOC_INPUTS:
        comparisons[voc] = compare_metrics(
            overall_metrics[voc]["reference"], overall_metrics[voc]["current"]
        )
        result_path = TEST_GEKO / "results" / f"{voc}.json"
        result_path.write_text(json.dumps(comparisons[voc], indent=2))

    summary_path = TEST_GEKO / "results" / "summary.json"
    summary_path.write_text(json.dumps(comparisons, indent=2))
    write_report(comparisons)
    print("Validation complete. See TestGeko directory for details.")


if __name__ == "__main__":
    main()

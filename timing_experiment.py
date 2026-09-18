"""Section 3.3 acceptance and timing experiment, using ordinary Python loops.

Run a small pilot: python3 timing_experiment.py --pilot
Run the original sample counts: python3 timing_experiment.py
Only the Python standard library and the adjacent sampler module are needed.
The existing samplers choose the proposal orientation by time reversal.
"""

import argparse
import csv
from datetime import datetime, timezone
import gc
import hashlib
import json
from pathlib import Path
import platform
import random
import statistics
import sys
import time

if __package__:
    from . import ICB_gh_samplers as samplers
else:
    import ICB_gh_samplers as samplers


# Experiment settings for Section 3.3 of the manuscript.
A = 2.0
T = 3.0
T_VALUES = [0.1, 1.0, 1.5, 2.0, 2.5, 2.9]
ENDPOINT_MARGIN = 0.05
GRID_POINTS = 5
B_ACCEPT = 2000
B_TIME = 10000
TIMING_REPS = 1
PILOT_B_ACCEPT = 100
PILOT_B_TIME = 500
ACCEPT_SEED = 100000
TIMING_SEED = 200000


def build_scenarios():
    grid = []
    for i in range(GRID_POINTS):
        grid.append(ENDPOINT_MARGIN + (A - 2.0 * ENDPOINT_MARGIN) * i / (GRID_POINTS - 1))
    scenarios = []
    for algorithm in (3, 4, 5, 6):
        x_values = grid if algorithm in (3, 5) else [0.0]
        z_values = grid if algorithm in (3, 4) else [A]
        scenario_id = 0
        # Time varies fastest, then the initial endpoint, then the final endpoint.
        for z in z_values:
            for x in x_values:
                for t in T_VALUES:
                    scenario_id += 1
                    scenarios.append(dict(algorithm=algorithm, scenario_id=scenario_id,
                                          t=t, T=T, x=x, z=z, a=A))
    return scenarios


def select_sampler(scenario):
    algorithm = scenario["algorithm"]
    t, horizon, x, z, a = (scenario[key] for key in ("t", "T", "x", "z", "a"))
    if algorithm == 3:
        return samplers.ICB_interpolate_BB_I, (t, horizon, x, z, a)
    if algorithm == 4:
        return samplers.ICB_interpolate_BB_II, (t, horizon, z, a)
    if algorithm == 5:
        return samplers.ICB_interpolate_BB_III, (t, horizon, x, a)
    if algorithm == 6:
        return samplers.ICB_interpolate_BB_IV, (t, horizon, a)
    raise ValueError("Algorithm must be 3, 4, 5 or 6")


def count_proposals(draw, args, n_draws, a):
    """Count outer bridge proposals, excluding internal g/h rejections.

    Count BOTH public proposal samplers because time reversal can switch
    between g and h. These wrappers are used only in this serial counting
    pass; restore the original functions before the separate timing pass.
    """
    original_g = samplers.ICB_sample_g
    original_h = samplers.ICB_sample_h
    counts = [0, 0]

    def counted_g(t, x, width):
        counts[0] += 1
        return original_g(t, x, width)

    def counted_h(t, width):
        counts[1] += 1
        return original_h(t, width)

    proposals_per_draw = []
    samplers.ICB_sample_g = counted_g
    samplers.ICB_sample_h = counted_h
    try:
        for i in range(n_draws):
            before = counts[0] + counts[1]
            value = draw(*args)
            if not 0.0 < value < a:
                raise ArithmeticError("A bridge draw is outside the interval")
            proposals_per_draw.append(counts[0] + counts[1] - before)
    finally:
        samplers.ICB_sample_g = original_g
        samplers.ICB_sample_h = original_h
    return proposals_per_draw, counts[0], counts[1]


def time_draws(draw, args, n_draws, repetitions, scenario_id, a):
    elapsed = []
    for repetition in range(1, repetitions + 1):
        random.seed(TIMING_SEED + 100 * scenario_id + repetition)
        gc.collect()
        values = [0.0] * n_draws
        started = time.perf_counter()
        for i in range(n_draws):
            values[i] = draw(*args)
        elapsed.append(time.perf_counter() - started)
        # Validation, seed setup, allocation and file output are not timed.
        for value in values:
            if not 0.0 < value < a:
                raise ArithmeticError("A timed bridge draw is outside the interval")
    return elapsed


def summarize_results(rows):
    summary = []
    for algorithm in (3, 4, 5, 6):
        for t in T_VALUES:
            acceptance_rates = []
            scenario_times = []
            for row in rows:
                if row["algorithm"] == algorithm and row["t"] == t:
                    acceptance_rates.append(row["acceptance_rate"])
                    scenario_times.append(row["elapsed_median_sec"])
            if acceptance_rates:
                summary.append(dict(
                    algorithm=algorithm, t=t, scenario_count=len(acceptance_rates),
                    acceptance_rate_mean=statistics.mean(acceptance_rates),
                    elapsed_median_sec=statistics.median(scenario_times)))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", action="store_true",
                        help="Use 100 acceptance draws and 500 timing draws per scenario")
    options = parser.parse_args()
    mode = "pilot" if options.pilot else "full"
    b_accept = PILOT_B_ACCEPT if options.pilot else B_ACCEPT
    b_time = PILOT_B_TIME if options.pilot else B_TIME
    scenarios = build_scenarios()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_%f")
    output = Path(__file__).resolve().parent / "timing_results" / f"{mode}_{stamp}"
    output.mkdir(parents=True)
    print(f"{mode.upper()}: {len(scenarios)} scenarios; {b_accept} acceptance draws, "
          f"{b_time} timed draws, {TIMING_REPS} timing repetition(s) each.", flush=True)
    print(f"Output: {output}", flush=True)

    metadata = dict(
        status="running", mode=mode, started_utc=datetime.now(timezone.utc).isoformat(),
        a=A, T=T, t_values=T_VALUES, endpoint_margin=ENDPOINT_MARGIN,
        endpoint_grid_points=GRID_POINTS, scenarios=len(scenarios),
        B_accept=b_accept, B_time=b_time, timing_reps=TIMING_REPS,
        automatic_time_reversal=True, python=sys.version, executable=sys.executable,
        platform=platform.platform(), machine=platform.machine(), processor=platform.processor(),
        timer="time.perf_counter", timing_execution="serial",
        acceptance_seed_rule=f"{ACCEPT_SEED} + scenario_id",
        timing_seed_rule=f"{TIMING_SEED} + 100*scenario_id + repetition (starting at 1)",
        scenario_id_rule="Restart for each algorithm; time varies fastest, then x, then z",
        acceptance_definition="Accepted draws / outer g-or-h proposals within each scenario",
        acceptance_summary="Arithmetic mean of per-scenario acceptance rates",
        timing_summary="Median across endpoint scenarios of within-scenario median elapsed seconds",
        counted_functions=["ICB_sample_g", "ICB_sample_h"],
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        sampler_sha256=hashlib.sha256(Path(samplers.__file__).read_bytes()).hexdigest())
    metadata_path = output / "metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")

    fields = ["algorithm", "scenario_id", "t", "T", "x", "z", "a",
              "B_accept", "g_proposals", "h_proposals", "total_proposals",
              "acceptance_rate", "mean_proposals_per_draw", "median_proposals_per_draw",
              "B_time", "timing_reps", "elapsed_median_sec", "elapsed_mean_sec",
              "elapsed_sd_sec", "sec_per_draw_median", "acceptance_seed"]
    timing_fields = ["algorithm", "scenario_id", "t", "x", "z", "repetition",
                     "B_time", "timing_seed", "elapsed_sec"]
    rows = []
    experiment_start = time.perf_counter()
    with (output / "scenarios.csv").open("w", newline="") as results_file, \
            (output / "timing_batches.csv").open("w", newline="") as timings_file:
        writer = csv.DictWriter(results_file, fieldnames=fields)
        timing_writer = csv.DictWriter(timings_file, fieldnames=timing_fields)
        writer.writeheader()
        timing_writer.writeheader()
        for number, scenario in enumerate(scenarios, 1):
            algorithm = scenario["algorithm"]
            scenario_id = scenario["scenario_id"]
            draw, args = select_sampler(scenario)
            random.seed(ACCEPT_SEED + scenario_id)
            counts, g_count, h_count = count_proposals(draw, args, b_accept, scenario["a"])
            total_proposals = g_count + h_count
            if min(counts) < 1 or sum(counts) != total_proposals:
                raise ArithmeticError("Invalid outer proposal counts")
            # count_proposals has restored both original functions at this point.
            elapsed = time_draws(draw, args, b_time, TIMING_REPS, scenario_id, scenario["a"])
            median_elapsed = statistics.median(elapsed)
            row = dict(scenario)
            row.update(B_accept=b_accept, g_proposals=g_count, h_proposals=h_count,
                       total_proposals=total_proposals, acceptance_rate=b_accept / total_proposals,
                       mean_proposals_per_draw=statistics.mean(counts),
                       median_proposals_per_draw=statistics.median(counts),
                       B_time=b_time, timing_reps=TIMING_REPS,
                       elapsed_median_sec=median_elapsed, elapsed_mean_sec=statistics.mean(elapsed),
                       elapsed_sd_sec=statistics.stdev(elapsed) if len(elapsed) > 1 else "",
                       sec_per_draw_median=median_elapsed / b_time,
                       acceptance_seed=ACCEPT_SEED + scenario_id)
            rows.append(row)
            writer.writerow(row)
            results_file.flush()
            for repetition, seconds in enumerate(elapsed, 1):
                timing_writer.writerow(dict(
                    algorithm=algorithm, scenario_id=scenario_id, t=scenario["t"],
                    x=scenario["x"], z=scenario["z"], repetition=repetition,
                    B_time=b_time, timing_seed=TIMING_SEED + 100 * scenario_id + repetition,
                    elapsed_sec=seconds))
            timings_file.flush()
            print(f"[{number}/{len(scenarios)}] Algorithm {algorithm}: "
                  f"t={scenario['t']:g}, x={scenario['x']:g}, z={scenario['z']:g}; "
                  f"acceptance={row['acceptance_rate']:.4f}, {b_time} draws={median_elapsed:.3f} s",
                  flush=True)

    summary = summarize_results(rows)
    with (output / "summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["algorithm", "t", "scenario_count",
                                                   "acceptance_rate_mean", "elapsed_median_sec",
                                                   "B_accept", "B_time", "timing_reps", "mode"])
        writer.writeheader()
        for row in summary:
            writer.writerow(dict(row, B_accept=b_accept, B_time=b_time,
                                 timing_reps=TIMING_REPS, mode=mode))

    # A compact table for inspection. Pilot times refer to 500 draws, not 10,000.
    table = [f"# Section 3.3 timing experiment: {mode}", "",
             f"Each acceptance estimate uses {b_accept} accepted draws. "
             f"Elapsed times are seconds for {b_time} draws, with {TIMING_REPS} repetition(s).",
             "", "| t | Alg. 3 acc. | Alg. 3 seconds | Alg. 4 acc. | Alg. 4 seconds | "
             "Alg. 5 acc. | Alg. 5 seconds | Alg. 6 acc. | Alg. 6 seconds |",
             "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for t in T_VALUES:
        cells = [f"{t:g}"]
        for algorithm in (3, 4, 5, 6):
            for row in summary:
                if row["algorithm"] == algorithm and row["t"] == t:
                    cells.append(f"{row['acceptance_rate_mean']:.4f}")
                    cells.append(f"{row['elapsed_median_sec']:.4f}")
        table.append("| " + " | ".join(cells) + " |")
    table.extend(["", "Acceptance is averaged over endpoint settings; timings are medians over "
                  "those settings. Proposal counting is disabled during timing. Automatic "
                  "time reversal is enabled."])
    if options.pilot:
        table.extend(["", "Pilot measurements are preliminary."])
    (output / "summary.md").write_text("\n".join(table) + "\n")
    metadata.update(status="complete", completed_utc=datetime.now(timezone.utc).isoformat(),
                    completed_scenarios=len(rows), total_elapsed_sec=time.perf_counter()-experiment_start)
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print("\n" + "\n".join(table), flush=True)
    print(f"Saved results in {output}", flush=True)


if __name__ == "__main__":
    main()
